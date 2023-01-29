#include "IOUtil.hpp"
#include "OptParse.hpp"
#include "mpi.h"
#include <cstdio>
#include <iostream>
#include <chrono>
#include <string>
#include <vector>
#include <algorithm>
#include <string>
#include <memory>
#include <stddef.h>
using namespace std;

const string EPS_OPTION = "-EPS";
const string MIN_PTS_OPTION = "-MIN_PTS";
const string ORIG_FILE_OPTION = "-ORIG_FILE";
const string INPUT_DIR_OPTION = "-INPUT_DIR";
const string VERBOSE_OPTION = "-VERBOSE";
const string OPTIMIZE_COMM_OPTION = "-OPTIMIZE_COMM";

typedef Index<L2<double>> fIndex;
typedef chrono::system_clock::time_point tp;

int _machines;
int _rank;

inline long long millisecondsElapsed(tp&);

class DSU {
private:
    int N;
    vector<int> par, deg;

public:
    DSU(int _N) : N(_N) {
        par.resize(N);
        deg.assign(N, 0);
        for (int i = 0; i < N; i++) {
            MakeSet(i);
        }
    }

    void unite(int x, int y) {
        x = find(x); y = find(y);
        if (deg[x] < deg[y]) {
            par[x] = y;
        } else {
            par[y] = x;
            if (deg[y] == deg[x]) {
                deg[x]++;
            }
        }
    }

    int find(int x) {
        int root = x, aux;
        while (par[root] != root) {
            root = par[root];
        }

        while (par[x] != x) {
            aux = par[x];
            par[x] = root;
            x = aux;
        }

        return root;
    }

private:
    inline void MakeSet(int node) {
        par[node] = node;
    }
};

enum Categories {
    CORE, NOISE, OTHER
};

struct foreignEdge {
    // Note that from is a core point
    int from, to;
};

struct localNodeInfo {
    // categ is 0, 1, 2 i.e. an instance of type Categories...but I am more confident in
    // encoding this as a normal int
    int node, par, categ;
};

class DBScan {
private:
    string fileName;
    double eps;
    int minPts;
    bool verbose, optimizeCommunication;
    tp startTime;
    shared_ptr<fIndex> flannIndex;
    shared_ptr<DSU> dsu;

    vector< pt<double> > dataset;
    Matrix<double> flannDataset;

    // Result
    vector<Categories> categories;

// This info needs to be visible to the outside so that it can be sent to master
public:
    vector<localNodeInfo> localInfo;
    vector<foreignEdge> foreignEdges;

public:
    DBScan(double _eps, int _minPts, string& _fileName, bool _verbose, bool _optimizeCommunication, tp& _startTime) :
        fileName(_fileName), eps(_eps), minPts(_minPts), verbose(_verbose),
        optimizeCommunication(_optimizeCommunication), startTime(_startTime) {
        // I / O ; note "true" means it is a partition file
        parseFromFile(dataset, fileName, true);
        for (size_t i = 0; i < dataset.size(); i++) {
            dataset[i].setLocalIndex(i);
        }
        // Remember ; those point have different global indexing ; flann uses a local indexing from 0
        getFlannMatrix(flannDataset, dataset);
        if (_rank == 0) {
            printf("(dbScanMulti) Time after local read in milliseconds: %lld\n", millisecondsElapsed(startTime));
        }
        flannIndex = shared_ptr<fIndex> (new fIndex(flannDataset, flann::KDTreeIndexParams(4)));
        dsu = shared_ptr<DSU>(new DSU(dataset.size()));

        categories.assign(dataset.size(), NOISE);
    }

    void runDBScan();

private:
    vector<int> regionQuery(pt<double>&, double);

    inline bool belongsToThisNode(pt<double>&);
};


vector<int> DBScan::regionQuery(pt<double>& p, double eps) {
    int cols = p.coords.size();

    Matrix<double> flannSingleQuery = Matrix<double>(new double[cols], 1, cols);
    copy(p.coords.begin(), p.coords.end(), flannSingleQuery[0]);

    vector<vector<int>> indices(1, vector<int>());
    vector<vector<double>> dists(1, vector<double>());

    flannIndex -> radiusSearch(flannSingleQuery, indices, dists, eps, flann::SearchParams(128));

    return indices[0];
}

inline bool DBScan::belongsToThisNode(pt<double>& p) {
    return p.partition == _rank;
}

void DBScan::runDBScan() {
    flannIndex -> buildIndex();
    if (_rank == 0) {
        printf("(dbScanMulti) Time after local index creation in milliseconds: %lld\n", millisecondsElapsed(startTime));
    }

    // crossEdges[dest] is a vector containing indices of core points from this partition
    // crossEdges[localIdx], where dataset[localIdx] belongs to this node, should be empty
    vector<vector<int>> crossEdges(dataset.size(), vector<int>());
    for (auto& p: dataset) {
        if (belongsToThisNode(p)) {
            vector<int> ngbIndices = regionQuery(p, eps);
            if (ngbIndices.size() >= minPts) {
                categories[p.localIdx] = CORE;
                for (auto& ngb: ngbIndices) {
                    if (!belongsToThisNode(dataset[ngb])) {
                        crossEdges[ngb].push_back(p.localIdx);
                    } else {
                        // is not part of a cluster
                        if (categories[ngb] == NOISE) {
                            // It can still become CORE later on
                            categories[ngb] = OTHER;
                            dsu -> unite(ngb, p.localIdx);
                        } else if (categories[ngb] == CORE) {
                            dsu -> unite(ngb, p.localIdx);
                        }
                    }   
                }
            }
        }
    }

    // Make sure to use global indices i.e. ngb is a local index
    if (optimizeCommunication) {
        vector<bool> markedNgbs(dataset.size(), false);
        for (size_t ghostLocalIdx = 0; ghostLocalIdx < dataset.size(); ghostLocalIdx++) {
            for (auto& coreIdx: crossEdges[ghostLocalIdx]) {
                int rootCoreidx = dsu -> find(coreIdx);
                if (!markedNgbs[rootCoreidx]) {
                    foreignEdges.push_back({dataset[rootCoreidx].idx, dataset[ghostLocalIdx].idx});
                    markedNgbs[rootCoreidx] = true;
                }
            }

            // erase visited for this point
            for (auto& coreIdx: crossEdges[ghostLocalIdx]) {
                markedNgbs[dsu -> find(coreIdx)] = false;
            }
        }
    } else {
        for (size_t ghostLocalIdx = 0; ghostLocalIdx < dataset.size(); ghostLocalIdx++) {
            for (auto& coreIdx: crossEdges[ghostLocalIdx]) {
                foreignEdges.push_back({dataset[coreIdx].idx, dataset[ghostLocalIdx].idx});
            }
        }
    }

    int localClusters = 0;
    for (size_t i = 0; i < dataset.size(); i++) {
        // i.e. "i" is a global root
        if (belongsToThisNode(dataset[i])) {
            if (categories[i] != NOISE && dsu -> find(i) == i) {
                localClusters++;
            }
            // Remember: global indices!
            localInfo.push_back({dataset[i].idx, 
                dataset[dsu -> find(i)].idx, 
                (int)categories[i]});
        }
    }

    if (verbose) {
        printf("NODE %d -> LOCAL clusters: %d\n", _rank, localClusters);
    }
}

void declareForeignEdge(MPI_Datatype& mpi_foreignEdge) {
    const int nitems = 2;
    int blocklengths[2] = {1, 1};
    MPI_Datatype types[2] = {MPI_INT, MPI_INT};
    MPI_Aint offsets[2];
    offsets[0] = offsetof(foreignEdge, from);
    offsets[1] = offsetof(foreignEdge, to);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_foreignEdge);
    MPI_Type_commit(&mpi_foreignEdge);
}

void declareLocalNodeInfo(MPI_Datatype& mpi_localNodeInfo) {
    const int nitems = 3;
    int blocklengths[3] = {1, 1, 1};
    MPI_Datatype types[3] = {MPI_INT, MPI_INT, MPI_INT};
    MPI_Aint offsets[3];
    offsets[0] = offsetof(localNodeInfo, node);
    offsets[1] = offsetof(localNodeInfo, par);
    offsets[2] = offsetof(localNodeInfo, categ);
    MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_localNodeInfo);
    MPI_Type_commit(&mpi_localNodeInfo);
}

void gatherLocalInfo(vector<localNodeInfo>& allLocalInfo, MPI_Datatype& mpi_localNodeInfo, DBScan& localSolver) {
    vector<int> countsLN(_machines, 0);
    int localCnt = (int)localSolver.localInfo.size();
    MPI_Gather(&localCnt, 1, MPI_INT, countsLN.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> dispsLN(_machines, 0);
    for (int i = 1; i < _machines; i++) {
        dispsLN[i] = dispsLN[i - 1] + countsLN[i - 1];
    }
    int totalCntLN = dispsLN[_machines - 1] + countsLN[_machines - 1];

    allLocalInfo.resize(totalCntLN);
    MPI_Gatherv(localSolver.localInfo.data(), localCnt, mpi_localNodeInfo,
            allLocalInfo.data(), countsLN.data(), dispsLN.data(), mpi_localNodeInfo, 
            0, MPI_COMM_WORLD);
}

// This function and the one above look similar. They can be templatized...I think
void gatherForeignEdges(vector<foreignEdge>& foreignEdges, MPI_Datatype& mpi_foreignEdge, DBScan& localSolver) {
    vector<int> countsFE(_machines, 0);
    int localFE = (int)localSolver.foreignEdges.size();
    MPI_Gather(&localFE, 1, MPI_INT, countsFE.data(), 1, MPI_INT, 0, MPI_COMM_WORLD);

    vector<int> dispsFE(_machines, 0);
    for (int i = 1; i < _machines; i++) {
        dispsFE[i] = dispsFE[i - 1] + countsFE[i - 1];
    }
    int totalCntFE = dispsFE[_machines - 1] + countsFE[_machines - 1];

    foreignEdges.resize(totalCntFE);
    MPI_Gatherv(localSolver.foreignEdges.data(), localFE, mpi_foreignEdge,
        foreignEdges.data(), countsFE.data(), dispsFE.data(), mpi_foreignEdge,
        0, MPI_COMM_WORLD);
}

void aggregateSolution(vector<localNodeInfo>& allLocalInfo, vector<foreignEdge>& foreignEdges, bool verbose) {
    // this is the total no. of points
    int N = (int)allLocalInfo.size();
    DSU dsu(N);
    vector<Categories> categories(N, NOISE);

    for (auto& nodeInfo: allLocalInfo) {
        categories[nodeInfo.node] = static_cast<Categories>(nodeInfo.categ);
        dsu.unite(nodeInfo.node, nodeInfo.par);
    }

    for (auto& edge: foreignEdges) {
        // we know edge.from is a core point ; if edge.to is NOISE it's not part of a cluster yet
        if (categories[edge.to] == CORE) {
            dsu.unite(edge.from, edge.to);
        } else if (categories[edge.to] == NOISE) {
            categories[edge.to] = OTHER;
            dsu.unite(edge.from, edge.to);
        }
    }

    int noClusters = 0;
    for (int i = 0; i < N; i++) {
        if (categories[i] != NOISE && dsu.find(i) == i) {
            noClusters++;
        }
    }

    if (verbose) {
        printf("GOT %d clusters after merging everything\n", noClusters);
    }
}

inline long long millisecondsElapsed(tp& startTime) {
    tp currTime = std::chrono::system_clock::now();

    return chrono::duration_cast<chrono::milliseconds>(currTime - startTime).count();
}

int main(int argc, char **argv) {
    tp startTime;

    // Parse options
    OptionParser optionParser(argc, argv);
    double eps = stod(optionParser.getOptionValue(EPS_OPTION));
    int minPts = stoi(optionParser.getOptionValue(MIN_PTS_OPTION));
    string intputDir = optionParser.getOptionValue(INPUT_DIR_OPTION);
    string originalFile = optionParser.getOptionValue(ORIG_FILE_OPTION);
    int verbose = stoi(optionParser.getOptionValue(VERBOSE_OPTION));
    int optimizeCommunication = stoi(optionParser.getOptionValue(OPTIMIZE_COMM_OPTION));

    if (MPI_Init(&argc, &argv) != MPI_SUCCESS) {
        printf("Could not initialise MPI!\n");
        return 1;
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &_rank);
    MPI_Comm_size(MPI_COMM_WORLD, &_machines);

    startTime = chrono::system_clock::now();
    // The file this machine should read from can be obtained using this encoding
    string file = intputDir + "/" + originalFile + to_string(_rank);
    DBScan localSolver(eps, minPts, file, verbose, optimizeCommunication, startTime);
    localSolver.runDBScan();
    if (_rank == 0) {
        printf("(dbScanMulti) Time after local step in milliseconds: %lld\n", millisecondsElapsed(startTime));
    }

    /* create a type for struct foreignEdge */
    MPI_Datatype mpi_foreignEdge;
    declareForeignEdge(mpi_foreignEdge);

    /* create a type for struct localNodeInfo */
    MPI_Datatype mpi_localNodeInfo;
    declareLocalNodeInfo(mpi_localNodeInfo);

    // Now the communication part
    /* int MPI_Gather(void *sendbuf, int sendcount, MPI_Datatype senddatatype,
        void *recvbuf, int recvcount, MPI_Datatype recvdatatype,
        int target, MPI_Comm comm)

       int MPI_Gatherv(void *sendbuf, int sendcount, MPI_Datatype senddatatype,
        void *recvbuf, int *recvcounts, int *displs, MPI_Datatype recvdatatype,
        int target, MPI_Comm comm)*/

    
    vector<localNodeInfo> allLocalInfo;
    gatherLocalInfo(allLocalInfo, mpi_localNodeInfo, localSolver);
    if (_rank == 0) {
        printf("(dbScanMulti) Time after get local clusters info in milliseconds: %lld\n", millisecondsElapsed(startTime));
    }

    vector<foreignEdge> allForeignEdges;
    gatherForeignEdges(allForeignEdges, mpi_foreignEdge, localSolver);
    if (_rank == 0) {
        printf("(dbScanMulti) Time after get foreign edges info in milliseconds: %lld\n", millisecondsElapsed(startTime));
    }

    if (_rank == 0) {
        if (verbose) {
            printf("LOCAL NODE INFO %d\n", (int)allLocalInfo.size());
        }

        printf("(dbScanMulti) FOREIGN EDGES NO: %d\n", (int)allForeignEdges.size());
        aggregateSolution(allLocalInfo, allForeignEdges, verbose);

        printf("(dbScanMulti) Parallel runtime in milliseconds: %lld\n", millisecondsElapsed(startTime));
    }

    if (MPI_Finalize() != MPI_SUCCESS) {
        printf("Could not finalize MPI properly!\n");
        return 1;
    }

    return 0;
}
