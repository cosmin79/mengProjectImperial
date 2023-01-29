#include "IOUtil.hpp"
#include "OptParse.hpp"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <queue>
#include <memory>
#include <chrono>

const string EPS_OPTION = "-EPS";
const string MIN_PTS_OPTION = "-MIN_PTS";
const string FILE_OPTION = "-FILE";
const string INDEX_TYPE = "-INDEX";
const string VERBOSE_OPTION = "-VERBOSE";

enum IndexTypes {
    LINEAR_ARRAY, KD_TREES, NO_OPTIONS
};

typedef Index<L2<double>> fIndex;

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

class DBScan {
private:
    string fileName;
    double eps;
    int minPts;
    bool verbose;
    shared_ptr<fIndex> flannIndex;
    shared_ptr<DSU> dsu;

    vector< pt<double> > dataset;
    Matrix<double> flannDataset;

    // Result
    vector<Categories> categories;

public:
    DBScan(double _eps, int _minPts, string& _fileName, int indexType, bool _verbose) : 
        fileName(_fileName), eps(_eps), minPts(_minPts), verbose(_verbose) {
        // I / O
        parseFromFile(dataset, fileName, false);
        getFlannMatrix(flannDataset, dataset);
        switch (indexType) {
            case LINEAR_ARRAY:
                flannIndex = shared_ptr<fIndex> (new fIndex(flannDataset, flann::LinearIndexParams()));
                break ;
            case KD_TREES:
                flannIndex = shared_ptr<fIndex> (new fIndex(flannDataset, flann::KDTreeIndexParams(4)));
                break ;
        }
        dsu = shared_ptr<DSU>(new DSU(dataset.size()));

        categories.assign(dataset.size(), NOISE);
    }

    void runDBScan();

private:
    vector<int> regionQuery(pt<double>&, double);
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

void DBScan::runDBScan() {
    flannIndex -> buildIndex();

    for (auto& p: dataset) {
        vector<int> ngbIndices = regionQuery(p, eps);
        if (ngbIndices.size() >= minPts) {
            categories[p.idx] = CORE;
            for (auto& ngb: ngbIndices) {
                // doesn't belong to a cluster yet
                if (categories[ngb] == NOISE) {
                    // It can still become CORE later on
                    categories[ngb] = OTHER;
                    dsu -> unite(ngb, p.idx);
                } else if (categories[ngb] == CORE) {
                    dsu -> unite(ngb, p.idx);
                }
            }
        }
    }

    int noClusters = 0;
    for (size_t i = 0; i < dataset.size(); i++) {
        // i.e. "i" is a global root
        if (categories[i] != NOISE && dsu -> find(i) == i) {
            noClusters++;
        }
    }

    printf("(dbScanDSU) Number clusters: %d\n", noClusters);
}

int main(int argc, char** argv) {
    chrono::system_clock::time_point startTime, endTime;
    startTime = chrono::system_clock::now();

    OptionParser optionParser(argc, argv);
    double eps = stod(optionParser.getOptionValue(EPS_OPTION));
    int minPts = stoi(optionParser.getOptionValue(MIN_PTS_OPTION));
    string file = optionParser.getOptionValue(FILE_OPTION);
    int verbose = stoi(optionParser.getOptionValue(VERBOSE_OPTION));

    int indexType = stoi(optionParser.getOptionValue(INDEX_TYPE));
    int noOptions = (int)(IndexTypes::NO_OPTIONS);
    if (indexType < 0 || indexType > noOptions) {
        throw invalid_argument("Expected an index type between 0 and " + to_string(noOptions));
    }

    DBScan clusteringAlgo(eps, minPts, file, indexType, verbose);
    clusteringAlgo.runDBScan();
    
    endTime = std::chrono::system_clock::now();
    printf("(dbScanDSU) Total running time in milliseconds: %lld\n",
           chrono::duration_cast<chrono::milliseconds>(endTime - startTime).count());
    return 0;
}
