#include "IOUtil.hpp"
#include "OptParse.hpp"
#include <cstdio>
#include <vector>
#include <algorithm>
#include <string>
#include <cassert>
#include <queue>
#include <memory>

const string EPS_OPTION = "-EPS";
const string MIN_PTS_OPTION = "-MIN_PTS";
const string FILE_OPTION = "-FILE";
const string INDEX_TYPE = "-INDEX";
enum IndexTypes {
    LINEAR_ARRAY, KD_TREES, NO_OPTIONS
};

const int UNDEF_CLUSTER = -1;
typedef Index<L2<double>> fIndex;

enum Categories {
    CORE, NOISE, OTHER
};

class DBScan {
private:
    string fileName;
    double eps;
    int minPts;
    shared_ptr<fIndex> flannIndex;

    vector< pt<double> > dataset;
    Matrix<double> flannDataset;

    vector<bool> visited;
    vector<bool> inQueue;

    // Result
    vector<int> whichCluster;
    vector<Categories> categories;

public:
    DBScan(double _eps, int _minPts, string& _fileName, int indexType) : fileName(_fileName), eps(_eps), minPts(_minPts) {
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

        visited.assign(dataset.size(), false);
        inQueue.assign(dataset.size(), false);
        whichCluster.assign(dataset.size(), UNDEF_CLUSTER);
        categories.assign(dataset.size(), NOISE);
    }

    void runDBScan();

private:
    vector<int> regionQuery(pt<double>&, double);

    inline void setClusterIfOrphan(int, int);

    void expandCluster(pt<double>&, vector<int>&, int);
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

inline void DBScan::setClusterIfOrphan(int ptIdx, int currCluster) {
    if (whichCluster[ptIdx] == UNDEF_CLUSTER) {
        whichCluster[ptIdx] = currCluster;
        // this point is on the fringe
        categories[ptIdx] = OTHER;
    }
}

void DBScan::expandCluster(pt<double>& p, vector<int>& ngbIndices, int currCluster) {
    // add curr point
    whichCluster[p.idx] = currCluster;
    categories[p.idx] = CORE;

    // setup queue
    queue<int> Q;
    for (auto& neighbourIdx: ngbIndices) {
        if (!visited[neighbourIdx]) {
            Q.push(neighbourIdx);
            inQueue[neighbourIdx] = true;
        } else {
            setClusterIfOrphan(neighbourIdx, currCluster);
        }
    }

    while (!Q.empty()) {
        int neighbourIdx = Q.front();
        visited[neighbourIdx] = true;
        Q.pop();
        inQueue[neighbourIdx] = false;

        vector<int> ngbNgbsIndices = regionQuery(dataset[neighbourIdx], eps);
        // expand ?
        if (ngbNgbsIndices.size() >= minPts) {
            categories[neighbourIdx] = CORE;
            for (auto& candIdx: ngbNgbsIndices) {
                if (!visited[candIdx]) {
                    // is it already in the neighbours set?
                    if (!inQueue[candIdx]) {
                        Q.push(candIdx);
                        inQueue[candIdx] = true;
                    }
                } else {
                    setClusterIfOrphan(candIdx, currCluster);
                }
            }
        }
    }
}

void DBScan:: runDBScan() {
    flannIndex -> buildIndex();

    int noClusters = 0;
    for (auto& p: dataset) {
        if (visited[p.idx]) {
            continue ;
        }
        visited[p.idx] = true;

        vector<int> ngbIndices = regionQuery(p, eps);
        if (ngbIndices.size() >= minPts) {
            expandCluster(p, ngbIndices, ++noClusters);
        }
    }

    printf("NO. clusters: %d\n", noClusters);
}

int main(int argc, char** argv) {
    OptionParser optionParser(argc, argv);
    double eps = stod(optionParser.getOptionValue(EPS_OPTION));
    int minPts = stoi(optionParser.getOptionValue(MIN_PTS_OPTION));
    string file = optionParser.getOptionValue(FILE_OPTION);

    int indexType = stoi(optionParser.getOptionValue(INDEX_TYPE));
    int noOptions = (int)(IndexTypes::NO_OPTIONS);
    if (indexType < 0 || indexType > noOptions) {
        throw invalid_argument("Expected an index type between 0 and " + to_string(noOptions));
    }

    DBScan clusteringAlgo(eps, minPts, file, indexType);
    clusteringAlgo.runDBScan();
    return 0;
}
