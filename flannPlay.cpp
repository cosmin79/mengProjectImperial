
#include <flann/flann.hpp>
#include <flann/io/hdf5.h>

#include <stdio.h>
#include <vector>
#include <algorithm>
#include <string>
#include <sstream>
#include <fstream>
#include <cassert>
#include <unordered_set>

using namespace flann;
using namespace std;

const string INPUT_FILE = "input/2DSimple.pts";

template <typename T>
struct pt {
    int idx;
    vector<T> coords;

    pt(int _idx, vector<T>& _coords) {
        idx = _idx;
        coords = _coords;
    }
};

inline void parseDimensions(string& s, int& rows, int& cols) {
    stringstream ss(s);
    ss >> rows >> cols;
}

template <typename T>
inline pt<T> parsePoint(string& line) {
    stringstream ss(line);
    int idx;
    T elem;
    vector<T> coords;
    
    ss >> idx;
    while (ss >> elem) {
        coords.push_back(elem);
    }

    return pt<T>(idx - 1, coords);
}

template <typename T>
void parseFromFile(vector< pt<T> >& dataset, const string& name) {
    ifstream in(name);
    string line;
    int noPoints, noDims;

    getline(in, line);
    parseDimensions(line, noPoints, noDims);
    dataset.reserve(noPoints);

    while (getline(in, line)) {
        dataset.push_back(parsePoint<T>(line));
        assert(dataset.back().coords.size() == noDims);
    }
    assert(dataset.size() == noPoints);
}

template<typename T>
void getFlannMatrix(Matrix<T>& flannDataset, vector< pt<T> >& dataset) {
    assert(!dataset.empty());
    int rows = dataset.size(), cols = dataset[0].coords.size();
    flannDataset = Matrix<T>(new T[rows * cols], rows, cols);

    for (size_t rowIdx = 0; rowIdx < dataset.size(); rowIdx++) {
        for (size_t colIdx = 0; colIdx < dataset[rowIdx].coords.size(); colIdx++) {
            flannDataset[rowIdx][colIdx] = dataset[rowIdx].coords[colIdx];
        }
    }
}

int main(int argc, char** argv)
{
    int nn = 3;
    float eps = 50.0;

    vector< pt<float> > datasetC;
    parseFromFile(datasetC, INPUT_FILE);

    Matrix<float> dataset;
    Matrix<float> query;
    getFlannMatrix(dataset, datasetC);
    query = Matrix<float>(new float[cols], 1, cols);

    vector<vector<int>> indices(query.rows, vector<int>());
    vector<vector<float>> dists(query.rows, vector<float>());
    //Matrix<int> indices(new int[query.rows*nn], query.rows, nn);
    //Matrix<float> dists(new float[query.rows*nn], query.rows, nn);

    // construct an randomized kd-tree index using 4 kd-trees
    Index<L2<float> > index(dataset, flann::LinearIndexParams());
    index.buildIndex();                                                                                               

    // do a knn search, using 128 checks
    index.radiusSearch(query, indices, dists, eps, flann::SearchParams(128));
    //index.knnSearch(query, indices, dists, nn, flann::SearchParams(128));

    printf("DISTS\n");
    printf("%d %d\n", query.rows, query.cols);
    for (int i = 0; i < query.rows; i++) {
        for (size_t j = 0; j < dists[i].size(); j++) {
            printf("%f ", dists[i][j]);
        }
        printf("\n");
    }
    printf("INDICES\n");
    for (int i = 0; i < query.rows; i++) {
        for (size_t j = 0; j < indices[i].size(); j++) {
            printf("%d ", indices[i][j]);
        }
        printf("\n");
    }

    delete[] dataset.ptr();
    delete[] query.ptr();
    
    return 0;
}
