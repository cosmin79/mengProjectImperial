#ifndef IO_UTILS_HPP_
#define IO_UTILS_HPP_

#include <vector>
#include <algorithm>
#include <flann/flann.hpp>
#include <sstream>
#include <string>
#include <fstream>
#include <cassert>

using namespace std;
using namespace flann;

const int UNDEF_PARTITION = -1;

template <typename T>
struct pt {
    int idx, partition, localIdx;
    vector<T> coords; 

    pt(int _idx, vector<T>& _coords) : idx(_idx), coords(_coords) {
        partition = UNDEF_PARTITION;
    }

    pt(const pt<T>& other, int _partition) : idx(other.idx), coords(other.coords) {
        partition = _partition;
    }

    pt(int _idx, vector<T>& _coords, int _partition) : idx(_idx), coords(_coords) {
        partition = _partition;
    }

    void setLocalIndex(int _localIdx) {
        localIdx = _localIdx;
    }

    string toString() {
        ostringstream out;
        out << idx; 
        if (partition != UNDEF_PARTITION) {
            out << " " << partition;
        }
        for (auto& coord: coords) {
            out << " " << coord;
        }

        return out.str();
    }
};

inline void parseDimensions(string& s, int& rows, int& cols) {
    stringstream ss(s);
    ss >> rows >> cols;
}

template <typename T>
inline pt<T> parsePoint(string& line, bool partitioned) {
    stringstream ss(line);
    int idx, partition;
    T elem;
    vector<T> coords;
    
    ss >> idx;
    if (partitioned) {
        ss >> partition;
    }
    while (ss >> elem) {
        coords.push_back(elem);
    }

    if (partitioned) {
        // idx already starts from 0
        return pt<T>(idx, coords, partition);
    }
    return pt<T>(idx - 1, coords);
}

template <typename T>
void parseFromFile(vector< pt<T> >& dataset, const string& name, bool partitioned) {
    ifstream in(name);
    string line;
    int noPoints, noDims;

    getline(in, line);
    parseDimensions(line, noPoints, noDims);
    dataset.reserve(noPoints);

    while (getline(in, line)) {
        dataset.push_back(parsePoint<T>(line, partitioned));
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

#endif