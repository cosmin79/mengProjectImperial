#ifndef UTIL_GEN_HPP_
#define UTIL_GEN_HPP_

#include "IOUtil.hpp"
#include <vector>
#include <algorithm>
#include <cassert>
#include <ctime>
#include <functional>
using namespace std;

typedef vector<pt<double>> Partition;
const double EPS = 1e-6;

class InputSplitter {
private:
    vector<pt<double>>& dataset;

    int depth, noDims;

    vector<function<bool(const pt<double>&, const pt<double>&)>> sortByDim;

    vector<Partition> partitionDescription;

public:
    InputSplitter(int _depth, vector<pt<double>>& _dataset) : dataset(_dataset), depth(_depth) {
        assert(!dataset.empty());
        noDims = dataset[0].coords.size();

        sortByDim.reserve(noDims);
        for (int i = 0; i < noDims; i++) {
            sortByDim.push_back([i](const pt<double>& a, const pt<double>& b) {
                return a.coords[i] < b.coords[i];
            });
        }
        partitionDescription.assign(1 << depth, vector<pt<double>>());

        // Answer to life
        srand(42);
    }

    vector<Partition> getPartitions() {
        partitionData(0, dataset.size() - 1, depth, 0);
        return partitionDescription;
    }

private:
    void pickBestMeanSplit(int&, double&, vector<double>&, vector<double>&);

    void pickBestMeanSplit(int, int, int&, double&, vector<double>&);

    void meanSplit(int, int, int&, double&);

    void partitionData(int, int, int, int);
};

#endif