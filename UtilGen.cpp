#include "UtilGen.hpp"

void InputSplitter::pickBestMeanSplit(int& dimSplit, double& splitV, vector<double>& stDev, vector<double>& medianV) {
    // init with first dim
    int bestStDev = stDev[0];
    vector<int> possibleCandidates;
    possibleCandidates.push_back(0);
    for (int dim = 1; dim < noDims; dim++) {
        if (stDev[dim] > bestStDev) {
            bestStDev = stDev[dim];
            possibleCandidates.clear();
            possibleCandidates.push_back(dim); 
        } else if (abs(stDev[dim] - bestStDev) < EPS) {
            possibleCandidates.push_back(dim);
        }
    }

    // in case of equality pick randomly
    int idx = rand() % possibleCandidates.size();
    dimSplit = possibleCandidates[idx];
    splitV = medianV[dimSplit];
}

void InputSplitter::pickBestMeanSplit(int fromIdx, int toIdx, int& dimSplit, double& splitV, vector<double>& medianV) {
    vector<double> stDev(noDims, 0.0);
    for (int dim = 0; dim < noDims; dim++) {
        for (int i = fromIdx; i <= toIdx; i++) {
            stDev[dim] += (dataset[i].coords[dim] - medianV[dim]) * (dataset[i].coords[dim] - medianV[dim]);
        }
    }

    pickBestMeanSplit(dimSplit, splitV, stDev, medianV);
}

void InputSplitter::meanSplit(int fromIdx, int toIdx, int& dimSplit, double& splitV) {
    int midIdx = (fromIdx + toIdx) >> 1;
    vector<double> medianV(noDims, 0);
    for (int dim = 0; dim < noDims; dim++) {
        nth_element(dataset.begin() + fromIdx, dataset.begin() + midIdx,
            dataset.begin() + toIdx + 1, sortByDim[dim]);
        medianV[dim] = dataset[midIdx].coords[dim];
    }

    pickBestMeanSplit(fromIdx, toIdx, dimSplit, splitV, medianV);
}

void InputSplitter::partitionData(int fromIdx, int toIdx, int remDepth, int idx) {
    if (remDepth == 0) {
        partitionDescription[idx].reserve(toIdx - fromIdx + 1);
        for (int i = fromIdx; i <= toIdx; i++) {
            partitionDescription[idx].push_back(pt<double>(dataset[i], idx));
        }
        return ;
    }

    int dimSplit;
    double splitV;
    meanSplit(fromIdx, toIdx, dimSplit, splitV);

    int nextFreePos = fromIdx;
    for (int i = fromIdx; i <= toIdx; i++) {
        if (dataset[i].coords[dimSplit] < splitV) {
            swap(dataset[i], dataset[nextFreePos++]);
        }
    }

    partitionData(fromIdx, nextFreePos - 1, remDepth - 1, idx * 2);
    partitionData(nextFreePos, toIdx, remDepth - 1, idx * 2 + 1);
}