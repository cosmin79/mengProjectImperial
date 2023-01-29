#include <cstdio>
#include "IOUtil.hpp"
#include "OptParse.hpp"
#include "UtilGen.hpp"
#include <memory>
#include <fstream>
#include <chrono>
typedef Index<L2<double>> fIndex;

const string EPS_OPTION = "-EPS";
const string DEPTH_OPTION = "-DEPTH";
const string FILE_OPTION = "-FILE";
const string OUTPUT_DIR_OPTION = "-OUT_DIR";
const string VERBOSE_OPTION = "-VERBOSE";
const string STRATEGY_OPTION = "-STRATEGY";
enum Strategy {
    SPEED, ACCURACY, NO_OPTIONS
};

struct coordInfo {
    double minCoord, maxCoord;
};

class HyperPlane {
private:
    vector<coordInfo> dimInfo;
    int dims;

public:
    HyperPlane(int _dims, const pt<double>& p) : dims(_dims) {
        dimInfo.resize(dims);
        for (int dim = 0; dim < dims; dim++) {
            dimInfo[dim] = {p.coords[dim], p.coords[dim]};
        }
    }

    void updateWithPoint(const pt<double>& p) {
        for (int dim = 0; dim < dims; dim++) {
            dimInfo[dim].minCoord = min(dimInfo[dim].minCoord, p.coords[dim]);
            dimInfo[dim].maxCoord = max(dimInfo[dim].maxCoord, p.coords[dim]);
        }
    }

    bool mayBeNear(const pt<double>& p, double eps) {
        double minDist = 0.0;
        for (int dim = 0; dim < dims; dim++) {
            if (p.coords[dim] < dimInfo[dim].minCoord) {
                minDist += (dimInfo[dim].minCoord - p.coords[dim]) * (dimInfo[dim].minCoord - p.coords[dim]);
            } else if (p.coords[dim] > dimInfo[dim].maxCoord) {
                minDist += (dimInfo[dim].maxCoord - p.coords[dim]) * (dimInfo[dim].maxCoord - p.coords[dim]);
            } // otherwise it's in between and nothing can be derived ; there may be a point with the same coord
        }

        return minDist <= eps;
    }
};

class InputGenerator {
private:
    double eps;
    int depth, noDims;
    string fileName, outputDir;
    bool verbose;
    chrono::system_clock::time_point startTime;

    string outputFile;

    vector< pt<double> > dataset;
    Matrix<double> flannDataset;
    shared_ptr<fIndex> flannIndex;
    shared_ptr<InputSplitter> inputSplit;

    vector<int> whichPartition;
    vector<Partition> partitionDescription;

    Strategy strategy;
public:
    InputGenerator(double _eps, int _depth, string _fileName, string _outputDir, int _strategy, bool _verbose) : 
        eps(_eps), depth(_depth), fileName(_fileName), outputDir(_outputDir), verbose(_verbose) {
        strategy = static_cast<Strategy>(_strategy);

        parseFromFile(dataset, fileName, false);
        assert(!dataset.empty());
        noDims = dataset[0].coords.size();

        inputSplit = shared_ptr<InputSplitter>(new InputSplitter(depth, dataset));
        // Make sure file name to append at the end doesn't contain '/'
        fileName = fileName.substr(fileName.find_last_of("/\\") + 1);
        outputFile = outputDir + "/" + fileName;

        if (strategy == ACCURACY) {
            // init flann
            getFlannMatrix(flannDataset, dataset);
            flannIndex = shared_ptr<fIndex> (new fIndex(flannDataset, flann::KDTreeIndexParams(4)));
            flannIndex -> buildIndex();
        }
    }

    void generatePartitions() {
        partitionDescription = inputSplit -> getPartitions();
    }

    void writePartitions();

private:
    vector<int> regionQuery(pt<double>&, double);

    void writePartition(int);
};

vector<int> InputGenerator::regionQuery(pt<double>& p, double eps) {
    int cols = p.coords.size();

    Matrix<double> flannSingleQuery = Matrix<double>(new double[cols], 1, cols);
    copy(p.coords.begin(), p.coords.end(), flannSingleQuery[0]);

    vector<vector<int>> indices(1, vector<int>());
    vector<vector<double>> dists(1, vector<double>());

    flannIndex -> radiusSearch(flannSingleQuery, indices, dists, eps, flann::SearchParams(128));

    return indices[0];
}

void InputGenerator::writePartition(int partitionNo) {
    string path = outputFile + to_string(partitionNo);
    if (verbose) {
        printf("IDX %d -> %s\n", partitionNo, path.c_str());
    }
    ofstream out(path.c_str());

    vector<bool> markedOther(dataset.size(), false);
    int cntOthers = 0;
    if (strategy == SPEED) {
        assert(!partitionDescription[partitionNo].empty());
        HyperPlane partitionHP(noDims, partitionDescription[partitionNo][0]);
        for (auto& p: partitionDescription[partitionNo]) {
            partitionHP.updateWithPoint(p);
        }
        for (auto& p: dataset) {
            if (whichPartition[p.idx] != partitionNo && partitionHP.mayBeNear(p, eps)) {
                markedOther[p.idx] = true;
                cntOthers++;
            }
        }
    } else { // strategy = ACCURACY
        for (auto& p: partitionDescription[partitionNo]) {
            vector<int> epsNgbs = regionQuery(p, eps);
            for (auto& ngbIdx: epsNgbs) {
                if (whichPartition[ngbIdx] != partitionNo && !markedOther[ngbIdx]) {
                    markedOther[ngbIdx] = true;
                    cntOthers++;
                }
            }
        }
    }

    int totalNo = partitionDescription[partitionNo].size() + cntOthers;
    printf("(inputGen) Partition IDX %d -> noElems: %d, total (with phantom): %d\n",
        partitionNo, (int)partitionDescription[partitionNo].size(), totalNo);
    out << totalNo << " " << noDims << "\n";
    // points eps-reachable from other points in this partition
    for (auto& p: dataset) {
        if (whichPartition[p.idx] == partitionNo || markedOther[p.idx]) {
            out << pt<double>(p, whichPartition[p.idx]).toString() << "\n";
        }
    }

    out.close();
}

void InputGenerator::writePartitions() {
    whichPartition.resize(dataset.size());
    for (int idx = 0; idx < (1 << depth); idx++) {
        for (auto& p: partitionDescription[idx]) {
            whichPartition[p.idx] = idx;
        }
    }

    #pragma omp parallel for
    for (int idx = 0; idx < (1 << depth); idx++) {
        writePartition(idx);
    }
}

int main(int argc, char** argv) {
    chrono::system_clock::time_point startTime, currTime;
    startTime = chrono::system_clock::now();

    OptionParser optionParser(argc, argv);
    double eps = stod(optionParser.getOptionValue(EPS_OPTION));
    int depth = stoi(optionParser.getOptionValue(DEPTH_OPTION));
    string fileName = optionParser.getOptionValue(FILE_OPTION);
    string outputDir = optionParser.getOptionValue(OUTPUT_DIR_OPTION);
    int verbose = stoi(optionParser.getOptionValue(VERBOSE_OPTION));

    int strategy = stoi(optionParser.getOptionValue(STRATEGY_OPTION));
    int noOptions = (int)(Strategy::NO_OPTIONS);
    if (strategy < 0 || strategy > noOptions) {
        throw invalid_argument("Expected an index type between 0 and " + to_string(noOptions));
    } 

    InputGenerator generator(eps, depth, fileName, outputDir, strategy, verbose);
    currTime = std::chrono::system_clock::now();
    printf("(inputGen) Time after reading in milliseconds: %lld\n",
        chrono::duration_cast<chrono::milliseconds>(currTime - startTime).count());

    generator.generatePartitions();
    currTime = std::chrono::system_clock::now();
    printf("(inputGen) Time after generating partition points in milliseconds: %lld\n",
           chrono::duration_cast<chrono::milliseconds>(currTime - startTime).count());
    generator.writePartitions();
    currTime = std::chrono::system_clock::now();
    printf("(inputGen) Time after writing partitions in milliseconds: %lld\n",
           chrono::duration_cast<chrono::milliseconds>(currTime - startTime).count());
    return 0;
}
