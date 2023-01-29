#include <cstdio>
#include <string>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <chrono>
#include <vector>
#include <random>
#include <vector>
#include <cassert>
using namespace std;
const string file = "input/10D_500K_synthetic.ds";
const int MAX_CLUSTERS = 750;
const int MIN_CLUSTERS = 250;
const int MAX_DIM = 10000;
typedef vector<vector<int>> Cluster;

vector<Cluster> clustersComp;

int main() {
    freopen(file.c_str(), "w", stdout);

    srand(time(0));
    unsigned seed = chrono::system_clock::now().time_since_epoch().count();
    default_random_engine generator(seed);
    uniform_int_distribution<int> unifDistribution(1, MAX_DIM);

    int N = 500000, D = 10;
    int clusters = MIN_CLUSTERS + rand() % (MAX_CLUSTERS - MIN_CLUSTERS + 1);
    clustersComp.reserve(clusters);

    // generate cluster centers
    printf("%d %d\n", N, D);
    for (int i = 0; i < clusters; i++) {
        Cluster newCluster;
        vector<int> clusterCenter;
        for (int j = 0; j < D; j++) {
            clusterCenter.push_back(unifDistribution(generator));
        }
        newCluster.push_back(clusterCenter);
        clustersComp.push_back(newCluster);

        printf("%d ", i + 1);
        for (auto& x: clusterCenter) {
            printf("%d ", x);
        }
        printf("\n");
    }

    normal_distribution<double> normalDistribution(0.0, 8.0);
    for (int i = clusters; i < N; i++) {
        int closestCluster = rand() % clustersComp.size();
        int closestPtCluster = rand() % clustersComp[closestCluster].size();

        vector<int> newPt;
        printf("%d ", i + 1);
        for (int j = 0; j < D; j++) {
            newPt.push_back((int)(clustersComp[closestCluster][closestPtCluster][j] + normalDistribution(generator)));
            printf("%d ", newPt.back());
        }
        printf("\n");

        clustersComp[closestCluster].push_back(newPt);
    }
    return 0;
}