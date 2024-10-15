#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>  // std::stringstream
#include <string>
#include <vector>

// TODO
// #include "reading.cpp"

// instead of `using namespace`
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

vector<int> randomTSP(vector<vector<int>> distanceMatrix, vector<int> costs, int)
{
    vector<int> solution;
    solution.reserve(distanceMatrix.size());
    for (int i = 0; i < distanceMatrix.size(); i++) {
        solution.push_back(i);
    }
    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(solution.begin(), solution.end(), g);
    return solution;
}

int nearestNeighbor(vector<vector<int>> distanceMatrix, vector<int> costs, int node,
                    vector<int> solution)
{
    vector<int> neighbors = distanceMatrix[node];
    int current_impact = 1000000;
    int nearest_neighbor = 0;
    int impact;
    for (int i = 0; i < neighbors.size(); i++) {
        // std::cout << "i: " << i << std::endl;
        if (std::count(solution.begin(), solution.end(), i)) {
            // std::cout << "in solution: " << i << std::endl;
            continue;
        }
        impact = distanceMatrix[node][i] + costs[i];
        // std::cout << impact << std::endl;
        if (impact < current_impact) {
            nearest_neighbor = i;
            current_impact = impact;
        }
    }
    // std::cout << nearest_neighbor << std::endl;
    return nearest_neighbor;
}

vector<int> nearestNeighborTSP(vector<vector<int>> distanceMatrix, vector<int> costs,
                               int starting_node)
{
    const int N = distanceMatrix.size();
    vector<int> solution;
    if (starting_node < 0) {
        starting_node = rand() % N;
    }
    solution.push_back(starting_node);
    int last = starting_node;
    int next;
    for (int i = 0; i < distanceMatrix.size() / 2; i++) {
        next = nearestNeighbor(distanceMatrix, costs, last, solution);
        solution.push_back(next);
        last = next;
    }
    return solution;
}
