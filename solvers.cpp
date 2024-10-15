#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <random>
#include <string>
#include <vector>

// instead of `using namespace`
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

vector<int> randomTSP(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
                      int)
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

int nearestNeighbor(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
                    int node, const vector<int>& solution)
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

vector<int> nearestNeighborTSP(const vector<vector<int>>& distanceMatrix,
                               const vector<int>& costs, int starting_node)
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

vector<int> nearestNeighborAnyTSP(const vector<vector<int>>& distanceMatrix,
                                  const vector<int>& costs, int starting_node)
{
    const int N = distanceMatrix.size();
    vector<int> solution;
    if (starting_node < 0) {
        starting_node = rand() % N;
    }
    solution.push_back(starting_node);
    int _last = starting_node;
    int next;
    for (int i = 0; i < N / 2; i++) {
        int best_candidate = 0;
        int best_impact = 1000000;
        int candidate_index = 0;
        for (int position = 0; position < solution.size(); position++) {
            int candidate =
                nearestNeighbor(distanceMatrix, costs, solution[position], solution);
            int impact = distanceMatrix[solution[position]][candidate] + costs[candidate];
            if (impact < best_impact) {
                best_candidate = candidate;
                best_impact = impact;
                candidate_index = position;
            }
        }
        next = best_candidate;
        solution.insert(solution.begin() + candidate_index, best_candidate);
        _last = next;
    }
    return solution;
}

vector<int> greedyCycleTSP(const vector<vector<int>>& distanceMatrix,
                           const vector<int>& costs, int starting_node)
{
    const int N = distanceMatrix.size();
    vector<int> solution;
    if (starting_node < 0) {
        starting_node = rand() % N;
    }
    solution.push_back(starting_node);
    int candidate = nearestNeighbor(distanceMatrix, costs, starting_node, solution);
    solution.push_back(candidate);
    // cout<<solution[0] << " " << solution[1] << endl;
    for (int i = 0; i < N / 2 - 1; i++) {
        int best_candidate = 0;
        int best_impact = 1000000;
        int candidate_index = 0;
        int impact;
        for (int candidate = 0; candidate < N; candidate++)  // For each candidate
        {
            if (std::count(solution.begin(), solution.end(),
                           candidate))  // If candidate already in solution, skip
            {
                // cout << "in solution: " << i << endl;
                continue;
            }
            for (int j = 0; j < solution.size();
                 j++)  // Check insertion for candidate at each index
            {
                if (j == 0) {
                    // cout <<"INSERT AT 0, BETWEEEN: "<<  solution[j] << ' ' <<
                    // candidate << solution[solution.size() - 1] << endl;

                    impact = distanceMatrix[solution[j]][candidate] +
                             distanceMatrix[solution.size() - 1][candidate] +
                             costs[candidate];
                }
                else {
                    // cout << solution[j - 1] << ' ' << candidate << solution[j] <<
                    // endl;

                    impact = distanceMatrix[solution[j]][candidate] +
                             distanceMatrix[solution[j - 1]][candidate] +
                             costs[candidate];
                }
                // cout<<"IMPACT: " << impact<<endl;
                if (impact < best_impact) {
                    best_impact = impact;
                    best_candidate = candidate;
                    // cout << "Better impact: " << best_impact << endl;

                    candidate_index = j;
                }
            }
        }

        solution.insert(solution.begin() + candidate_index, best_candidate);
    }
    return solution;
}
