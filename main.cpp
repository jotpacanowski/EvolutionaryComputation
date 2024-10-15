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

// TODO: proper linking
#include "reading.cpp"
#include "solvers.cpp"

// instead of `using namespace`
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

vector<int> nearestNeighborAnyTSP(vector<vector<int>> distanceMatrix, vector<int> costs,
                                  int starting_node)
{
    const int N = distanceMatrix.size();
    vector<int> solution;
    if (starting_node < 0) {
        starting_node = rand() % N;
    }
    solution.push_back(starting_node);
    int _last = starting_node;
    int next;
    for (int i = 0; i < distanceMatrix.size() / 2; i++) {
        // vector<int> candidates;
        int best_candidate = 0;
        int best_impact = 1000000;
        int candidate_index = 0;
        for (int j = 0; j < solution.size(); j++) {
            int candidate = nearestNeighbor(distanceMatrix, costs, solution[j], solution);
            int impact = distanceMatrix[solution[j]][candidate] + costs[candidate];
            if (impact < best_impact) {
                best_candidate = candidate;
                best_impact = impact;
                candidate_index = j;
            }
        }
        next = best_candidate;
        solution.insert(solution.begin() + candidate_index, best_candidate);
        _last = next;
    }
    return solution;
}

void printResults(vector<int> sol, bool toFile = false,
                  const char *filename = "solution.txt")
{
    if (toFile) {
        string path = "results/";
        path.append(filename);
        std::ofstream out(path);
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(out, "\n"));
    }
    else {
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(std::cout, "\n"));
    }
}

void randomTSPTask(vector<vector<int>> distanceMatrix, vector<int> costs)
{
    vector<int> sol;
    int value;

    vector<int> best_sol;
    int best_sol_value = 1000000;

    vector<int> worst_sol;
    int worst_sol_value = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        sol = randomTSP(distanceMatrix, costs);
        value = evaluateSolution(distanceMatrix, costs, sol);
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
    }

    cout << "BEST: " << best_sol_value << endl;
    printResults(best_sol, true, "random/best.txt");
    cout << "WORST: " << worst_sol_value << endl;
    printResults(worst_sol, true, "random/worst.txt");
}

void nearestNeighborTSPTask(vector<vector<int>> distanceMatrix, vector<int> costs)
{
    vector<int> sol;
    int value;

    vector<int> best_sol;
    int best_sol_value = 1000000;

    vector<int> worst_sol;
    int worst_sol_value = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        sol = nearestNeighborTSP(distanceMatrix, costs, i);
        value = evaluateSolution(distanceMatrix, costs, sol);
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
    }

    cout << "BEST: " << best_sol_value << endl;
    printResults(best_sol, true, "nearest/best.txt");
    cout << "WORST: " << worst_sol_value << endl;
    printResults(worst_sol, true, "nearest/worst.txt");
}

void nearestNeighborAnyTSPTask(vector<vector<int>> distanceMatrix, vector<int> costs)
{
    vector<int> sol;
    int value;

    vector<int> best_sol;
    int best_sol_value = 1000000;

    vector<int> worst_sol;
    int worst_sol_value = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        sol = nearestNeighborAnyTSP(distanceMatrix, costs, i);
        value = evaluateSolution(distanceMatrix, costs, sol);
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
    }

    cout << "BEST: " << best_sol_value << endl;
    printResults(best_sol, true, "nearest_any/best.txt");
    cout << "WORST: " << worst_sol_value << endl;
    printResults(worst_sol, true, "nearest_any/worst.txt");
}

int main()
{
    vector<vector<int>> instance = readInstance("TSPA.csv");
    vector<vector<int>> D = distanceMatrix(instance);
    vector<int> costs = getCosts(instance);
    // for (int i = 0; i < D.size(); i++)
    // {
    //     for (int j = 0; j < D.size(); j++)
    //     {
    //         std::cout << D[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    randomTSPTask(D, costs);
    nearestNeighborTSPTask(D, costs);
    nearestNeighborAnyTSPTask(D, costs);
    // vector<int> sol = nearestNeighborAnyTSP(D, costs, 0);
    // // vector<int>
    // //     sol = nearestNeighborTSP(D, costs, 0);
    // printResults(sol);
    // std::cout << evaluateSolution(D, costs, sol) << std::endl;

    return 0;
}
