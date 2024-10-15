#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>
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

// func(distanceMatrix, costs) returning sequence of indices
using TSPSolver = vector<int>(vector<vector<int>>, vector<int>);
using TSPSolverStarting = vector<int>(vector<vector<int>>, vector<int>, int);

void printResults(vector<int> sol, bool toFile = false,
                  const char* filename = "solution.txt")
{
    if (toFile) {
        string path = "data/results/";
        path.append(filename);
        std::ofstream out(path);
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(out, "\n"));
    }
    else {
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(std::cout, "\n"));
    }
}

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
std::vector<int> greedyCycleTSP(std::vector<std::vector<int>> distanceMatrix,
                                std::vector<int> costs, int starting_node)
{
    const int N = distanceMatrix.size();
    std::vector<int> solution;
    if (starting_node < 0) {
        starting_node = rand() % N;
    }
    solution.push_back(starting_node);
    int candidate = nearestNeighbor(distanceMatrix, costs, starting_node, solution);
    solution.push_back(candidate);
    // std::cout<<solution[0] << " " << solution[1] << std::endl;
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
                // std::cout << "in solution: " << i << std::endl;
                continue;
            }
            for (int j = 0; j < solution.size();
                 j++)  // Check insertion for candidate at each index
            {
                if (j == 0) {
                    // std::cout <<"INSERT AT 0, BETWEEEN: "<<  solution[j] << ' ' <<
                    // candidate << solution[solution.size() - 1] << std::endl;

                    impact = distanceMatrix[solution[j]][candidate] +
                             distanceMatrix[solution.size() - 1][candidate] +
                             costs[candidate];
                }
                else {
                    // std::cout << solution[j - 1] << ' ' << candidate << solution[j] <<
                    // std::endl;

                    impact = distanceMatrix[solution[j]][candidate] +
                             distanceMatrix[solution[j - 1]][candidate] +
                             costs[candidate];
                }
                // std::cout<<"IMPACT: " << impact<<std::endl;
                if (impact < best_impact) {
                    best_impact = impact;
                    best_candidate = candidate;
                    // std::cout << "Better impact: " << best_impact << std::endl;

                    candidate_index = j;
                }
            }
        }

        solution.insert(solution.begin() + candidate_index, best_candidate);
    }
    return solution;
}

void evalWithStarting(vector<vector<int>> distanceMatrix, vector<int> costs,
                      TSPSolverStarting solver, const char* name)
{
    vector<int> sol;
    int value;

    vector<int> best_sol;
    int best_sol_value = 1'000'000;

    vector<int> worst_sol;
    int worst_sol_value = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        // sol = nearestNeighborTSP(distanceMatrix, costs, i);
        // sol = nearestNeighborAnyTSP(distanceMatrix, costs, i);
        // {       sol = greedyCycleTSP(distanceMatrix, costs, i);
        sol = solver(distanceMatrix, costs, i);
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
    std::stringstream ss;
    ss << name << "/best.txt";
    printResults(best_sol, true, ss.str().c_str());
    ss.clear();
    ss << name << "/worst.txt";
    cout << "WORST: " << worst_sol_value << endl;
    printResults(worst_sol, true, ss.str().c_str());
}

int main()
{
    vector<vector<int>> instance = readInstance("TSPA.csv");
    vector<vector<int>> D = distanceMatrix(instance);
    vector<int> costs = getCosts(instance);

    // for (int i = 0; i < D.size(); i++) {
    //     for (int j = 0; j < D.size(); j++) {
    //         cout << D[i][j] << ' ';
    //     }
    //     cout << std::endl;
    // }

    cerr << "\x1b[32m random solution \x1b[0m" << endl;
    // randomTSPTask(D, costs, randomTSP);
    evalWithStarting(D, costs, randomTSP, "random");
    cerr << "\x1b[32m  NN task\x1b[0m" << endl;
    // nearestNeighborTSPTask(D, costs);
    evalWithStarting(D, costs, nearestNeighborTSP, "NN-end");
    cerr << "\x1b[32m  NN any TSP\x1b[0m" << endl;
    // nearestNeighborAnyTSPTask(D, costs);
    evalWithStarting(D, costs, nearestNeighborAnyTSP, "NN-any");

    cerr << "\x1b[32m  greedy cycle\x1b[0m" << endl;
    evalWithStarting(D, costs, greedyCycleTSP, "greedyCycle");

    // vector<int> sol;
    // sol = nearestNeighborAnyTSP(D, costs, 0);
    // sol = nearestNeighborTSP(D, costs, 0);
    // printResults(sol);
    // cout << evaluateSolution(D, costs, sol) << endl;

    return 0;
}
