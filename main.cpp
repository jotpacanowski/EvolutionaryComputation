#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <filesystem>
#include <format>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <string>
#include <string_view>
#include <system_error>
#include <vector>

// TODO: proper linking
#include "reading.cpp"
#include "solvers.cpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::string_view;
using std::vector;

namespace fs = std::filesystem;

// func(distanceMatrix, costs) returning sequence of indices
using TSPSolverStarting = vector<int>(const vector<vector<int>>&, /*D*/
                                      const vector<int>&,         /* costs*/
                                      int                         /*starting*/
);

void printResults(const vector<int>& sol, bool toFile = false,
                  string_view filename = "solution.txt")
{
    if (toFile) {
        fs::path path = "data/results/";
        path = path / filename;
        std::ofstream out(path);
        if (!out.is_open()) {
            // perror("Error");
            std::error_code ec(errno, std::generic_category());
            cerr << "Error: " << ec.message() << " (code: " << ec.value() << ")\n";

            fs::path parent_dir = fs::path(path).parent_path();
            if (!fs::exists(parent_dir)) {
                cerr << "Creating missing directory: " << parent_dir << endl;
                std::error_code ec;
                fs::create_directories(parent_dir, ec);
                // Check error without exceptions
                if (ec.value() != 0) {
                    cerr << "Failed: " << ec.message() << endl;
                    exit(1);
                }
            }
        }
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(out, "\n"));
    }
    else {
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(std::cout, "\n"));
    }
}

void evalWithStarting(const TSPInstance& instance, TSPSolverStarting solver,
                      const char* name)
{
    vector<int> sol;

    vector<int> best_sol;
    int best_sol_value = 1'000'000;

    vector<int> worst_sol;
    int worst_sol_value = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        // sol = nearestNeighborTSP(distanceMatrix, costs, i);
        // sol = nearestNeighborAnyTSP(distanceMatrix, costs, i);
        // {       sol = greedyCycleTSP(distanceMatrix, costs, i);
        sol = solver(instance.distances, instance.costs, i);
        int value = instance.evaluateSolution(sol);
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
    }

    cout << "BEST : " << best_sol_value << endl;
    cout << "WORST: " << worst_sol_value << endl;
    printResults(best_sol, true, std::format("{}_best.txt", name));
    printResults(worst_sol, true, std::format("{}_worst.txt", name));
}

int main()
{
    const auto inst = TSPInstance::readFromFile("TSPA.csv");

    cerr << "\x1b[32m random solution \x1b[0m" << endl;
    evalWithStarting(inst, randomTSP, "random");
    cerr << "\x1b[32m  NN task\x1b[0m" << endl;
    evalWithStarting(inst, nearestNeighborTSP, "NN-end");
    cerr << "\x1b[32m  NN any TSP\x1b[0m" << endl;
    evalWithStarting(inst, nearestNeighborAnyTSP, "NN-any");
    cerr << "\x1b[32m  greedy cycle\x1b[0m" << endl;
    evalWithStarting(inst, greedyCycleTSP, "greedyCycle");

    return 0;
}
