#include "headers.hpp"
#include "reading.hpp"
#include "solvers.hpp"

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
                // Try again
                out.open(path);
                if (!out.is_open()) {
                    perror("second open");
                    exit(1);
                }
            }
            else {
                exit(1);
            }
        }
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(out, "\n"));
    }
    else {
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(std::cout, "\n"));
    }
}

void evalWithStarting(const TSPInstance& instance, TSPSolverStarting solver,
                      string_view solver_name, string_view instance_name)
{
    // vector<int> sol;

    vector<int> best_sol;
    int best_sol_value = LARGE_SCORE;

    vector<int> worst_sol;
    int worst_sol_value = 0;

    long long average_numerator = 0;
    long long average_denominator = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        auto sol = solver(instance.distances, instance.costs, i);
        int value = instance.evaluateSolution(sol);
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
        average_numerator += value;
        average_denominator += 1;
    }

    cout << std::format(
        "   best: {:>7}\n"
        "average: {:>11.3f}\n"
        "  worst: {:>7}\n",
        best_sol_value, (double)(average_numerator) / average_denominator,
        worst_sol_value);
    if (instance_name.ends_with(".csv")) {
        instance_name.remove_suffix(4);
    }
    printResults(best_sol, true,
                 std::format("{}_{}_best.txt", instance_name, solver_name));
    printResults(worst_sol, true,
                 std::format("{}_{}_worst.txt", instance_name, solver_name));
}

int main(int argc, char* argv[])
{
    const char* input_file_name = "TSPA.csv";
    if (argc >= 2) {
        const auto arg = argv[1];
        // if (fs::exists(arg)) {
        input_file_name = arg;
    }
    const auto inst = TSPInstance::readFromFile(input_file_name);

    const std::pair<TSPSolverStarting*, const char*> WHAT_TO_RUN[] = {
        {randomTSP, "random"},
        {nearestNeighborTSP, "NN-end"},
        {nearestNeighborAnyTSP, "NN-any"},
        {greedyCycleTSP, "greedyCycle"},
    };

    for (auto it : WHAT_TO_RUN) {
        auto [solver, name] = it;
        cerr << std::format("\x1b[32m {} \x1b[0m", name) << endl;
        evalWithStarting(inst, solver, name, input_file_name);
    }
    return 0;
}
