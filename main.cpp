#include <initializer_list>
#include <sstream>
#include <vector>

#include "headers.hpp"
#include "local_search.hpp"
#include "reading.hpp"
#include "solvers.hpp"

// func(distanceMatrix, costs) returning sequence of indices
using TSPSolverStarting = vector<int>(const vector<vector<int>>&, /*D*/
                                      const vector<int>&,         /* costs*/
                                      int                         /*starting*/
);

ofstream openOutFile(string_view filename, string_view prefix = "")
{
    fs::path path;
    if (prefix.empty()) {
        path = fs::path(filename);
    }
    else {
        path = fs::path(prefix) / filename;
    }
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
    return out;
}

void printResults(const vector<int>& sol, bool toFile = false,
                  string_view filename = "solution.txt",
                  string_view prefix = "data/results")
{
    if (toFile) {
        ofstream out = openOutFile(filename, prefix);
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(out, "\n"));
    }
    else {
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(std::cout, "\n"));
    }
}

namespace run {

static std::stringstream latextables;

std::string format_with_spaces(long long number)
{
    auto str = std::to_string(number);
    int n = str.length() - 3;
    while (n > 0) {
        str.insert(n, " ");
        n -= 3;
    }
    return str;
}

std::string format_with_spaces(double number)
{
    auto [integral, fractional] = std::div(static_cast<long long>(number * 1000), 1000ll);
    auto int_part = format_with_spaces(integral);
    return std::format("{}.{:03}", int_part, std::abs(fractional));
}

struct SolutionStats {
    vector<int> best_sol;
    int best_sol_value = LARGE_SCORE;

    vector<int> worst_sol;
    int worst_sol_value = 0;

    long long average_numerator = 0;
    long long average_denominator = 0;

    [[nodiscard]] double average() const
    {
        return (double)(average_numerator) / average_denominator;
    }

    void track(const vector<int>& sol, int value)
    {
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

    [[nodiscard]] string format_latex_3() const
    {
        return std::format("{:>7} & {:>12} & {:>7}",
                           format_with_spaces((long long)this->best_sol_value),
                           format_with_spaces(this->average()),
                           format_with_spaces((long long)this->worst_sol_value));
    }
    [[nodiscard]] string format_latex_one_field() const
    {
        return std::format("{:>12} ({:>7} - {:>7})", format_with_spaces(this->average()),
                           format_with_spaces((long long)this->best_sol_value),
                           format_with_spaces((long long)this->worst_sol_value));
    }
};

void evalWithStarting(const TSPInstance& instance, TSPSolverStarting solver,
                      string_view solver_name, string_view instance_name)
{
    SolutionStats stats;
    Stopwatch tic;

    for (int i = 0; i < 200; i++) {
        // starting node index i
        auto sol = solver(instance.distances, instance.costs, i);
        int value = instance.evaluateSolution(sol);
        stats.track(sol, value);
    }

    cout << std::format("Took {}\n", tic.pretty_print());

    if (instance_name.ends_with(".csv")) {
        instance_name.remove_suffix(4);
    }
    latextables << std::format("{:12} & {} \\\\\n", solver_name, stats.format_latex_3());

    cout << std::format(
        "   best: {:>7}\n"
        "average: {:>11.3f}\n"
        "  worst: {:>7}\n",
        stats.best_sol_value, stats.average(), stats.worst_sol_value);
    printResults(stats.best_sol, true,
                 std::format("{}_{}_best.txt", instance_name, solver_name));
    printResults(stats.worst_sol, true,
                 std::format("{}_{}_worst.txt", instance_name, solver_name));
}

void weightedExperiments(const TSPInstance& inst)
{
    for (int param = 1; param < 10; param++) {
        float objective_weight = static_cast<float>(param) / 10;
        cerr << std::format("\x1b[32m {} with weight = {} \x1b[0m", "weighted",
                            objective_weight)
             << endl;

        vector<int> best_sol;
        int best_sol_value = LARGE_SCORE;

        vector<int> worst_sol;
        int worst_sol_value = 0;

        long long average_numerator = 0;
        long long average_denominator = 0;

        // Generating 200 random solutions, keeping the best and the worst ones.
        for (int i = 0; i < 200; i++) {
            auto sol =
                weightedSum2RegretTSP(inst.distances, inst.costs, i, objective_weight);
            int value = inst.evaluateSolution(sol);
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
        printResults(best_sol, true,
                     std::format("{}_{}_{}_best.txt", "TSPA", "weighted", param));
        printResults(worst_sol, true,
                     std::format("{}_{}_{}_worst.txt", "TSPA", "weighted", param));
    }
}

void main_assignment_2(const TSPInstance& inst, string_view input_file_name)
{
    // note: std::array<> does not deduce its length here
    const initializer_list<pair<TSPSolverStarting*, const char*>> WHAT_TO_RUN{
        {randomTSP, "random"},
        {nearestNeighborTSP, "NN-end"},
        {nearestNeighborAnyTSP, "NN-any"},
        {greedyCycleTSP, "greedyCycle"},
        {regret2TSP, "regret"},
        {[](const vector<vector<int>>& d, const vector<int>& costs, int starting_node) {
             // hard-coded weight
             return weightedSum2RegretTSP(d, costs, starting_node, 0.5);
         },
         "w-regret"},
    };

    latextables.clear();
    for (auto it : WHAT_TO_RUN) {
        auto [solver, name] = it;
        cerr << std::format("\x1b[32m {} \x1b[0m", name) << endl;
        // latextables << "  method " << name << endl;
        evalWithStarting(inst, solver, name, input_file_name);
    }
    cerr << "Results for " << input_file_name << ":\n" << latextables.view() << endl;
    weightedExperiments(inst);
}

void eval_local_search(const TSPInstance& inst, string_view input_file_name,
                       const vector<int>& sol1)
{
    const initializer_list<pair<decltype(&steepestLocalSearch), const char*>> LS_TYPES{
        {steepestLocalSearch, "steepest"},
        {greedyLocalSearch, "greedy"},
    };
    const initializer_list<pair<bool, const char*>> INTRA_SWAP_TYPES{{false, "nodes"},
                                                                     {true, "edges"}};

    Stopwatch timer;
    latextables.clear();

    for (auto [func, funcname] : LS_TYPES) {
        // check types using a compile error:
        // func.XD;
        for (auto [swaptype, movename] : INTRA_SWAP_TYPES) {
            timer.reset();
            auto solution = func(sol1, inst.distances, inst.costs, swaptype);
            cout << std::format(
                "\x1b[32m\t {} with {} type: {} \x1b[0m\n", funcname, movename,
                format_with_spaces((long long)inst.evaluateSolution(solution)));
            cout << std::format("Took {}\n", timer.pretty_print());
        }
    }
}

void main_local_search(const TSPInstance& inst, string_view input_file_name)
{
    Stopwatch tic;

    const auto initial_random = randomTSP(inst.distances, inst.costs, 123);

    const auto initial_gr = weightedSum2RegretTSP(inst.distances, inst.costs, 0, 0.5);

    cout << "Old score (weighted 2-regret): " << inst.evaluateSolution(initial_gr)
         << endl;
    eval_local_search(inst, input_file_name, initial_gr);

    cout << std::format("Everything took {}\n", tic.pretty_print());

    // printResults(sol2_nodes, true, "local_test.txt");
}

}  // namespace run

int main(int argc, char* argv[])
{
    const char* input_file_name = "TSPA.csv";
    if (argc >= 2) {
        const auto arg = argv[1];
        // if (fs::exists(arg)) {
        input_file_name = arg;
    }
    const auto inst = TSPInstance::readFromFile(input_file_name);

    // run::main_assignment_2(inst, input_file_name);

    run::main_local_search(inst, input_file_name);
    return 0;
}
