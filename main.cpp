#include <cstring>
#include <functional>
#include <initializer_list>
#include <sstream>
#include <utility>
#include <vector>

#include "experiment_utils.hpp"
#include "headers.hpp"
#include "local_search.hpp"
#include "reading.hpp"
#include "solvers.hpp"

// func(distanceMatrix, costs) returning sequence of indices
using TSPSolverStarting = vector<int>(const vector<vector<int>>&, /*D*/
                                      const vector<int>&,         /* costs*/
                                      int                         /*starting*/
);

namespace run {

static std::stringstream latextables;

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
    // latextables << std::format("{:12} & {} \\\\\n", solver_name,
    // stats.format_latex_3());
    latextables << std::format(" {} & {} & \n", solver_name,
                               stats.format_latex_one_field());

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

// note: std::array<> does not deduce its length here
const initializer_list<pair<TSPSolverStarting*, const char*>> GREEDY_SOLVERS{
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

void main_assignment_2(const TSPInstance& inst, string_view input_file_name)
{
    latextables.clear();
    for (auto it : GREEDY_SOLVERS) {
        auto [solver, name] = it;
        cerr << std::format("\x1b[32m {} \x1b[0m", name) << endl;
        // latextables << "  method " << name << endl;
        evalWithStarting(inst, solver, name, input_file_name);
    }
    cerr << "Results for " << input_file_name << ":\n" << latextables.view() << endl;
}

void main_local_search(const TSPInstance& inst, string_view input_file_name)
{
    // Prepare lookup tables for nearest neighbors
    const SteepestLocalSearchWithCandidateMoves s_ls_candidate_moves(inst.distances,
                                                                     inst.costs, 10);

    using _Func = std::function<vector<int>(vector<int>, const vector<vector<int>>&,
                                            const vector<int>&, bool)>;
    const initializer_list<pair<_Func, const char*>> LS_TYPES{
        {[&s_ls_candidate_moves](vector<int> solution,
                                 const vector<vector<int>>& distanceMatrix,
                                 const vector<int>& costs, bool edges) {
             // call member function
             return s_ls_candidate_moves.do_local_search(std::move(solution),
                                                         distanceMatrix, costs, edges);
         },
         "St-LS-candidates"},
        {steepestLocalSearch, "steepest"},
        {greedyLocalSearch, "greedy"},
    };
    const initializer_list<pair<bool, const char*>> INTRA_SWAP_TYPES{{false, "nodes"},
                                                                     {true, "edges"}};

    Stopwatch timer;
    Stopwatch timer_all;
    latextables.clear();

    // for each Greedy / Steepest
    for (auto [localsearchfunc, funcname] : LS_TYPES) {
        for (auto [swaptype, movename] : INTRA_SWAP_TYPES) {
            for (int _is_random = 0; _is_random < 2; _is_random++) {
                bool random = (_is_random == 0);
                SolutionStats stats;
                for (auto [initsolver, initname] : GREEDY_SOLVERS) {
                    if (strncmp(initname, "random", 6) == 0) {
                        if (!random) continue;
                    }
                    else if (random || strcmp(initname, "w-regret") != 0) {
                        continue;
                    }

                    SolutionStats local_stats;

                    for (int starting = 0; starting < 200; starting++) {
                        const auto initial =
                            initsolver(inst.distances, inst.costs, starting);

                        timer.reset();
                        auto solution = localsearchfunc(initial, inst.distances,
                                                        inst.costs, swaptype);
                        auto t = timer.count_nanos() / 1000.0;
                        auto value = inst.evaluateSolution(solution);
                        stats.track(solution, value);
                        stats.add_time(t);
                    }
                    cout << std::format("{{initial {}\\\\ {} LS, {}}}  & {} & \n",
                                        initname, funcname, movename,
                                        stats.format_latex_one_field());
                    cout << std::format(
                        "Took {:.3f} s\n",
                        // sum
                        std::accumulate(stats.timings.begin(), stats.timings.end(), 0LL)
                            / 1e6);
                    // cout << "Timing summary (\u03BCs):\n";
                    cout << "Timing summary in us:\n";
                    print_summary(describe_vec(stats.timings));

                    // save to different folder
                    printResults(stats.best_sol, true,
                                 std::format("{}_{}_{}_{}_best.txt", input_file_name,
                                             funcname, initname, movename),
                                 "data/results-LS");
                }
                // cout << "count = " << stats.average_denominator << endl;
            }
        }
    }

    cout << std::format("Everything took {}\n", timer_all.pretty_print());

    // printResults(sol2_nodes, true, "local_test.txt");
}

}  // namespace run

int main(int argc, char* argv[])
{
    vector<string_view> input_file_names;
    if (argc >= 2) {
        const auto arg = argv[1];
        // if (fs::exists(arg)) {
        input_file_names.emplace_back(arg);
    }
    else {
        input_file_names.emplace_back("TSPA.csv");
        input_file_names.emplace_back("TSPB.csv");
    }

    for (auto input_file_name : input_file_names) {
        const auto inst = TSPInstance::readFromFile(input_file_name);
        cout << "\n";

        // run::main_assignment_2(inst, input_file_name);

        // run::weightedExperiments(inst);

        run::main_local_search(inst, input_file_name);

        cout << "\n\n";
    }
    return 0;
}
