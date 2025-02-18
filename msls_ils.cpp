#include "msls_ils.hpp"

#include <algorithm>
#include <chrono>
#include <cstdlib>
#include <numeric>
#include <set>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "headers.hpp"
#include "local_search.hpp"
#include "solvers.hpp"

static inline int _evaluate_solution(const vector<int>& solution,
                                     const vector<vector<int>>& D, const vector<int>& C)
{
    int result = C[solution[0]];
    for (int i = 1; i < solution.size(); i++) {
        int node = solution[i];
        int prev_node = solution[i - 1];
        result += D[prev_node][node] + C[node];
    }
    result += D[solution.front()][solution.back()];
    return result;
}

vector<int> generate_random_solution_sliding_window(const vector<vector<int>>& D,
                                                    const vector<int>& C, int seed)
{
    const int N = D.size();
    // solution size
    const int S = (N + 1) / 2;
    vector<int> solution;
    solution.reserve(N);
    for (int i = 0; i < N; i++) {
        solution.push_back(i);
    }
    // std::random_device rd;
    // std::mt19937 g(rd());
    // deterministic:
    std::mt19937 g(3 * seed);

    std::shuffle(solution.begin(), solution.end(), g);

    // int score = _evaluate_solution(solution, distanceMatrix, costs);
    int score = C[solution[0]];
    for (int i = 1; i < S; i++) {
        int node = solution[i];
        int prev_node = solution[i - 1];
        score += D[prev_node][node] + C[node];
    }
    score += D[solution[0]][solution[S - 1]];
    // try to find better slice of "solution" than indices 0 .. 99
    int besti = 0;
    int bestscore = score;
    for (int i = 1; i < N - S; i++) {
        // without first
        score -= C[solution[0 + i - 1]];
        score -= D[solution[0 + i - 1]][solution[i]];
        score -= D[solution[0 + i - 1]][solution[S - 1 + i - 1]];
        // add next
        auto next = i + S - 1;
        score += C[solution[next]];
        score += D[solution[next]][solution[S - 1 + i - 1]];
        score += D[solution[next]][solution[i]];

        if (score < bestscore) {
            besti = i;
            bestscore = score;
        }
    }

    if (besti > 0) {
        for (int i = 0; i < S; i++) {
            solution[i] = solution[besti + i];
        }
    }
    solution.resize(S);
    solution.shrink_to_fit();
    return solution;
}

vector<int> multiple_start_steepestLS(const vector<vector<int>>& distanceMatrix,
                                      const vector<int>& costs, int seed)
{
    vector<int> best;
    int bestscore = LARGE_SCORE;
    int improvement = 0;
    int* ls_iterations = new int();
    int avg_ls_iterations = 0;
    for (int iters = 0; iters < 200; iters++) {
        auto initial =
            generate_random_solution_sliding_window(distanceMatrix, costs, seed * iters);

        auto sol = steepestLocalSearch(std::move(initial), distanceMatrix, costs, true,
                                       ls_iterations);
        //
        avg_ls_iterations += (*ls_iterations);
        int score = _evaluate_solution(sol, distanceMatrix, costs);
        if (score < bestscore) {
            improvement++;
            bestscore = score;
            best = std::move(sol);
        }
    }
    cout << 200 << "\t& " << improvement << "\t& " << (double)(avg_ls_iterations) / 200.0
         << "\\\\ \\hline" << endl;

    return best;
}

vector<int> get_n_nodes_not_in_sol(const vector<int>& solution, int n, int seed)
{
    unordered_set<int> solutionSet(solution.begin(), solution.end());
    vector<int> randomNumbers;
    randomNumbers.reserve(n);
    srand(seed);
    while (randomNumbers.size() < n) {
        int randomNum = rand() % 200;
        if (solutionSet.find(randomNum) == solutionSet.end()) {
            solutionSet.insert(randomNum);
            if (solutionSet.size() != randomNumbers.size()) {
                randomNumbers.push_back(randomNum);
            }
        }
    }

    return randomNumbers;
}

// Perturbation - random choice of random few methods
vector<int> perturb(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
                    vector<int>& solution, int seed)
{
    srand(seed);
    int n = rand() % 3;
    vector<int> perturbed = solution;

    // Swap k (3-8)  random nodes from outside solution
    if (n == 0) {
        int k = rand() % 6 + 3;
        vector<int> v = get_n_nodes_not_in_sol(perturbed, k, seed);
        for (auto node : v) {
            int location = rand() % solution.size();
            perturbed[location] = node;
        }
    }
    // Swap k (3-8)  random nodes from inside solution
    else if (n == 1) {
        int k = rand() % 6 + 3;
        for (int i = 0; i < k; i++) {
            int pos1 = rand() % perturbed.size();
            int pos2 = rand() % perturbed.size();
            int temp = perturbed[pos1];
            perturbed[pos1] = perturbed[pos2];
            perturbed[pos2] = temp;
        }
    }
    // Swap k (3-8) edges inside solution
    else if (n == 2) {
        int k = rand() % 6 + 3;
        for (int i = 0; i < k; i++) {
            int pos1 = rand() % perturbed.size();
            int pos2 = rand() % perturbed.size();
            if (pos1 < pos2) {
                // past-the-end iterator: adding +1
                reverse(perturbed.begin() + pos1, perturbed.begin() + pos2 + 1);
            }
            else if (pos2 < pos1) {
                reverse(perturbed.begin() + pos2, perturbed.begin() + pos1 + 1);
            }
        }
    }

    return perturbed;
}

vector<int> iterative_steepest_LS(const vector<vector<int>>& distanceMatrix,
                                  const vector<int>& costs, int seed)
{
    vector<int> sol =
        generate_random_solution_sliding_window(distanceMatrix, costs, seed + 1);
    sol = steepestLocalSearch(std::move(sol), distanceMatrix, costs, true);
    int bestscore = _evaluate_solution(sol, distanceMatrix, costs);
    int improvement = 0;
    auto now = std::chrono::steady_clock::now;
    using namespace std::chrono_literals;
    auto work_duration = 1336688.95us;
    auto start = now();
    int a = 1;
    int* ls_iterations = new int();
    int avg_ls_iterations = 0;
    while ((now() - start) < work_duration) {
        auto perturbed = perturb(distanceMatrix, costs, sol, seed * a++);
        // return sol;
        auto sol2 = steepestLocalSearch(std::move(perturbed), distanceMatrix, costs, true,
                                        ls_iterations);
        int score = _evaluate_solution(sol2, distanceMatrix, costs);
        // cout<<*ls_iterations<<endl;
        avg_ls_iterations += (*ls_iterations);
        // break;
        if (score < bestscore) {
            improvement++;
            bestscore = score;
            sol = std::move(sol2);
        }
    }

    cout << a << "\t& " << improvement << "\t& "
         << (double)(avg_ls_iterations) / (double)(a) << "\\\\ \\hline" << endl;

    return sol;
}

vector<int> get_slice(vector<int>& vec, int start_index, int length)
{
    int end_index = min((int)(start_index + length), (int)vec.size());
    int get_from_beginning = 0;
    if ((start_index + length) > vec.size()) {
        get_from_beginning = (start_index + length) - vec.size();
    }

    vector<int> slice = vector<int>(vec.begin() + start_index, vec.begin() + end_index);
    slice.insert(slice.end(), vec.begin(), vec.begin() + get_from_beginning);

    return slice;
}
int weightedRandom(const vector<pair<int, int>>& arr, int seed)
{
    // Step 1: Calculate the prefix sum of weights
    vector<int> prefixSum(arr.size());
    iota(prefixSum.begin(), prefixSum.end(), 1);  // Weights are 1, 2, 3, ..., n
    for (size_t i = 1; i < prefixSum.size(); ++i) {
        if (arr[i].first == 0) {
            prefixSum[i] = 0;
            continue;
        }
        prefixSum[i] += prefixSum[i - 1];
    }

    // Step 2: Generate a random number in the range [1, total weight]
    int totalWeight = prefixSum.back();
    mt19937 gen(seed);  // Mersenne Twister RNG
    uniform_int_distribution<> dis(1, totalWeight);
    int randomValue = dis(gen);

    // Step 3: Find the index corresponding to the random value
    auto it = lower_bound(prefixSum.begin(), prefixSum.end(), randomValue);
    return distance(prefixSum.begin(), it);  // Convert iterator to index
}

int get_worst_slice_randomized(const vector<vector<int>>& distanceMatrix,
                               const vector<int>& costs, vector<int>& solution,
                               int length, vector<int> taken_indices, int seed)
{
    int score;
    vector<pair<int, int>> all_scores;
    for (int start_index = 0; start_index < solution.size(); start_index++) {
        if (taken_indices.size() > 0) {
            if (std::find(taken_indices.begin(), taken_indices.end(), start_index)
                != taken_indices.end()) {
                score = 0;
            }
            else {
                vector<int> slice = get_slice(solution, start_index, length);
                score = _evaluate_solution(slice, distanceMatrix, costs);
            }
        }
        else {
            vector<int> slice = get_slice(solution, start_index, length);
            score = _evaluate_solution(slice, distanceMatrix, costs);
        }
        pair<int, int> p = {score, start_index};
        all_scores.push_back(p);
    }
    sort(all_scores.begin(), all_scores.end(),
         [](const pair<int, int>& a, const pair<int, int>& b) {
             return a.first < b.first;
         });
    // for (const auto& p : all_scores) {
    //     std::cout << "Score: " << p.first << ", Index: " << p.second << '\n';
    // }
    int index = weightedRandom(all_scores, seed);
    // cout << "index " << index << endl;
    return all_scores[index].second;
}

// nn heuristic - adapted for ils
vector<int> greedyPath(vector<int>& vec, const vector<vector<int>>& D,
                       const vector<int>& C)
{
    vector<int> sol;
    const int N = D.size();
    // stupid hack, trick to not consider all nodes except candidates
    vector<uint8_t> is_in_sol(N, 1);
    for (auto i : vec) {
        is_in_sol[i] = 0;
    }
    int starting_node = vec[0];
    sol.push_back(starting_node);
    is_in_sol[starting_node] = 1;

    int last = starting_node;
    for (int _ = 1; _ < vec.size(); _++) {
        int next = findNearestNeighbor(D, C, last, sol, is_in_sol);
        sol.push_back(next);
        is_in_sol[next] = 1;
        last = next;
    }
    return sol;
}

// Greedy Cycle Heuristic - adapted for ils
vector<int> greedyCycleRepair(vector<int> vec, const vector<vector<int>>& D,
                              const vector<int>& C, int desired_length)
{
    const int N = D.size();
    vector<uint8_t> is_in_sol(N, 0);
    for (auto i : vec) {
        is_in_sol[i] = 1;
    }

    while (vec.size() < desired_length) {
        int best_candidate = 0;
        int best_impact = LARGE_SCORE;
        int candidate_index = 0;
        int impact;
        for (int candidate = 0; candidate < N; candidate++)  // For each candidate
        {
            if (is_in_sol[candidate]) continue;
            for (int j = 0; j < vec.size(); j++)
            // Check insertion for candidate at each index
            {
                auto before = vec[(j == 0) ? (vec.size() - 1) : (j - 1)];  // node before
                auto after = vec[j];                                       // node after
                impact = D[after][candidate] + D[before][candidate] - D[before][after]
                         + C[candidate];

                if (impact < best_impact) {
                    best_impact = impact;
                    best_candidate = candidate;
                    candidate_index = j;
                }
            }
        }
        // cout << "BEST CANDIDATE " << best_candidate << " AT INDEX" << candidate_index
        //      << endl;
        vec.insert(vec.begin() + candidate_index, best_candidate);
        is_in_sol[best_candidate] = 1;
    }
    return vec;
}

vector<int> destroy_repair_non_deterministic_multiple_chains(
    const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
    vector<int>& solution, int seed)
{
    srand(seed);

    vector<int> perturbed = solution;
    vector<int> taken;
    // Total nodes changed - 10-40
    // Random number of chains, 2-4
    int k = rand() % 3 + 2;
    // Random length of chain, 5-10
    int n = rand() % 6 + 5;
    vector<int> to_remove;
    for (int i = 0; i < k; i++) {
        int start_index = get_worst_slice_randomized(distanceMatrix, costs, perturbed, n,
                                                     taken, seed * i);
        for (int j = -n; j < n; j++) {
            int idx = (start_index + j);
            if (idx < 0) {
                idx = solution.size() + idx;
            }
            if (idx > solution.size()) {
                idx = idx % solution.size();
            }
            taken.push_back(idx);
        }
        to_remove.push_back(start_index);
    }
    sort(to_remove.begin(), to_remove.end());
    for (int start_index : to_remove) {
        for (int i = start_index; i < start_index + n; i++) {
            int idx = i;
            if (i >= solution.size()) {
                idx = i % solution.size();
            }
            perturbed[idx] = -1;
        }
    }
    perturbed.erase(remove(perturbed.begin(), perturbed.end(), -1), perturbed.end());
    perturbed = greedyCycleRepair(perturbed, distanceMatrix, costs, solution.size());

    return perturbed;
}

vector<int> large_scale_neighbourhood_LS(const vector<vector<int>>& distanceMatrix,
                                         const vector<int>& costs, int seed)
{
    vector<int> sol =
        generate_random_solution_sliding_window(distanceMatrix, costs, seed + 1);
    sol = steepestLocalSearch(std::move(sol), distanceMatrix, costs, true);
    int bestscore = _evaluate_solution(sol, distanceMatrix, costs);
    int improvement = 0;
    auto now = std::chrono::steady_clock::now;
    using namespace std::chrono_literals;
    auto work_duration = 1336688.95us;
    auto start = now();
    int a = 1;
    int* ls_iterations = new int();
    int avg_ls_iterations = 0;
    while ((now() - start) < work_duration) {
        auto perturbed = destroy_repair_non_deterministic_multiple_chains(
            distanceMatrix, costs, sol, seed * a++);
        auto sol2 = steepestLocalSearch(std::move(perturbed), distanceMatrix, costs, true,
                                        ls_iterations);
        int score = _evaluate_solution(sol2, distanceMatrix, costs);
        // cout<<*ls_iterations<<endl;
        avg_ls_iterations += (*ls_iterations);
        // break;
        if (score < bestscore) {
            improvement++;
            bestscore = score;
            sol = std::move(sol2);
        }
    }
    cout << a << "\t& " << improvement << "\t& "
         << (double)(avg_ls_iterations) / (double)(a) << "\\\\ \\hline" << endl;

    return sol;
}

vector<int> large_scale_neighbourhood_LS2(const vector<vector<int>>& distanceMatrix,
                                          const vector<int>& costs, int seed)
{
    vector<int> sol =
        generate_random_solution_sliding_window(distanceMatrix, costs, seed + 1);
    sol = steepestLocalSearch(std::move(sol), distanceMatrix, costs, true);
    int bestscore = _evaluate_solution(sol, distanceMatrix, costs);
    int improvement = 0;
    auto now = std::chrono::steady_clock::now;
    using namespace std::chrono_literals;
    auto work_duration = 1336688.95us;
    auto start = now();
    int a = 1;
    int* ls_iterations = new int();
    int avg_ls_iterations = 0;
    while ((now() - start) < work_duration) {
        auto sol2 = destroy_repair_non_deterministic_multiple_chains(
            distanceMatrix, costs, sol, seed * a++);
        // auto sol2 = steepestLocalSearch(std::move(perturbed), distanceMatrix, costs,
        // true, ls_iterations);
        int score = _evaluate_solution(sol2, distanceMatrix, costs);
        // cout<<*ls_iterations<<endl;
        avg_ls_iterations += (*ls_iterations);
        // break;
        if (score < bestscore) {
            improvement++;
            bestscore = score;
            sol = std::move(sol2);
        }
    }
    cout << a << "\t& " << improvement << "\t& "
         << (double)(avg_ls_iterations) / (double)(a) << "\\\\ \\hline" << endl;

    return sol;
}
