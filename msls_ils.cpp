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
    for (int iters = 0; iters < 200; iters++) {
        auto initial =
            generate_random_solution_sliding_window(distanceMatrix, costs, seed * iters);

        auto sol = steepestLocalSearch(std::move(initial), distanceMatrix, costs, true);
        //
        int score = _evaluate_solution(sol, distanceMatrix, costs);
        if (score < bestscore) {
            bestscore = score;
            best = std::move(sol);
        }
    }

    return best;
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

int get_worst_slice(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
                    vector<int>& solution, int length)
{
    int best_id = 0;
    int WORST_SCORE = 0;
    int score;
    for (int start_index = 0; start_index < solution.size(); start_index++) {
        vector<int> slice = get_slice(solution, start_index, length);
        score = _evaluate_solution(slice, distanceMatrix, costs);
        if (score > WORST_SCORE) {
            WORST_SCORE = score;
            best_id = start_index;
        }
    }
    // cout << "best id: " << best_id << " score: " << WORST_SCORE;
    return best_id;
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
    // n = 1;
    vector<int> perturbed = solution;

    // Swap k (2-5) random nodes from outside solution
    if (n == 0) {
        int k = rand() % 4 + 2;
        vector<int> v = get_n_nodes_not_in_sol(perturbed, k, seed);
        for (auto node : v) {
            int location = rand() % solution.size();
            perturbed[location] = node;
        }
    }
    // Swap k (2-5) random nodes from inside solution
    else if (n == 1) {
        int k = rand() % 4 + 2;
        for (int i = 0; i < k; i++) {
            int pos1 = rand() % perturbed.size();
            int pos2 = rand() % perturbed.size();
            int temp = perturbed[pos1];
            perturbed[pos1] = perturbed[pos2];
            perturbed[pos2] = temp;
        }
    }
    // Swap k (2-5) edges inside solution
    else if (n == 2) {
        int k = rand() % 4 + 2;
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

    else if (n == 3) {
    }
    else if (n == 4) {
    }
    return perturbed;

    // return solution;
    vector<int> new_chain = get_n_nodes_not_in_sol(perturbed, n, seed);

    int start_index = get_worst_slice(distanceMatrix, costs, perturbed, n);
    int end_index = min((int)perturbed.size(), start_index + (int)new_chain.size());

    // change as much as we can from start to end
    copy(new_chain.begin(), new_chain.begin() + (end_index - start_index),
         perturbed.begin() + start_index);

    // if there is anything left, continue at the beginning (it's a chain)
    if (new_chain.size() > (end_index - start_index)) {
        copy(new_chain.begin() + (end_index - start_index), new_chain.end(),
             perturbed.begin());
    }

    return perturbed;
}

static inline int _evaluate_partial_path(const vector<int>& solution,
                                         const vector<vector<int>>& D)
{
    int result = 0;
    for (int i = 1; i < solution.size(); i++) {
        int node = solution[i];
        int prev_node = solution[i - 1];
        result += D[prev_node][node];
    }
    return result;
}
vector<int> exhaustiveSearchTSP(const vector<vector<int>>& distanceMatrix,
                                vector<int> vec)
{
    int minDistance = LARGE_SCORE;
    vector<int> best_sol;

    do {
        int currentDistance = _evaluate_partial_path(vec, distanceMatrix);
        if (currentDistance < minDistance) {
            minDistance = currentDistance;
            best_sol = vec;
        }
    } while (next_permutation(vec.begin(), vec.end()));

    return best_sol;
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
vector<int> greedyCycleLimited(vector<int> vec, const vector<vector<int>>& D,
                               const vector<int>& C)
{
    const int N = D.size();
    vector<int> sol;
    vector<uint8_t> is_in_sol(N, 1);
    for (auto i : vec) {
        is_in_sol[i] = 0;
    }
    int starting_node = vec[0];
    sol.push_back(starting_node);
    is_in_sol[starting_node] = 1;

    int candidate = findNearestNeighbor(D, C, starting_node, sol, is_in_sol);
    sol.push_back(candidate);
    is_in_sol[candidate] = 1;

    while (sol.size() < vec.size()) {
        int best_candidate = 0;
        int best_impact = LARGE_SCORE;
        int candidate_index = 0;
        int impact;
        for (auto candidate : vec)  // For each candidate
        {
            if (is_in_sol[candidate]) continue;
            for (int j = 0; j < sol.size(); j++)
            // Check insertion for candidate at each index
            {
                auto before = sol[(j == 0) ? (sol.size() - 1) : (j - 1)];  // node before
                auto after = sol[j];                                       // node after
                impact = D[after][candidate] + D[before][candidate] - D[before][after]
                         + C[candidate];

                if (impact < best_impact) {
                    best_impact = impact;
                    best_candidate = candidate;
                    candidate_index = j;
                }
            }
        }
        sol.insert(sol.begin() + candidate_index, best_candidate);
        is_in_sol[best_candidate] = 1;
    }
    return sol;
}

// Perturbation - rearrange worst n consecutive nodes with greedy cycle heuristic

vector<int> destroy_repair(const vector<vector<int>>& distanceMatrix,
                           const vector<int>& costs, vector<int>& solution, int seed)
{
    srand(seed);
    // pick n consecutive nodes (10, 20...90)
    int n = rand() % 5 + 1 * 10;
    vector<int> perturbed = solution;

    int start_index = get_worst_slice(distanceMatrix, costs, perturbed, n);
    vector<int> new_chain = get_slice(perturbed, start_index, n);
    new_chain = greedyPath(new_chain, distanceMatrix, costs);
    int end_index = min((int)perturbed.size(), start_index + (int)new_chain.size());

    // change as much as we can from start to end
    copy(new_chain.begin(), new_chain.begin() + (end_index - start_index),
         perturbed.begin() + start_index);

    // if there is anything left, continue at the beginning (it's a chain)
    if (new_chain.size() > (end_index - start_index)) {
        copy(new_chain.begin() + (end_index - start_index), new_chain.end(),
             perturbed.begin());
    }

    return perturbed;
}

bool isunique(vector<int> vec)
{
    set<int> uniqueElements(vec.begin(), vec.end());
    return uniqueElements.size() == vec.size();
}

vector<int> iterative_steepest_LS(const vector<vector<int>>& distanceMatrix,
                                  const vector<int>& costs, int seed)
{
    vector<int> sol =
        generate_random_solution_sliding_window(distanceMatrix, costs, seed + 1);
    int bestscore = _evaluate_solution(sol, distanceMatrix, costs);
    int improvement = 0;
    auto now = std::chrono::steady_clock::now;
    using namespace std::chrono_literals;
    auto work_duration = 1163518.25us;
    auto start = now();
    int a = 1;
    while ((now() - start) < work_duration) {
        auto perturbed = perturb(distanceMatrix, costs, sol, seed * a++);
        // return sol;
        auto sol2 =
            steepestLocalSearch(std::move(perturbed), distanceMatrix, costs, true);
        int score = _evaluate_solution(sol2, distanceMatrix, costs);

        if (score < bestscore) {
            improvement++;
            bestscore = score;
            sol = std::move(sol2);
        }
    }
    cout << "Improvement found: " << improvement << " times." << endl;
    cout << "Perturbation count: " << a << endl;
    return sol;
}
