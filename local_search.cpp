#include "local_search.hpp"

#include <algorithm>
#include <cstdint>
#include <vector>

#include "headers.hpp"

// TODO: g++ -DDEBUG=1 in terminal
#define DEBUG
#ifdef DEBUG
#include <cassert>
#else
#define assert(fact)
#endif

#include "local_search_helpers_inc.hpp"

// Intra-route move: swap two nodes within solution
int intraSwapTwoNodesImpact(const vector<int> &solution, const vector<vector<int>> &D,
                            int id1, int id2)
{
    int node1 = solution[id1];
    int node2 = solution[id2];

    int prev1 = cyclePrev(solution, id1);
    int next1 = cycleNext(solution, id1);
    int prev2 = cyclePrev(solution, id2);
    int next2 = cycleNext(solution, id2);

    if (next1 == node2) {
        int delta =
            // new
            +D[prev1][node2]
            + D[node1][next2]
            // remove
            - D[prev1][node1] - D[node2][next2];
        return delta;
    }
    if (next2 == node1) {
        int delta =
            // add
            +D[prev2][node1]
            + D[node2][next1]
            // remove
            - D[node1][next1] - D[prev2][node2];
        return delta;
    }
    int delta =
        // new edges
        +D[prev1][node2] + D[node2][next1] + D[prev2][node1]
        + D[node1][next2]
        // without node1
        - D[prev1][node1]
        - D[node1][next1]
        // without node2
        - D[prev2][node2] - D[node2][next2];
    return delta;
}

// Intra-route move: swap two edges within solution
int intraSwapTwoEdgesImpact(const vector<int> &solution, const vector<vector<int>> &D,
                            int id1, int id2)
{
    if (id1 > id2) {
        swap(id1, id2);
    }
    assert(id1 < id2);
    assert(id2 - id1 >= 1);
    assert((solution.size() + id1) - id2 >= 1);
    int node1 = solution[id1];
    int node2 = solution[id2];

    // assuming symmetric D
    int prev1 = cyclePrev(solution, id1);
    int next2 = cycleNext(solution, id2);

    int delta = -D[prev1][node1] + D[prev1][node2] - D[node2][next2] + D[node1][next2];
    return delta;
}

// Inter-route move:
// Exchange internal node with an external (i.e. outside the solution) node
int interSwapTwoNodesImpact(const vector<int> &solution, const vector<vector<int>> &D,
                            const vector<int> &costs, int idx, int external_node)
{
    int internal_node = solution[idx];

    int prev = cyclePrev(solution, idx);
    int next = cycleNext(solution, idx);

    int delta =
        // new edges
        +D[prev][external_node]
        + D[external_node][next]
        // existing edges
        - D[prev][internal_node]
        - D[internal_node][next]
        // cost
        - costs[internal_node] + costs[external_node];
    return delta;
}

int findIndex(const vector<int> &arr, int item)
{
    auto ret = std::find(arr.begin(), arr.end(), item);
    if (ret == arr.end()) return -1;
    return ret - arr.begin();
}

vector<int> steepestLocalSearch(vector<int> solution,
                                const vector<vector<int>> &distanceMatrix,
                                const vector<int> &costs, bool edges)
{
    const auto N = distanceMatrix.size();
    const auto solution_size = solution.size();
    assert(solution.size() == (N + 1) / 2);
    vector<int> in_sol;
    vector<int> not_in_sol;
    in_sol.reserve(solution_size);
    not_in_sol.reserve(N - solution_size);
    for (int i = 0; i < N; i++) {
        if (std::find(solution.begin(), solution.end(), i) != solution.end()) {
            in_sol.push_back(i);
        }
        else {
            not_in_sol.push_back(i);
        }
    }

    int pos1 = -5;
    int pos2 = -5;
    int highest_delta = 0;
    int best_external, best_internal;
    bool found = true;
    int _iters = 0;
    for (; found; _iters++) {
        found = false;
        bool intra_moves = true;
        {  // Do intra moves
            for (int i1 = 0; i1 < solution_size; i1++) {
                for (int i2 = 0; i2 < i1; i2++) {
                    int delta;
                    assert(i1 > i2);
                    if (edges) {
                        if (i2 == 0 && i1 == solution_size - 1)
                            continue;  // skip (0,99) - considering single edge
                        delta = intraSwapTwoEdgesImpact(solution, distanceMatrix, i1, i2);
                    }
                    else {
                        delta = intraSwapTwoNodesImpact(solution, distanceMatrix, i1, i2);
                    }

                    delta = -delta;
                    if (delta > highest_delta) {
                        highest_delta = delta;
                        pos1 = i1;
                        pos2 = i2;
                        found = true;
                    }
                }
            }
        }
        {  // Do inter-route moves
            for (int internal : in_sol) {
                int idx = findIndex(solution, internal);
                for (int external : not_in_sol) {
                    int delta = interSwapTwoNodesImpact(solution, distanceMatrix, costs,
                                                        idx, external);
                    delta = -delta;
                    if (delta > highest_delta) {
                        highest_delta = delta;
                        best_external = external;
                        best_internal = internal;
                        pos1 = idx;
                        found = true;
                        // also:
                        intra_moves = false;
                    }
                }
            }
        }
        if (found) {
            if (intra_moves == true) {
                if (edges) {
                    assert(pos1 != pos2);
                    if (pos1 < pos2) {
                        // past-the-end iterator: adding +1
                        reverse(solution.begin() + pos1, solution.begin() + pos2 + 1);
                    }
                    else if (pos2 < pos1) {
                        reverse(solution.begin() + pos2, solution.begin() + pos1 + 1);
                    }
                }
                else {
                    assert(pos1 != pos2);
                    iter_swap(solution.begin() + pos1, solution.begin() + pos2);
                }
            }
            else {
                assert(best_internal != best_external);
                solution[pos1] = best_external;
                std::erase(in_sol, best_internal);
                in_sol.push_back(best_external);
                std::erase(not_in_sol, best_external);
                not_in_sol.push_back(best_internal);
            }
            // Reset
            // found = false;
            highest_delta = 0;
        }
    }

#if 0
    cerr << "Did " << _iters << " iterations.\n";
#endif
    return solution;
}

vector<int> greedyLocalSearch(vector<int> solution,
                              const vector<vector<int>> &distanceMatrix,
                              const vector<int> &costs, bool edges)
{
    const auto solution_size = solution.size();
    assert(solution.size() == (distanceMatrix.size() + 1) / 2);
    int variable = 0;
    int pos1, pos2;
    bool found = false;
    int best_external, best_internal;
    vector<int> in_sol;
    vector<int> not_in_sol;
    not_in_sol.reserve(distanceMatrix.size() - solution_size);
    in_sol.reserve(solution_size);
    for (int i = 0; i < distanceMatrix.size(); i++) {
        if (std::find(solution.begin(), solution.end(), i) != solution.end()) {
            in_sol.push_back(i);
        }
        else {
            not_in_sol.push_back(i);
        }
    }
    int patience = 10;
    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(not_in_sol.begin(), not_in_sol.end(), g);
    std::shuffle(in_sol.begin(), in_sol.end(), g);

    vector<int> id1_random;
    vector<int> id2_random;
    id1_random.reserve(solution_size);
    id2_random.reserve(solution_size);
    for (int i = 0; i < solution_size; i++) {
        id1_random.push_back(i);
        id2_random.push_back(i);
    }
    std::shuffle(id1_random.begin(), id1_random.end(), g);
    std::shuffle(id2_random.begin(), id2_random.end(), g);

    for (int _ = 0; _ < 1000000; _++) {  // some arbitrary limit
        variable = rand() % 2;
        if (variable == 0) {  // Do intra moves
            for (int i1 : id1_random) {
                if (found) break;
                for (int i2 : id2_random) {
                    if (i1 == i2) continue;
                    int delta;
                    if (edges) {
                        if (i2 == 0 && i1 == solution_size - 1)
                            continue;  // skip (0,99) - considering single edge
                        if (i1 == 0 && i2 == solution_size - 1) continue;
                        delta = intraSwapTwoEdgesImpact(solution, distanceMatrix, i1, i2);
                    }
                    else {
                        delta = intraSwapTwoNodesImpact(solution, distanceMatrix, i1, i2);
                    }

                    delta = -delta;   // WE WANT IMPROVEMENT -> smaller score
                    if (delta > 0) {  // GREEDY variation
                        // cout << "GREEDY IMRPOVEMENT " << delta << endl;
                        pos1 = i1;
                        pos2 = i2;
                        found = true;
                        break;
                    }
                }
            }
        }
        else {  // Do inter moves
            for (int internal : in_sol) {
                if (found) break;
                int idx = findIndex(solution, internal);
                for (int external : not_in_sol) {
                    int delta = interSwapTwoNodesImpact(solution, distanceMatrix, costs,
                                                        idx, external);
                    delta = -delta;   // WE WANT IMPROVEMENT -> smaller score
                    if (delta > 0) {  // Steepest variation
                        // cout << "GREEDY IMRPOVEMENT " << delta << endl;
                        best_external = external;
                        best_internal = internal;
                        pos1 = idx;
                        found = true;
                        break;
                    }
                }
            }

            std::shuffle(not_in_sol.begin(), not_in_sol.end(), g);
            std::shuffle(in_sol.begin(), in_sol.end(), g);
        }
        if (found) {
            found = false;
            if (variable == 0) {
                if (edges) {
                    assert(pos1 != pos2);
                    if (pos1 < pos2) {
                        // past-the-end iterator: adding +1
                        reverse(solution.begin() + pos1, solution.begin() + pos2 + 1);
                    }
                    else if (pos2 < pos1) {
                        reverse(solution.begin() + pos2, solution.begin() + pos1 + 1);
                    }
                }
                else {
                    iter_swap(solution.begin() + pos1, solution.begin() + pos2);
                }
            }
            else {
                solution[pos1] = best_external;
                std::erase(in_sol, best_internal);
                in_sol.push_back(best_external);
                std::erase(not_in_sol, best_external);
                not_in_sol.push_back(best_internal);
            }
            patience = 10;
        }
        else {
            patience--;
            if (patience == 0) {
#if 0
                cout << "Run out of patience, no improvement for 10 moves. Got to "
                        "iterations: "
                     << _ << endl;
#endif
                return solution;
            }
        }
    }
    return solution;
}

SteepestLocalSearchWithCandidateMoves::SteepestLocalSearchWithCandidateMoves(
    const vector<vector<int>> &distanceMatrix, const vector<int> &costs,
    // default: 10
    int candidate_count)
{
    const int N = distanceMatrix.size();
    nearest10_only_distance.resize(N);
    assert(candidate_count < N);
    vector<pair<int, int>> nearest;
    nearest.reserve(N - 1);

    for (int i = 0; i < N; i++) {
        // reuse already allocated memory
        nearest.clear();
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            nearest.emplace_back(distanceMatrix[i][j], j);
        }
        // std::sort(nearest.begin(), nearest.end());
        std::partial_sort(nearest.begin(), nearest.begin() + candidate_count,
                          nearest.end());
        for (int j = 0; j < candidate_count; j++) {
            nearest10_only_distance[i].push_back(nearest[j].second);
        }
    }
    nearest10_objective.resize(N);
    for (int i = 0; i < N; i++) {
        nearest.clear();
        for (int j = 0; j < N; j++) {
            if (i == j) continue;
            nearest.emplace_back(distanceMatrix[i][j] + costs[j], j);
        }
        // std::sort(nearest.begin(), nearest.end());
        std::partial_sort(nearest.begin(), nearest.begin() + candidate_count,
                          nearest.end());
        for (int j = 0; j < candidate_count; j++) {
            nearest10_objective[i].push_back(nearest[j].second);
        }
    }
}

vector<int> SteepestLocalSearchWithCandidateMoves::do_local_search(
    vector<int> solution, const vector<vector<int>> &distanceMatrix,
    const vector<int> &costs, bool edges) const
{
    const auto N = distanceMatrix.size();
    const auto solution_size = solution.size();
    assert(solution.size() == (N + 1) / 2);

    vector<uint8_t> is_in_sol(N, 0);
    vector<int> in_sol;
    vector<int> not_in_sol;
    in_sol.reserve(solution_size);
    not_in_sol.reserve(N - solution_size);
    for (int i = 0; i < N; i++) {
        if (std::find(solution.begin(), solution.end(), i) != solution.end()) {
            in_sol.push_back(i);
            is_in_sol[i] = 1;
        }
        else {
            not_in_sol.push_back(i);
        }
    }

    int pos1 = -5;
    int pos2 = -5;
    int highest_delta = 0;
    int best_external, best_internal;
    bool found = true;
    int _iters = 0;
    for (; found; _iters++) {
        found = false;
        bool intra_moves = true;
        {  // Do intra-route moves
            for (int i1 = 0; i1 < solution_size; i1++) {
                // Interate only over the candidates
                // considering only nearest 10 vertices without the cost,
                // because score of solution already contains these costs
                for (int i2_node_id : nearest10_only_distance[solution[i1]]) {
                    int i2 = findIndex(solution, i2_node_id);
                    // if (is_in_sol[i2] == 0) continue;
                    if (i2 == -1) continue;
                    int delta;

                    // Candidate edge that we would like to introduce is:
                    // solution[i1] -> i2

                    if (edges) {
                        int a = min(i1, i2);
                        int b = max(i1, i2);
                        // assert(a < b);

                        int i3 = cycleIndexBefore(solution, b);
                        if ((i3 != a) && (!(i3 == 0 && a == solution_size - 1))
                            && (!(a == 0 && i3 == solution_size - 1))) {
                            delta =
                                intraSwapTwoEdgesImpact(solution, distanceMatrix, a, i3);
                            delta = -delta;
                            if (delta > highest_delta) {
                                highest_delta = delta;
                                pos1 = a;
                                pos2 = i3;
                                found = true;
                            }
                        }
                        int i4 = cycleIndexAfter(solution, a);
                        if ((i4 != b) && (!(i4 == 0 && b == solution_size - 1))
                            && (!(b == 0 && i4 == solution_size - 1))) {
                            delta =
                                intraSwapTwoEdgesImpact(solution, distanceMatrix, i4, b);
                            delta = -delta;
                            if (delta > highest_delta) {
                                highest_delta = delta;
                                pos1 = i4;
                                pos2 = b;
                                found = true;
                            }
                        }
                    }
                    else {  // nodes
                        int i3 = cycleIndexBefore(solution, i1);
                        delta = intraSwapTwoNodesImpact(solution, distanceMatrix, i2, i3);
                        delta = -delta;
                        if (delta > highest_delta) {
                            highest_delta = delta;
                            pos1 = i2;
                            pos2 = i3;
                            found = true;
                        }
                        int i4 = cycleIndexAfter(solution, i1);
                        delta = intraSwapTwoNodesImpact(solution, distanceMatrix, i2, i4);
                        delta = -delta;
                        if (delta > highest_delta) {
                            highest_delta = delta;
                            pos1 = i2;
                            pos2 = i4;
                            found = true;
                        }
                    }
                }
            }
        }
        {  // Do inter-route moves
            for (int internal : in_sol) {
                int idx = findIndex(solution, internal);
                assert(idx != -1);
                assert(is_in_sol[internal] == 1);
                // iterate only over the candidates
                for (int external : nearest10_objective[internal]) {
                    if (is_in_sol[external]) continue;
                    // Candidate edge: insert "external" before or after "internal"

                    int idx_bef = cycleIndexBefore(solution, idx);

                    int delta = interSwapTwoNodesImpact(solution, distanceMatrix, costs,
                                                        idx_bef, external);
                    delta = -delta;
                    if (delta > highest_delta) {
                        highest_delta = delta;
                        best_external = external;
                        best_internal = cyclePrev(solution, idx);
                        pos1 = idx_bef;
                        found = true;
                        // also:
                        intra_moves = false;
                    }

                    int idx_aft = cycleIndexAfter(solution, idx);
                    delta = interSwapTwoNodesImpact(solution, distanceMatrix, costs,
                                                    idx_aft, external);
                    delta = -delta;
                    if (delta > highest_delta) {
                        highest_delta = delta;
                        best_external = external;
                        best_internal = cycleNext(solution, idx);
                        pos1 = idx_aft;
                        found = true;
                        // also:
                        intra_moves = false;
                    }
                }
            }
        }
        if (found) {
            if (intra_moves == true) {
                if (edges) {
                    assert(pos1 != pos2);
                    if (pos1 < pos2) {
                        assert(pos2 - pos1 >= 1);
                        // past-the-end iterator: adding +1
                        reverse(solution.begin() + pos1, solution.begin() + pos2 + 1);
                    }
                    else if (pos2 < pos1) {
                        assert(pos1 - pos2 >= 1);
                        reverse(solution.begin() + pos2, solution.begin() + pos1 + 1);
                    }
                }
                else {
                    assert(pos1 != pos2);
                    iter_swap(solution.begin() + pos1, solution.begin() + pos2);
                }
            }
            else {
                assert(best_internal != best_external);
                solution[pos1] = best_external;
                std::erase(in_sol, best_internal);
                in_sol.push_back(best_external);
                std::erase(not_in_sol, best_external);
                not_in_sol.push_back(best_internal);
                is_in_sol[best_internal] = 0;
                is_in_sol[best_external] = 1;
            }
            // Reset
            // found = false;
            highest_delta = 0;
        }
    }
    return solution;
}
