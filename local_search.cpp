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


constexpr int cyclePrev(const vector<int> &solution, int index)
{
    if (index == 0) return solution[solution.size() - 1];
    return solution[index - 1];
}

constexpr int cycleNext(const vector<int> &solution, int index)
{
    if (index == solution.size() - 1) return solution[0];
    return solution[index + 1];
}

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
        tie(id1, id2) = tie(id2, id1);
    }
    assert(id1 < id2);
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
                                const vector<int> &costs, bool edges = false)
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
                    if (edges) {
                        if (i2 - i1 < 2) continue;                    // skip (0, 1)
                        if ((solution_size + i1) - i2 < 2) continue;  // skip (0,99)
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
        {  // Do inter moves
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
    cerr << "Did " << _iters << " iterations.\n";
    return solution;
}

vector<int> greedyLocalSearch(vector<int> solution,
                              const vector<vector<int>> &distanceMatrix,
                              const vector<int> &costs, bool edges = false)
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
                        if (i2 - i1 < 2) continue;                    // skip (0, 1)
                        if ((solution_size + i1) - i2 < 2) continue;  // skip (0,99)
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
                cout << "Run out of patience, no improvement for 10 moves. Got to "
                        "iterations: "
                     << _ << endl;
                return solution;
            }
        }
    }
    return solution;
}