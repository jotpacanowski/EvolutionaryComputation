#include <algorithm>
#include <cstdint>
#include <vector>

#include "headers.hpp"

int intraSwapTwoNodesImpact(const vector<int> &solution, int id1, int id2,
                            const vector<vector<int>> &D, const int solution_size)
{
    int node1 = solution[id1];
    int node2 = solution[id2];

    int delta;

    int prev1 = id1 == 0 ? solution[(solution_size - 1)] : solution[id1 - 1];
    int next1 = id1 == (solution_size - 1) ? solution[0] : solution[id1 + 1];
    int prev2 = id2 == 0 ? solution[(solution_size - 1)] : solution[id2 - 1];
    int next2 = id2 == (solution_size - 1) ? solution[0] : solution[id2 + 1];

    if (next1 == node2) {  // id1+1=id2
        delta = -D[prev1][node1] + D[prev1][node2] - D[node2][next2] + D[node1][next2];
        return delta;
    }
    if (next2 == node1) {
        delta = -D[node1][next1] + D[node2][next1] - D[prev2][node2] + D[prev2][node1];
        return delta;
    }
    delta = -D[prev1][node1] + D[prev1][node2] - D[node1][next1] + D[node2][next1]
            - D[prev2][node2] + D[prev2][node1] - D[node2][next2] + D[node1][next2];
    return delta;
}

int intraSwapTwoEdgesImpact(const vector<int> &solution, int id1, int id2,
                            const vector<vector<int>> &D, const int solution_size)
{
    int node1 = solution[id1];
    int node2 = solution[id2];
    int delta, prev1, prev2, next1, next2;

    if (id1 < id2) {  //
        prev1 = id1 == 0 ? solution[(solution_size - 1)] : solution[id1 - 1];
        next2 = id2 == (solution_size - 1) ? solution[0] : solution[id2 + 1];

        delta = -D[prev1][node1] + D[prev1][node2] - D[node2][next2] + D[node1][next2];
    }
    else {  //(id1 > id2)
        next1 = id1 == (solution_size - 1) ? solution[0] : solution[id1 + 1];
        prev2 = id2 == 0 ? solution[(solution_size - 1)] : solution[id2 - 1];

        delta = -D[prev2][node2] + D[prev2][node1] - D[node1][next1] + D[node2][next1];
    }
    return delta;
}

int interSwapTwoNodesImpact(const vector<int> &solution, int idx, int external_node,
                            const vector<vector<int>> &D, const vector<int> &costs,
                            const int solution_size)
{
    int internal_node = solution[idx];

    int delta;

    int prev = idx == 0 ? solution[(solution_size - 1)] : solution[idx - 1];
    int next = idx == (solution_size - 1) ? solution[0] : solution[idx + 1];

    delta = -D[prev][internal_node] + D[prev][external_node] - D[internal_node][next]
            + D[external_node][next] - costs[internal_node] + costs[external_node];
    return delta;
}

int findIndex(const vector<int> &arr, int item)
{
    auto ret = std::find(arr.begin(), arr.end(), item);
    if (ret == arr.end()) return -1;
    return ret - arr.begin();
}

vector<int> localSearch(vector<int> solution, const vector<vector<int>> &distanceMatrix,
                        const vector<int> &costs, const int solution_size,
                        bool edges = false)
{
    int variable = 0;
    int delta;
    int highest_delta = 0;
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
    cout << "SIZE: " << solution_size << endl;

    for (int _ = 0; _ < 1000000; _++) {  // some arbitrary limit
        variable = rand() % 2;
        if (variable == 0) {  // Do intra moves
            for (int i1 = 0; i1 < solution_size; i1++) {
                for (int i2 = 0; i2 < i1; i2++) {
                    if (edges) {  // Swap two edges variation
                        delta = intraSwapTwoEdgesImpact(solution, i1, i2, distanceMatrix,
                                                        solution_size);
                    }
                    else {  // Swap two nodes variation
                        delta = intraSwapTwoNodesImpact(solution, i1, i2, distanceMatrix,
                                                        solution_size);
                    }

                    delta = -delta;               // WE WANT IMPROVEMENT -> smaller score
                    if (delta > highest_delta) {  // Steepest variation
                        highest_delta = delta;
                        pos1 = i1;
                        pos2 = i2;
                        found = true;
                    }
                }
            }
        }
        else {  // Do inter moves
            for (int internal : in_sol) {
                int idx = findIndex(solution, internal);
                for (int external : not_in_sol) {
                    delta = interSwapTwoNodesImpact(solution, idx, external,
                                                    distanceMatrix, costs, solution_size);
                    delta = -delta;               // WE WANT IMPROVEMENT -> smaller score
                    if (delta > highest_delta) {  // Steepest variation
                        highest_delta = delta;
                        best_external = external;
                        best_internal = internal;
                        pos1 = idx;
                        found = true;
                    }
                }
            }

            std::shuffle(not_in_sol.begin(), not_in_sol.end(), g);
            std::shuffle(in_sol.begin(), in_sol.end(), g);
        }
        if (found) {
            found = false;
            highest_delta = 0;

            if (variable == 0) {
                if (edges) {
                    reverse(solution.begin() + pos1, solution.begin() + pos2);
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

vector<int> localSearchGreedy(vector<int> solution,
                              const vector<vector<int>> &distanceMatrix,
                              const vector<int> &costs, const int solution_size,
                              bool edges = false)
{
    int variable = 0;
    int delta;
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

    cout << "SIZE: " << solution_size << endl;

    for (int _ = 0; _ < 1000000; _++) {  // some arbitrary limit
        variable = rand() % 2;
        if (variable == 0) {  // Do intra moves
            for (int i1 : id1_random) {
                if (found) break;
                for (int i2 : id2_random) {
                    if (i1 == i2) continue;
                    if (edges) {  // Swap two edges variation
                        delta = intraSwapTwoEdgesImpact(solution, i1, i2, distanceMatrix,
                                                        solution_size);
                    }
                    else {  // Swap two nodes variation
                        delta = intraSwapTwoNodesImpact(solution, i1, i2, distanceMatrix,
                                                        solution_size);
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
                    delta = interSwapTwoNodesImpact(solution, idx, external,
                                                    distanceMatrix, costs, solution_size);
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
                    reverse(solution.begin() + pos1, solution.begin() + pos2);
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