#include "local_search.hpp"

#include <algorithm>
#include <cstdint>
#include <unordered_set>
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

vector<int> steepestLocalSearch(vector<int> solution,
                                const vector<vector<int>> &distanceMatrix,
                                const vector<int> &costs, bool edges, int *iterations)
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
    if (iterations != nullptr) {
        *iterations = _iters;
    }
    return solution;
}

vector<int> greedyLocalSearch(vector<int> solution,
                              const vector<vector<int>> &distanceMatrix,
                              const vector<int> &costs, bool edges, int *iterations)
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

// all ints in this struct are node IDs for the problem instance
// (the same representation for different solutions of an instance)
struct LSMove {
    int score_delta;
    bool is_edge_swap;
    // for internal-external:
    int pos1;
    int pos2;
    // meta-data for edge-swap validity:
    // (as in intraSwapTwoEdgesImpact)
    int before_a;
    int after_a;
    int before_b;
    int after_b;
};

inline bool operator<(const LSMove &lhs, const LSMove &rhs)
{
    return lhs.score_delta < rhs.score_delta;
}

inline bool operator==(const LSMove &lhs, const LSMove &rhs)
{
    // https://stackoverflow.com/questions/5740310/no-operator-found-while-comparing-structs-in-c
    return tie(lhs.score_delta, lhs.is_edge_swap, lhs.pos1, lhs.pos2, lhs.before_a,
               lhs.after_b)
           == tie(rhs.score_delta, rhs.is_edge_swap, rhs.pos1, rhs.pos2, rhs.before_a,
                  rhs.after_b);
}

struct CustomMoveHash {
    std::size_t operator()(const LSMove &x) const noexcept
    {
        return std::hash<int>()(x.score_delta) ^ std::hash<bool>()(x.is_edge_swap)
               ^ std::hash<int>()(x.pos1);
    }
};

#include "msls_ils.hpp"

vector<int> steepest_LS_LM(vector<int> solution,
                           const vector<vector<int>> &distanceMatrix,
                           const vector<int> &costs, bool edges)
{
    assert(edges);
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

    vector<int> node_id_to_solution(N);
    std::fill(node_id_to_solution.begin(), node_id_to_solution.end(), -1);
    for (int i = 0; i < solution_size; i++) {
        node_id_to_solution[solution[i]] = i;
    }

    vector<LSMove> improving_moves_list;
    improving_moves_list.reserve(15'000 * 2);

    int lowest_delta = 0;
    bool found = true;

    // phase 1 - init LM (list of moves)
    found = false;
    LSMove best_move;
    best_move.score_delta = 0;
    {
        {  // Do intra moves
            for (int i1 = 0; i1 < solution_size; i1++) {
                for (int i2 = 0; i2 < i1; i2++) {
                    // assert(i1 > i2);
                    // if (edges) {
                    if (i2 == 0 && i1 == solution_size - 1)
                        continue;  // skip (0,99) - considering single edge
                    int delta = intraSwapTwoEdgesImpact(solution, distanceMatrix, i1, i2);

                    LSMove currentmove = {.score_delta = delta,
                                          .is_edge_swap = true,
                                          .pos1 = solution[i2],
                                          .pos2 = solution[i1],
                                          .before_a = cyclePrev(solution, i2),
                                          .after_a = cycleNext(solution, i2),
                                          .before_b = cyclePrev(solution, i1),
                                          .after_b = cycleNext(solution, i1)

                    };
                    if (delta < 0) {
                        improving_moves_list.push_back(currentmove);
                    }
                    if (delta < lowest_delta) {
                        found = true;
                        best_move = currentmove;
                        lowest_delta = delta;
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
                    LSMove currentmove = {.score_delta = delta,
                                          .is_edge_swap = false,
                                          .pos1 = internal,
                                          .pos2 = external,

                                          .before_a = cyclePrev(solution, idx),
                                          .after_a = cycleNext(solution, idx),
                                          .before_b = -1,
                                          .after_b = -1};

                    if (delta < 0) {
                        improving_moves_list.push_back(currentmove);
                    }
                    if (delta < lowest_delta) {
                        found = true;
                        best_move = currentmove;
                        lowest_delta = delta;
                    }
                }
            }
        }

        std::sort(improving_moves_list.begin(), improving_moves_list.end());

        if (found) {
            // assert(found && improving_moves_list.front() == best_move);
            auto _i = std::find(improving_moves_list.begin(), improving_moves_list.end(),
                                best_move);
            assert(improving_moves_list.end() != _i);
            for (auto it = improving_moves_list.begin();
                 it != improving_moves_list.end() && it != _i; ++it) {
                assert(it->score_delta == best_move.score_delta);
            }
        }
    }

    // best_move = improving_moves_list.front();
    assert(best_move.score_delta == improving_moves_list[0].score_delta);
    best_move = improving_moves_list[0];

    // found = true;
    // lowest_delta = 0;
    // phase 2
    int _iters = 0;
    unordered_set<int> affected_node_ids;
    affected_node_ids.reserve(N);
    do {
        ++_iters;
        affected_node_ids.clear();

        /*
            auto try_edge_swap = [&](int i1, int i2) {
                // Do intra moves
                if (i2 > i1) swap(i1, i2);

                assert(i2 < i1);
                if (i2 == 0 && i1 == solution_size - 1)
                    return;  // skip (0,99) - considering single edge
                int delta = intraSwapTwoEdgesImpact(solution, distanceMatrix, i1, i2);

                LSMove currentmove = {.score_delta = delta,
                                      .is_edge_swap = true,
                                      .pos1 = solution[i2],
                                      .pos2 = solution[i1],
                                      .before_a = cyclePrev(solution, i2),
                                      .after_a = cycleNext(solution, i2),
                                      .before_b = cyclePrev(solution, i1),
                                      .after_b = cycleNext(solution, i1)

                };
                if (delta < 0) {
                    improving_moves_list.push_back(currentmove);
                }
                // if (delta < lowest_delta) {
                //     found = true;
                //     best_move = currentmove;
                //     lowest_delta = delta;
                // }
            };

            auto try_internal_external = [&](int internal, int external) {
                assert(findIndex(solution, internal) >= 0);
                assert(findIndex(solution, external) < 0);
                int idx = findIndex(solution, internal);
                int delta =
                    interSwapTwoNodesImpact(solution, distanceMatrix, costs, idx,
           external); LSMove currentmove = {.score_delta = delta, .is_edge_swap = false,
                                      .pos1 = internal,
                                      .pos2 = external,

                                      .before_a = cyclePrev(solution, idx),
                                      .after_a = cycleNext(solution, idx),
                                      .before_b = -1,
                                      .after_b = -1};

                if (delta < 0) {
                    improving_moves_list.push_back(currentmove);
                }
                // if (delta < lowest_delta) {
                //     found = true;
                //     best_move = currentmove;
                //     lowest_delta = delta;
                // }
            };
            */

        // int score_before = _evaluate_solution(solution, distanceMatrix, costs);

        // Apply best_move
        if (found) {
            if (best_move.is_edge_swap == true) {
                auto pos1 = node_id_to_solution[best_move.pos1];
                auto pos2 = node_id_to_solution[best_move.pos2];
                assert(pos1 != pos2);
                if (pos1 < pos2) {
                    int prev1 = cyclePrev(solution, pos1);
                    int next2 = cycleNext(solution, pos2);
                    assert(prev1 == best_move.before_a);
                    assert(next2 == best_move.after_b);
                    // past-the-end iterator: adding +1
                    reverse(solution.begin() + pos1, solution.begin() + pos2 + 1);

                    affected_node_ids.insert(prev1);
                    affected_node_ids.insert(best_move.pos1);
                    affected_node_ids.insert(best_move.pos2);
                    affected_node_ids.insert(next2);
                    assert(prev1 < 200);
                    assert(best_move.pos1 < 200);
                    assert(best_move.pos2 < 200);
                    assert(next2 < 200);
                }
                else if (pos2 < pos1) {
                    int prev1 = cyclePrev(solution, pos2);
                    int next2 = cycleNext(solution, pos1);
                    assert(prev1 == best_move.before_a);
                    assert(next2 == best_move.after_b);
                    reverse(solution.begin() + pos2, solution.begin() + pos1 + 1);

                    affected_node_ids.insert(prev1);
                    affected_node_ids.insert(best_move.pos2);
                    affected_node_ids.insert(best_move.pos1);
                    affected_node_ids.insert(next2);
                    assert(prev1 < 200);
                    assert(best_move.pos1 < 200);
                    assert(best_move.pos2 < 200);
                    assert(next2 < 200);
                }
            }
            else {
                int best_internal = best_move.pos1;
                int pos1 = node_id_to_solution[best_move.pos1];
                int best_external = best_move.pos2;
                assert(best_internal != best_external);
                solution[pos1] = best_external;
                std::erase(in_sol, best_internal);
                in_sol.push_back(best_external);
                std::erase(not_in_sol, best_external);
                not_in_sol.push_back(best_internal);

                // affected_node_ids.insert(best_internal);  // was removed
                affected_node_ids.insert(best_external);
                affected_node_ids.insert(cyclePrev(solution, pos1));
                affected_node_ids.insert(cycleNext(solution, pos1));
                assert(best_external < 200);
                assert(cyclePrev(solution, pos1) < 200);
                assert(cycleNext(solution, pos1) < 200);
            }
            // Re-compute
            //
            std::fill(node_id_to_solution.begin(), node_id_to_solution.end(), -1);
            for (int i = 0; i < solution_size; i++) {
                node_id_to_solution[solution[i]] = i;
            }

            // Now forget about best_move:
            std::erase(improving_moves_list, best_move);
        }

        // int score_after = _evaluate_solution(solution, distanceMatrix, costs);
        // assert(score_after == score_before + best_move.score_delta);
        // assert(score_after == score_before + lowest_delta);
        //
        lowest_delta = 0;
        found = false;

        // compute new moves
        // vector<decltype(improving_moves_list)::iterator> moves_to_remove;
        // vector<bool> moves_to_remove;
        // vector<uint8_t> moves_to_remove(improving_moves_list.size(), 0);
        std::unordered_set<LSMove, CustomMoveHash> moves_to_remove;
        // for (auto move : improving_moves_list) {
        // for (auto it = improving_moves_list.begin(); it != improving_moves_list.end();
        // ++it) {
        for (auto moveidx = 0; moveidx < improving_moves_list.size(); ++moveidx) {
            auto move = improving_moves_list[moveidx];
            if (move.is_edge_swap == false) {
                auto internal = move.pos1;
                auto external = move.pos2;

                auto ii = node_id_to_solution[internal];
                if (ii < 0) {
                    // moves_to_remove[moveidx] = 1;
                    moves_to_remove.insert(move);
                }

                auto i_bef = cycleIndexBefore(solution, ii);
                auto i_aft = cycleIndexAfter(solution, ii);

                if (solution[i_bef] != move.before_a || solution[i_aft] != move.after_a) {
                    moves_to_remove.insert(move);
                }

                // removed: before_a -- internal -- after_b
                // added:   ... -- external -- ...
            }
            else {
                // edge 1
                auto pa = move.before_a;
                auto a = move.pos1;

                auto i_pa = node_id_to_solution[pa];
                auto i_a = node_id_to_solution[a];
                if (i_pa < 0 || i_a < 0) {  // removed from solution
                    // moves_to_remove[moveidx] = 1;
                    moves_to_remove.insert(move);
                }
                if (i_a - i_pa != 1 || (i_a != 0 && i_pa != solution_size - 1)) {
                    moves_to_remove.insert(move);
                }

                // edge 2

                auto ab = move.after_b;
                auto b = move.pos2;

                auto i_ab = node_id_to_solution[ab];
                auto i_b = node_id_to_solution[b];
                if (i_ab < 0 || i_b < 0) {  // removed from solution
                    // moves_to_remove[moveidx] = 1;
                    moves_to_remove.insert(move);
                }
                if (i_ab - i_b != 1 || (i_ab != 0 && i_b != solution_size - 1)) {
                    moves_to_remove.insert(move);
                }
            }
        }
        auto _end =
            std::remove_if(improving_moves_list.begin(), improving_moves_list.end(),
                           [&](auto value) { return moves_to_remove.contains(value); }

            );

        {  // Do intra moves
            // for (int i1 = 0; i1 < solution_size; i1++) {
            for (auto i1_id : affected_node_ids) {
                int i1 = node_id_to_solution[i1_id];
                // assert(i1 >= 0);
                if (i1 == -1) continue;
                for (int i2 = 0; i2 < i1; i2++) {
                    // assert(i1 > i2);
                    // if (edges) {
                    if (i2 == 0 && i1 == solution_size - 1)
                        continue;  // skip (0,99) - considering single edge
                    int delta = intraSwapTwoEdgesImpact(solution, distanceMatrix, i1, i2);

                    LSMove currentmove = {.score_delta = delta,
                                          .is_edge_swap = true,
                                          .pos1 = solution[i2],
                                          .pos2 = solution[i1],
                                          .before_a = cyclePrev(solution, i2),
                                          .after_a = cycleNext(solution, i2),
                                          .before_b = cyclePrev(solution, i1),
                                          .after_b = cycleNext(solution, i1)

                    };
                    if (delta < 0) {
                        improving_moves_list.push_back(currentmove);
                    }
                    // if (delta < lowest_delta) {
                    //     found = true;
                    //     best_move = currentmove;
                    //     lowest_delta = delta;
                    // }
                }
            }
        }
        {  // Do inter-route moves
            // for (int internal : in_sol) {
            for (int internal : affected_node_ids) {
                int idx = node_id_to_solution[internal];
                if (idx < 0) continue;
                for (int external : not_in_sol) {
                    int delta = interSwapTwoNodesImpact(solution, distanceMatrix, costs,
                                                        idx, external);
                    LSMove currentmove = {.score_delta = delta,
                                          .is_edge_swap = false,
                                          .pos1 = internal,
                                          .pos2 = external,

                                          .before_a = cyclePrev(solution, idx),
                                          .after_a = cycleNext(solution, idx),
                                          .before_b = -1,
                                          .after_b = -1};
                    if (delta < 0) {
                        improving_moves_list.push_back(currentmove);
                    }
                    // if (delta < lowest_delta) {
                    //     found = true;
                    //     best_move = currentmove;
                    //     lowest_delta = delta;
                    // }
                }
            }
        }

        // std::sort(improving_moves_list.begin(), improving_moves_list.end());
        // best_move = improving_moves_list.front();
        // if (lowest_delta != best_move.score_delta) {
        //     lowest_delta = best_move.score_delta;
        // }
        found = false;
        best_move.score_delta = 0;
        lowest_delta = 0;
        for (auto move : improving_moves_list) {
            // check if move is valid
            if (move.is_edge_swap == false) {
                auto internal = move.pos1;
                auto external = move.pos2;

                auto ii = node_id_to_solution[internal];
                auto ei = node_id_to_solution[external];
                if (ii < 0 || ei >= 0) continue;

                auto i_bef = cycleIndexBefore(solution, ii);
                auto i_aft = cycleIndexAfter(solution, ii);

                if (solution[i_bef] != move.before_a || solution[i_aft] != move.after_a)
                    continue;

                int delta = interSwapTwoNodesImpact(solution, distanceMatrix, costs, ii,
                                                    external);
                if (delta < 0 && delta < lowest_delta) {
                    found = true;
                    lowest_delta = delta;
                    best_move = {.score_delta = delta,
                                 .is_edge_swap = false,
                                 .pos1 = internal,
                                 .pos2 = external,

                                 .before_a = cyclePrev(solution, ii),
                                 .after_a = cycleNext(solution, ii),
                                 .before_b = -1,
                                 .after_b = -1};
                }
            }
            else {
                // edge 1
                auto pa = move.before_a;
                auto a = move.pos1;

                auto i_pa = node_id_to_solution[pa];
                auto i_a = node_id_to_solution[a];
                if (i_pa < 0 || i_a < 0) {
                    continue;
                }
                // if (i_a - i_pa != 1 || (i_a != 0 && i_pa != solution_size - 1)) {
                //     moves_to_remove.insert(move);
                // }

                auto ab = move.after_b;
                auto b = move.pos2;

                auto i_ab = node_id_to_solution[ab];
                auto i_b = node_id_to_solution[b];
                if (i_ab < 0 || i_b < 0) {  // removed from solution
                    continue;
                }
                // if (i_ab - i_b != 1 || (i_ab != 0 && i_b != solution_size - 1)) {
                //     moves_to_remove.insert(move);
                // }

                int i1 = node_id_to_solution[a];
                int i2 = node_id_to_solution[b];
                if (i1 == 0 && i2 == solution_size - 1) continue;
                if (i2 == 0 && i1 == solution_size - 1) continue;
                if (i2 > i1) swap(i1, i2);
                int delta = intraSwapTwoEdgesImpact(solution, distanceMatrix, i1, i2);
                if (delta < 0 && delta < lowest_delta) {
                    found = true;
                    lowest_delta = delta;
                    best_move = {.score_delta = delta,
                                 .is_edge_swap = true,
                                 .pos1 = solution[i2],
                                 .pos2 = solution[i1],
                                 .before_a = cyclePrev(solution, i2),
                                 .after_a = cycleNext(solution, i2),
                                 .before_b = cyclePrev(solution, i1),
                                 .after_b = cycleNext(solution, i1)

                    };
                }
            }
        }
    } while (found);

    return solution;
}
