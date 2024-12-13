#include "convexity.hpp"

#include <algorithm>
#include <cerrno>
#include <cmath>
#include <cstdio>
#include <fstream>
#include <iterator>
#include <set>
#include <string_view>
#include <utility>
#include <vector>

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

int similarity_common_edges(vector<int> local_optimum, vector<int> best_solution)
{
    set<pair<int, int>> v1;
    set<pair<int, int>> v2;
    for (int i = 0; i < local_optimum.size() - 1; i++) {
        v1.insert(make_pair(local_optimum[i], local_optimum[i + 1]));
        v1.insert(make_pair(local_optimum[i + 1], local_optimum[i]));

        v2.insert(make_pair(best_solution[i + 1], best_solution[i]));
        v2.insert(make_pair(best_solution[i], best_solution[i + 1]));
    }
    v1.insert(make_pair(local_optimum[local_optimum.size() - 1], local_optimum[0]));
    v1.insert(make_pair(local_optimum[0], local_optimum[local_optimum.size() - 1]));

    v2.insert(make_pair(best_solution[best_solution.size() - 1], best_solution[0]));
    v2.insert(make_pair(best_solution[0], best_solution[best_solution.size() - 1]));

    vector<pair<int, int>> v3;
    // cout << "V1 size " << v1.size() << endl;
    // cout << "V2 size " << v2.size() << endl;
    ranges::set_intersection(v1, v2, std::back_inserter(v3));
    // for (auto c : v3) {
    //     cout << "(" << c.first << ", " << c.second << ")" << endl;
    // }
    return v3.size() / 2;
}

int similarity_common_nodes(vector<int> local_optimum, vector<int> best_solution)
{
    set<int> v1(local_optimum.begin(), local_optimum.end());
    set<int> v2(best_solution.begin(), best_solution.end());

    vector<int> v3;
    ranges::set_intersection(v1, v2, std::back_inserter(v3));
    return v3.size();
}

void task8(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
           const vector<int>& best_solution, string_view instance)
{
    // vector<int> a = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
    // vector<int> b = {10, 9, 8, 7, 6, 5, 4, 3, 2, 1};

    // cout << "NODES" << similarity_common_nodes(a, b) << endl;
    // cout << "EDGES" << similarity_common_edges(a, b) << endl;
    // return;
    vector<float> all_edges;
    vector<float> all_nodes;
    vector<vector<int>> all_local_optima;
    vector<float> all_eval;
    all_eval.push_back(_evaluate_solution(best_solution, distanceMatrix, costs));
    for (int i = 0; i < 1000; i++) {
        vector<int> local_optimum = greedyLocalSearch(
            randomTSP(distanceMatrix, costs, i + 42), distanceMatrix, costs, true);
        int similarity_edges = similarity_common_edges(local_optimum, best_solution);
        int similarity_nodes = similarity_common_nodes(local_optimum, best_solution);
        all_edges.push_back(similarity_edges);
        all_nodes.push_back(similarity_nodes);
        all_local_optima.push_back(local_optimum);
        all_eval.push_back(_evaluate_solution(local_optimum, distanceMatrix, costs));

        // cout << similarity_edges << '\t' << similarity_nodes << endl;
    }
    // cout << "edges\tnodes\t\n";
    // cout << sum(all_nodes / 1000.0 << '\t' << (float)all_edges / 1000.0 << endl;

    vector<float> all_intra_nodes;
    vector<float> all_intra_edges;
    for (int i = 0; i < all_local_optima.size(); i++) {
        float avg_intra_nodes = 0;
        float avg_intra_edges = 0;
        float c = 0;
        for (int j = 0; j < all_local_optima.size(); j++) {
            if (i == j) {
                continue;
            }
            c++;
            avg_intra_edges +=
                similarity_common_edges(all_local_optima[i], all_local_optima[j]);
            avg_intra_nodes +=
                similarity_common_nodes(all_local_optima[i], all_local_optima[j]);
        }
        all_intra_nodes.push_back((float)avg_intra_nodes / c);
        all_intra_edges.push_back((float)avg_intra_edges / c);
    }

    vector<string_view> s1 = {"best_edges", "best_nodes", "intra_edges", "intra_nodes",
                              "solutions"};
    vector<vector<float>> s2 = {all_edges, all_nodes, all_intra_edges, all_intra_nodes,
                                all_eval};
    for (int i = 0; i < s1.size(); i++) {
        string_view prefix = "data/convexity";
        ofstream out = openOutFile(std::format("{}_{}.txt", instance, s1[i]), prefix);
        std::copy(s2[i].begin(), s2[i].end(), std::ostream_iterator<float>(out, "\n"));
    }

    // cout << "Intra edges\tIntra nodes\n";
    // cout << (float)all_intra_edges / c << '\t' << (float)all_intra_nodes / c << endl;
}