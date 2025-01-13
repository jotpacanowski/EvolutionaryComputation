#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstdlib>
#include <numeric>
#include <set>
#include <type_traits>
#include <unordered_set>
#include <utility>
#include <vector>

#include "local_search.hpp"
#include "msls_ils.hpp"
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

vector<int> getCommonNodes(vector<int> v1, vector<int> v2)
{
    set<int> s1(v1.begin(), v1.end());
    set<int> s2(v2.begin(), v2.end());

    vector<int> v3;
    ranges::set_intersection(s1, s2, std::back_inserter(v3));
    return v3;
}

int lookForward(const vector<int>& v, int idx)
{
    if (idx < v.size() - 1) {
        return idx + 1;
    }
    return 0;
}

int lookBack(const vector<int>& v, int idx)
{
    if (idx == 0) {
        return v.size() - 1;
    }
    return idx - 1;
}
vector<vector<int>> getCommonSubsequences(vector<int> v1, vector<int> v2,
                                          vector<int> common_nodes)
{
    vector<vector<int>> common_subsequences;
    if (common_nodes.size() == v1.size()) {
        common_subsequences.push_back(v1);
        return common_subsequences;
    }
    while (common_nodes.size() > 0) {
        int node = common_nodes[0];
        // cout << node << endl;
        vector<int> subseq{node};
        auto it1 = find(v1.begin(), v1.end(), node);
        int idx1 = it1 - v1.begin();

        auto it2 = find(v2.begin(), v2.end(), node);
        int idx2 = it2 - v2.begin();

        bool found_forward = false;
        int idx1_forward = idx1;
        int idx2_forward = idx2;
        vector<int> nodes_forward;

        if (v1[lookForward(v1, idx1)] == v2[lookForward(v2, idx2)]) {
            found_forward = true;
        }

        while (found_forward) {
            idx1_forward = lookForward(v1, idx1_forward);
            idx2_forward = lookForward(v2, idx2_forward);
            // cout<<"Looking forward "<<v1[idx1_forward]<<endl;
            nodes_forward.push_back(v1[idx1_forward]);
            if (nodes_forward.size() == common_nodes.size()) {
                break;
            }
            if (v1[lookForward(v1, idx1_forward)] != v2[lookForward(v2, idx2_forward)]) {
                found_forward = false;
            }
        }

        bool found_back = false;
        int idx1_back = idx1;
        int idx2_back = idx2;
        vector<int> nodes_back;

        if (v1[lookBack(v1, idx1)] == v2[lookBack(v2, idx2)]) {
            found_back = true;
        }

        while (found_back) {
            // cout<<"Looking back"<<endl;
            idx1_back = lookBack(v1, idx1_back);
            idx2_back = lookBack(v2, idx2_back);
            nodes_back.push_back(v1[idx1_back]);
            if (nodes_back.size() == common_nodes.size()) {
                break;
            }
            if (v1[lookBack(v1, idx1_back)] != v2[lookBack(v2, idx2_back)]) {
                found_back = false;
            }
        }

        bool found_forward_back = false;
        idx1_forward = idx1;
        idx2_back = idx2;
        vector<int> nodes_forward_back;

        if (v1[lookForward(v1, idx1)] == v2[lookBack(v2, idx2)]) {
            found_forward_back = true;
        }

        while (found_forward_back) {
            // cout<<"Looking back"<<endl;
            idx1_forward = lookForward(v1, idx1_forward);
            idx2_back = lookBack(v2, idx2_back);
            // cout << "Looking forward " << v1[idx1_forward] << " back " << v2[idx2_back]
            //      << endl;

            nodes_forward_back.push_back(v1[idx1_forward]);
            if (nodes_forward_back.size() == common_nodes.size()) {
                break;
            }
            if (v1[lookForward(v1, idx1_forward)] != v2[lookBack(v2, idx2_back)]) {
                found_forward_back = false;
            }
        }
        bool found_back_forward = false;
        idx1_back = idx1;
        idx2_forward = idx2;
        vector<int> nodes_back_forward;

        if (v1[lookBack(v1, idx1)] == v2[lookForward(v2, idx2)]) {
            found_back_forward = true;
        }

        while (found_back_forward) {
            idx1_back = lookBack(v1, idx1_back);
            idx2_forward = lookForward(v2, idx2_forward);
            // cout << "Looking back " << v1[idx1_back] << " forward " << v2[idx2_forward]
            //      << endl;
            nodes_back_forward.push_back(v1[idx1_back]);
            if (nodes_back_forward.size() == common_nodes.size()) {
                break;
            }
            if (v1[lookBack(v1, idx1_back)] != v2[lookForward(v2, idx2_forward)]) {
                found_back_forward = false;
            }
        }

        subseq.insert(subseq.begin(), nodes_back.rbegin(), nodes_back.rend());
        subseq.insert(subseq.end(), nodes_forward.begin(), nodes_forward.end());
        subseq.insert(subseq.begin(), nodes_back_forward.rbegin(),
                      nodes_back_forward.rend());
        subseq.insert(subseq.end(), nodes_forward_back.begin(), nodes_forward_back.end());

        for (int n : subseq) {
            erase(common_nodes, n);
        }
        // if (subseq.size() == 1) {
        //     cout<<"subseq of lenght 1"<< endl;
        //     // continue;
        // }
        common_subsequences.push_back(subseq);
    }
    return common_subsequences;
}

vector<int> recombination1(const vector<int>& parent1, const vector<int>& parent2)
{
    set<int> ALL_NUMBERS;

    for (int i = 0; i < 200; ++i) {
        ALL_NUMBERS.insert(i);
    }
    vector<int> common_nodes = getCommonNodes(parent1, parent2);

    vector<vector<int>> common_subseq =
        getCommonSubsequences(parent1, parent2, common_nodes);

    vector<int> offspring;
    for (const auto& vec : common_subseq) {
        // Concatenate by inserting each vector's elements into the result
        offspring.insert(offspring.end(), vec.begin(), vec.end());
    }
    set<int> in_offspring(offspring.begin(), offspring.end());

    set<int> possible;
    std::set_difference(ALL_NUMBERS.begin(), ALL_NUMBERS.end(), in_offspring.begin(),
                        in_offspring.end(), std::inserter(possible, possible.end()));
    vector<int> chosen;
    random_device rd;
    mt19937 g(rd());

    sample(possible.begin(), possible.end(), std::back_inserter(chosen),
           parent1.size() - offspring.size(), g);
    shuffle(chosen.begin(), chosen.end(), g);
    offspring.insert(offspring.end(), chosen.begin(), chosen.end());
    return offspring;
}

vector<int> recombination2(const vector<vector<int>>& distanceMatrix,
                           const vector<int>& costs, const vector<int>& parent1,
                           const vector<int>& parent2)
{
    vector<int> common_nodes = getCommonNodes(parent1, parent2);

    vector<vector<int>> common_subseq =
        getCommonSubsequences(parent1, parent2, common_nodes);

    vector<int> offspring;
    for (const auto& vec : common_subseq) {
        // Concatenate by inserting each vector's elements into the result
        offspring.insert(offspring.end(), vec.begin(), vec.end());
    }
    offspring = greedyCycleRepair(offspring, distanceMatrix, costs, parent1.size());
    return offspring;
}

vector<int> recombination3(const vector<vector<int>>& distanceMatrix,
                           const vector<int>& costs, const vector<int>& parent1,
                           const vector<int>& parent2)
{
    set<int> ALL_NUMBERS;

    for (int i = 0; i < 200; ++i) {
        ALL_NUMBERS.insert(i);
    }

    vector<int> common_nodes = getCommonNodes(parent1, parent2);
    vector<int> offspring(parent1.size(), -1);
    vector<vector<int>> common_subseq =
        getCommonSubsequences(parent1, parent2, common_nodes);
    int eh = 0;
    for (const auto& subseq : common_subseq) {
        if (subseq.size() == 1) {
            eh++;
            continue;
        }
        for (int i : subseq) {
            for (int idx = 0; idx < parent1.size(); idx++) {
                if (parent1[idx] == i) {
                    offspring[idx] = parent1[idx];
                }
            }
        }
    }
    set<int> in_offspring(offspring.begin(), offspring.end());
    in_offspring.erase(-1);

    set<int> possible;
    std::set_difference(ALL_NUMBERS.begin(), ALL_NUMBERS.end(), in_offspring.begin(),
                        in_offspring.end(), std::inserter(possible, possible.end()));
    vector<int> chosen;
    random_device rd;
    mt19937 g(rd());
    int new_nodes = parent1.size() - common_nodes.size();
    new_nodes = parent1.size() - (common_nodes.size() - eh);
    // cout << "NEW NODES " << new_nodes;
    sample(possible.begin(), possible.end(), std::back_inserter(chosen), new_nodes, g);
    shuffle(chosen.begin(), chosen.end(), g);
    int a = 0;
    for (int i = 0; i < offspring.size(); i++) {
        if (offspring[i] == -1) {
            offspring[i] = chosen[a];
            a++;
        }
    }
    return offspring;
}

vector<int> recombination4(const vector<vector<int>>& distanceMatrix,
                           const vector<int>& costs, const vector<int>& parent1,
                           const vector<int>& parent2)
{
    set<int> ALL_NUMBERS;

    for (int i = 0; i < 200; ++i) {
        ALL_NUMBERS.insert(i);
    }
    vector<int> common_nodes = getCommonNodes(parent1, parent2);
    vector<int> offspring(parent1.size(), -1);
    vector<vector<int>> common_subseq =
        getCommonSubsequences(parent1, parent2, common_nodes);
    for (const auto& subseq : common_subseq) {
        if (subseq.size() == 1) {
            // eh++;
            continue;
        }
        for (int i : subseq) {
            for (int idx = 0; idx < parent1.size(); idx++) {
                if (parent1[idx] == i) {
                    offspring[idx] = parent1[idx];
                }
            }
        }
    }
    std::erase_if(offspring, [](int x) { return x == -1; });
    // if ((parent1.size() - offspring.size()) > 10) {
    //     cout << "How many nodes change: " << parent1.size() - offspring.size() << endl;
    // }
    offspring = greedyCycleRepair(offspring, distanceMatrix, costs, parent1.size());

    return offspring;
}

std::random_device random_dev;

vector<int> evolutionarySolver(const vector<vector<int>>& distanceMatrix,
                               const vector<int>& costs, int seed,
                               bool localSearch = false)
{
    using namespace std::chrono;
    using namespace std::chrono_literals;

    // std::mt19937_64 rng(seed);
    std::mt19937_64 rng(random_dev());

    // 1. Generate initial population
    const int POPULATION_SIZE = 100;
    static_assert(POPULATION_SIZE > 5, "minimal POPULATION_SIZE");

    vector<vector<int>> population;
    for (int i = 0; i < POPULATION_SIZE; i++) {
        vector<int> p1 = randomTSP(distanceMatrix, costs, (rng() + seed + (9 * (i + 3))));
        p1 = greedyLocalSearch(p1, distanceMatrix, costs, true);
        population.push_back(p1);
    }

    vector<int> population_scores;
    population_scores.reserve(POPULATION_SIZE);
    int val = 0;
    for (int i = 0; i < POPULATION_SIZE; i++) {
        vector<int> p1;
        p1.reserve(100);
        int p1_score = 1000000;
        do {
            p1 = randomTSP(distanceMatrix, costs, (seed + (3 * (i + 3+val))));
            p1 = steepestLocalSearch(p1, distanceMatrix, costs, true);
            p1_score = _evaluate_solution(p1, distanceMatrix, costs);
            if (std::find(population_scores.begin(), population_scores.end(), p1_score)
                != population_scores.end()) {
                cout << "Repeat in population" << endl;
                val++;
                continue;
            }
            break;
        } while (true);
        population_scores.push_back(p1_score);
        population.push_back(p1);
    }

    int worst_idx = std::distance(
        population_scores.begin(),
        std::max_element(population_scores.begin(), population_scores.end()));

    std::uniform_int_distribution<> dist100(0, POPULATION_SIZE - 1);
    std::uniform_int_distribution<> dist99(0, POPULATION_SIZE - 2);

    int patience = 0;
    int a = 0;
    auto now = chrono::steady_clock::now;
    // auto work_duration = 1'350ms;
    auto work_duration_ms = 1'350;
    auto start = now();
    int i = 0;
    int repeats = 0;
    while ( chrono::duration_cast<chrono::milliseconds>(now() - start).count() < work_duration_ms) {
        i++;

        // 2. Select Parents
        // int p1 = dist100(random_dev);
        int p1 = dist100(rng);
        int p2 = dist99(rng);
        if(p2 >= p1) p2 = p2 + 1;
        vector<int> parent1 = population[p1];
        vector<int> parent2 = population[p2];

        // 3. Construct child
        // vector<int> child = recombination1(parent1, parent2);
        vector<int> child;
        child.reserve(parent1.size());
        if (rand() % 2 == 0) {
            // cout<<"v1"<<endl;
            child = recombination4(distanceMatrix, costs, parent1, parent2);
        }
        else {
            // cout<<"v2"<<endl;
            // child = recombination1(parent1, parent2);
            child = recombination3(distanceMatrix, costs, parent1, parent2);
        }
        if (localSearch) {
            if (rand() % 2 == 0) {
                child = greedyLocalSearch(child, distanceMatrix, costs, true);
            }
        }

        // vector<int> child2 = recombination4(distanceMatrix, costs, parent1, parent2);
        // child2 = steepestLocalSearch(child2, distanceMatrix, costs, true);
        // cout << "parent1: " << _evaluate_solution(parent1, distanceMatrix, costs) <<
        // endl;

        // for (auto t : parent1) {
        //     cout << t << ',';
        // }
        // cout << endl;
        // cout << "parent2: " << _evaluate_solution(parent2, distanceMatrix, costs) <<
        // endl; for (auto t : parent2) {
        //     cout << t << ',';
        // }
        // cout << endl;
        // cout << "child1 (insertion random): "
        //      << _evaluate_solution(child, distanceMatrix, costs) << endl;
        // for (auto t : child) {
        //     cout << t << ',';
        // }
        // cout << endl;
        // cout << "child2 (repair): " << _evaluate_solution(child2, distanceMatrix,
        // costs)
        //      << endl;
        // for (auto t : child2) {
        //     cout << t << ',';
        // }
        // cout << endl;
        // return child;

        int child_score = _evaluate_solution(child, distanceMatrix, costs);
        if (std::find(population_scores.begin(), population_scores.end(), child_score)
            != population_scores.end()) {
            repeats++;
            patience++;
            continue;
        }
        // if (i % 200 == 0) {
        //     cout << "P1: " << p1 << " P2: " << p2 << endl;
        // cout << "Parent 1: " << population_scores[p1] << endl;
        //     // for (auto i : parent1) {
        //     //     cout << i << '\t';
        //     // }
        //     // cout << endl;
        // cout << "Parent 2: " << population_scores[p2] << endl;
        //     // for (auto i : parent2) {
        //     //     cout << i << '\t';
        //     // }
        //     // cout << endl;
        // cout << "Child1: " << child_score << endl;
        //     // for (auto i : child) {
        //     //     cout << i << '\t';
        //     // }
        //     // cout << endl;
        // }

        // 4. If child is better than worst element of population, add child to population
        if (child_score < population_scores[worst_idx]) {
            a++;
            patience = 0;
            // cout << "Found improvement" << endl;
            population_scores[worst_idx] = child_score;
            population[worst_idx] = child;
            worst_idx = std::distance(
                population_scores.begin(),
                std::max_element(population_scores.begin(), population_scores.end()));
        }
        patience++;
        // if (patience > 1000) {
        //     cout << "Run out of patience, loop: " << i << endl;
        //     break;
        // }
    }
    cout << "Iterations: " << i << "\tGood children: " << a << "\tRepeats: " << repeats
         << endl;
    // for (int idx = 0; idx < POPULATION_SIZE; idx++) {
    //     population[idx] =
    //         steepestLocalSearch(population[idx], distanceMatrix, costs, true);
    //     population_scores[idx] =
    //         _evaluate_solution(population[idx], distanceMatrix, costs);
    // }
    int best_idx = std::distance(
        population_scores.begin(),
        std::min_element(population_scores.begin(), population_scores.end()));

    // cout << population_scores[best_idx] << endl;
    return population[best_idx];
}