#include "headers.hpp"

vector<int> randomTSP(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
                      int seed)
{
    const int N = distanceMatrix.size();
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
    solution.resize((N + 1) / 2);
    // solution.shrink_to_fit();
    return solution;
}

static int findNearestNeighbor(const vector<vector<int>>& distanceMatrix,
                               const vector<int>& costs, int node,
                               const vector<int>& solution,
                               const vector<uint8_t>& is_in_sol)
{
    const int N = distanceMatrix.size();

    int current_impact = LARGE_SCORE;
    int nearest_neighbor = -1;

    for (int i = 0; i < N; i++) {
        if (is_in_sol[i] != 0) continue;
        int impact = distanceMatrix[node][i] + costs[i];

        if (impact < current_impact) {
            nearest_neighbor = i;
            current_impact = impact;
        }
    }

    // check for correctness (debug assertion?)
    if (nearest_neighbor == -1) {
        cerr << "error: nearest_neighbor not found" << endl;
        exit(1);
    }
    return nearest_neighbor;
}

vector<int> nearestNeighborTSP(const vector<vector<int>>& distanceMatrix,
                               const vector<int>& costs, int starting_node)
{
    const int N = distanceMatrix.size();
    vector<int> solution;
    solution.reserve((N + 1) / 2);
    vector<uint8_t> is_in_sol(N, 0);
    if (starting_node < 0) {
        starting_node = rand() % N;
    }
    solution.push_back(starting_node);
    is_in_sol[starting_node] = 1;

    int last = starting_node;
    while (solution.size() < ((N + 1) / 2)) {
        int next = findNearestNeighbor(distanceMatrix, costs, last, solution, is_in_sol);
        solution.push_back(next);
        is_in_sol[next] = 1;
        last = next;
    }
    return solution;
}

vector<int> nearestNeighborAnyTSP(const vector<vector<int>>& distanceMatrix,
                                  const vector<int>& costs, int starting_node)
{
    const int N = distanceMatrix.size();
    vector<int> solution;
    solution.reserve((N + 1) / 2);
    vector<uint8_t> is_in_sol(N, 0);
    if (starting_node < 0) {
        starting_node = rand() % N;
    }
    solution.push_back(starting_node);
    is_in_sol[starting_node] = 1;

    while (solution.size() < ((N + 1) / 2)) {
        int best_candidate = 0;
        int best_impact = LARGE_SCORE;
        int candidate_index = 0;
        for (int position = 0; position < solution.size(); position++) {
            int candidate = findNearestNeighbor(distanceMatrix, costs, solution[position],
                                                solution, is_in_sol);
            int impact = distanceMatrix[solution[position]][candidate] + costs[candidate];
            if (impact < best_impact) {
                best_candidate = candidate;
                best_impact = impact;
                candidate_index = position;
            }
        }
        solution.insert(solution.begin() + candidate_index, best_candidate);
        is_in_sol[best_candidate] = 1;
    }
    return solution;
}

vector<int> greedyCycleTSP(const vector<vector<int>>& distanceMatrix,
                           const vector<int>& costs, int starting_node)
{
    const int N = distanceMatrix.size();
    vector<int> solution;
    solution.reserve((N + 1) / 2);
    vector<uint8_t> is_in_sol(N, 0);
    if (starting_node < 0) {
        starting_node = rand() % N;
    }
    solution.push_back(starting_node);
    is_in_sol[starting_node] = 1;

    int candidate =
        findNearestNeighbor(distanceMatrix, costs, starting_node, solution, is_in_sol);
    solution.push_back(candidate);
    is_in_sol[candidate] = 1;

    while (solution.size() < ((N + 1) / 2)) {
        int best_candidate = 0;
        int best_impact = LARGE_SCORE;
        int candidate_index = 0;
        int impact;
        for (int candidate = 0; candidate < N; candidate++)  // For each candidate
        {
            // If candidate already in solution, skip
            if (is_in_sol[candidate]) continue;
            for (int j = 0; j < solution.size(); j++)
            // Check insertion for candidate at each index
            {
                if (j == 0) {
                    // cout <<"INSERT AT 0, BETWEEEN: "<<  solution[j] << ' ' <<
                    // candidate << solution[solution.size() - 1] << endl;

                    impact = distanceMatrix[solution[j]][candidate] +
                             distanceMatrix[solution.size() - 1][candidate] +
                             costs[candidate];
                }
                else {
                    // cout << solution[j - 1] << ' ' << candidate << solution[j] <<
                    // endl;

                    impact = distanceMatrix[solution[j]][candidate] +
                             distanceMatrix[solution[j - 1]][candidate] +
                             costs[candidate];
                }
                // cout<<"IMPACT: " << impact<<endl;
                if (impact < best_impact) {
                    best_impact = impact;
                    best_candidate = candidate;
                    // cout << "Better impact: " << best_impact << endl;

                    candidate_index = j;
                }
            }
        }

        solution.insert(solution.begin() + candidate_index, best_candidate);
        is_in_sol[best_candidate] = 1;
    }
    return solution;
}
