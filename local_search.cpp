#include <algorithm>
#include <cstdint>
#include <vector>

#include "headers.hpp"

int intraSwapTwoNodesImpact(const vector<int> &solution, int id1, int id2,
                            const vector<vector<int>> &distanceMatrix,
                            const int solution_size)
{
    int node1 = solution[id1];
    int node2 = solution[id2];

    int delta;

    int prev1 = id1 == 0 ? solution[(solution_size - 1)] : solution[id1 - 1];
    int next1 = id1 == (solution_size - 1) ? solution[0] : solution[id1 + 1];
    int prev2 = id2 == 0 ? solution[(solution_size - 1)] : solution[id2 - 1];
    int next2 = id2 == (solution_size - 1) ? solution[0] : solution[id2 + 1];

    if (next1 == node2) {  // id1+1=id2
        delta = -distanceMatrix[prev1][node1] + distanceMatrix[prev1][node2] -
                distanceMatrix[node2][next2] + distanceMatrix[node1][next2];
        return delta;
    }
    if (next2 == node1) {
        delta = -distanceMatrix[node1][next1] + distanceMatrix[node2][next1] -
                distanceMatrix[prev2][node2] + distanceMatrix[prev2][node1];
        return delta;
    }
    delta = -distanceMatrix[prev1][node1] + distanceMatrix[prev1][node2] -
            distanceMatrix[node1][next1] + distanceMatrix[node2][next1] -
            distanceMatrix[prev2][node2] + distanceMatrix[prev2][node1] -
            distanceMatrix[node2][next2] + distanceMatrix[node1][next2];
    return delta;
}

int intraSwapTwoEdgesImpact(const vector<int> &solution, int id1, int id2,
                            const vector<vector<int>> &distanceMatrix,
                            const int solution_size)
{
    int node1 = solution[id1];
    int node2 = solution[id2];
    int delta, prev1, prev2, next1, next2;

    if (id1 < id2) {  //
        prev1 = id1 == 0 ? solution[(solution_size - 1)] : solution[id1 - 1];
        next2 = id2 == (solution_size - 1) ? solution[0] : solution[id2 + 1];

        delta = -distanceMatrix[prev1][node1] + distanceMatrix[prev1][node2] -
                distanceMatrix[node2][next2] + distanceMatrix[node1][next2];
    }
    else {  //(id1 > id2)
        next1 = id1 == (solution_size - 1) ? solution[0] : solution[id1 + 1];
        prev2 = id2 == 0 ? solution[(solution_size - 1)] : solution[id2 - 1];

        delta = -distanceMatrix[prev2][node2] + distanceMatrix[prev2][node1] -
                distanceMatrix[node1][next1] + distanceMatrix[node2][next1];
    }
    return delta;
}

int interSwapTwoNodesImpact(const vector<int> &solution, int idx,
                            const vector<vector<int>> &distanceMatrix,
                            const int solution_size, const vector<uint8_t> &is_in_sol)
{
}

vector<int> localSearch(vector<int> solution, const vector<vector<int>> &distanceMatrix,
                        const vector<int> &costs)
{
    int variable = 0;
    int delta;
    int highest_delta = 0;
    int pos1;
    int pos2;
    bool found = false;
    vector<uint8_t> is_in_sol(costs.size(), 0);
    for (const int node : solution) {
        is_in_sol[node] = 1;
    }
    cout << "SIZE: " << solution.size() << endl;

    //     variable = rand() % 2;
    for (int _ = 0; _ < 1000; _++) {
        // if (variable == 0) {  // Do only intra moves
            for (int i1 = 0; i1 < solution.size(); i1++) {
                for (int i2 = 0; i2 < solution.size(); i2++) {
                    delta = intraSwapTwoNodesImpact(
                        solution, i1, i2, distanceMatrix,
                        solution.size());  // Swap two nodes variation
                    // delta = intraSwapTwoEdgesImpact(
                    //     solution, i1, i2, distanceMatrix,
                    //     solution.size());  // Swap two edges variation

                    delta = -delta;               // WE WANT IMPROVEMENT -> smaller score
                    if (delta > highest_delta) {  // Steepest variation
                        cout << "NEW HIGHEST DELTA: " << delta << endl;
                        cout << i1 << '\t' << i2 << endl;
                        highest_delta = delta;
                        pos1 = i1;
                        pos2 = i2;
                        found = true;
                    }
                }
            }
            cout << "FOUND: " << found << endl;
            if (found) {
                found = false;
                // reverse(solution.begin()+pos1,solution.begin()+pos2+1);
                iter_swap(solution.begin() + pos1, solution.begin() + pos2);
                continue;
            }
            else {
                cout << "No further improvement found" << endl;
                return solution;
            }
        // }
        // else {
        // }
    }
    cout << "HUH???" << endl;
    return solution;
}