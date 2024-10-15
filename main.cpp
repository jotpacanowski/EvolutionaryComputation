#include <algorithm>
#include <cmath>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>  // std::stringstream
#include <string>
#include <vector>

std::vector<std::vector<int>> readInstance(char *filename)
{
    std::string path = "data/instances/";
    path.append(filename);
    std::ifstream data(path);
    std::vector<std::vector<int>> instance;
    std::string line;
    while (std::getline(data, line)) {
        std::stringstream lineStream(line);
        std::vector<int> row;
        std::string cell;
        while (std::getline(lineStream, cell, ';')) {
            row.push_back(std::stoi(cell));
        }
        instance.push_back(row);
    }
    return instance;
}

std::vector<std::vector<int>> distanceMatrix(std::vector<std::vector<int>> instance)
{
    int instance_size = instance.size();

    std::vector<std::vector<int>> D(instance_size, std::vector<int>(instance_size, 0));

    for (int i = 0; i < instance_size; i++) {
        for (int j = 0; j < instance_size; j++) {
            int dist = (int)nearbyint(std::hypot(instance[i][0] - instance[j][0],
                                                 instance[i][1] - instance[j][1]));
            D[i][j] = dist;
        }
    }
    return D;
}

std::vector<int> getCosts(std::vector<std::vector<int>> instance)
{
    std::vector<int> costs(instance.size(), 0);
    for (int i = 0; i < instance.size(); i++) {
        costs[i] = instance[i][2];
    }
    return costs;
}

int evaluateSolution(std::vector<std::vector<int>> distanceMatrix, std::vector<int> costs,
                     std::vector<int> solution)
{
    int result = costs[solution[0]];
    for (int i = 1; i < solution.size(); i++) {
        int node = solution[i];
        int prev_node = solution[i - 1];
        result += distanceMatrix[prev_node][node] + costs[node];
    }
    result += distanceMatrix[solution[0]][solution[solution.size() - 1]];
    return result;
}

std::vector<int> randomTSP(std::vector<std::vector<int>> distanceMatrix,
                           std::vector<int> costs)
{
    std::vector<int> solution;
    for (int i = 0; i < distanceMatrix.size(); i++) {
        solution.push_back(i);
    }
    std::random_device rd;
    std::mt19937 g(rd());

    std::shuffle(solution.begin(), solution.end(), g);
    return solution;
}

int nearestNeighbor(std::vector<std::vector<int>> distanceMatrix, std::vector<int> costs,
                    int node, std::vector<int> solution)
{
    std::vector<int> neighbors = distanceMatrix[node];
    int current_impact = 1000000;
    int nearest_neighbor = 0;
    int impact;
    for (int i = 0; i < neighbors.size(); i++) {
        // std::cout << "i: " << i << std::endl;
        if (std::count(solution.begin(), solution.end(), i)) {
            // std::cout << "in solution: " << i << std::endl;
            continue;
        }
        impact = distanceMatrix[node][i] + costs[i];
        // std::cout << impact << std::endl;
        if (impact < current_impact) {
            nearest_neighbor = i;
            current_impact = impact;
        }
    }
    // std::cout << nearest_neighbor << std::endl;
    return nearest_neighbor;
}

std::vector<int> nearestNeighborTSP(std::vector<std::vector<int>> distanceMatrix,
                                    std::vector<int> costs, int starting_node)
{
    std::vector<int> solution;
    if (starting_node == NULL) {
        starting_node = rand() % 200;
    }
    solution.push_back(starting_node);
    int last = starting_node;
    int next;
    for (int i = 0; i < distanceMatrix.size() / 2; i++) {
        next = nearestNeighbor(distanceMatrix, costs, last, solution);
        solution.push_back(next);
        last = next;
    }
    return solution;
}

std::vector<int> nearestNeighborAnyTSP(std::vector<std::vector<int>> distanceMatrix,
                                       std::vector<int> costs, int starting_node)
{
    std::vector<int> solution;
    if (starting_node == NULL) {
        starting_node = rand() % 200;
    }
    solution.push_back(starting_node);
    int last = starting_node;
    int next;
    for (int i = 0; i < distanceMatrix.size() / 2; i++) {
        // std::vector<int> candidates;
        int best_candidate = 0;
        int best_impact = 1000000;
        int candidate_index = 0;
        for (int j = 0; j < solution.size(); j++) {
            int candidate = nearestNeighbor(distanceMatrix, costs, solution[j], solution);
            int impact = distanceMatrix[solution[j]][candidate] + costs[candidate];
            if (impact < best_impact) {
                best_candidate = candidate;
                best_impact = impact;
                candidate_index = j;
            }
        }
        next = best_candidate;
        solution.insert(solution.begin() + candidate_index, best_candidate);
        last = next;
    }
    return solution;
}

void printResults(std::vector<int> sol, bool toFile = false,
                  char *filename = "solution.txt")
{
    if (toFile) {
        std::string path = "results/";
        path.append(filename);
        std::ofstream out(path);
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(out, "\n"));
    }
    else {
        std::copy(sol.begin(), sol.end(), std::ostream_iterator<int>(std::cout, "\n"));
    }
}

void randomTSPTask(std::vector<std::vector<int>> distanceMatrix, std::vector<int> costs)
{
    std::vector<int> sol;
    int value;

    std::vector<int> best_sol;
    int best_sol_value = 1000000;

    std::vector<int> worst_sol;
    int worst_sol_value = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        sol = randomTSP(distanceMatrix, costs);
        value = evaluateSolution(distanceMatrix, costs, sol);
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
    }

    std::cout << "BEST: " << best_sol_value << std::endl;
    printResults(best_sol, true, "random/best.txt");
    std::cout << "WORST: " << worst_sol_value << std::endl;
    printResults(worst_sol, true, "random/worst.txt");
}

void nearestNeighborTSPTask(std::vector<std::vector<int>> distanceMatrix,
                            std::vector<int> costs)
{
    std::vector<int> sol;
    int value;

    std::vector<int> best_sol;
    int best_sol_value = 1000000;

    std::vector<int> worst_sol;
    int worst_sol_value = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        sol = nearestNeighborTSP(distanceMatrix, costs, i);
        value = evaluateSolution(distanceMatrix, costs, sol);
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
    }

    std::cout << "BEST: " << best_sol_value << std::endl;
    printResults(best_sol, true, "nearest/best.txt");
    std::cout << "WORST: " << worst_sol_value << std::endl;
    printResults(worst_sol, true, "nearest/worst.txt");
}

void nearestNeighborAnyTSPTask(std::vector<std::vector<int>> distanceMatrix,
                               std::vector<int> costs)
{
    std::vector<int> sol;
    int value;

    std::vector<int> best_sol;
    int best_sol_value = 1000000;

    std::vector<int> worst_sol;
    int worst_sol_value = 0;

    // Generating 200 random solutions, keeping the best and the worst ones.
    for (int i = 0; i < 200; i++) {
        sol = nearestNeighborAnyTSP(distanceMatrix, costs, i);
        value = evaluateSolution(distanceMatrix, costs, sol);
        if (value < best_sol_value) {
            best_sol_value = value;
            best_sol = sol;
        }
        if (value > worst_sol_value) {
            worst_sol_value = value;
            worst_sol = sol;
        }
    }

    std::cout << "BEST: " << best_sol_value << std::endl;
    printResults(best_sol, true, "nearest_any/best.txt");
    std::cout << "WORST: " << worst_sol_value << std::endl;
    printResults(worst_sol, true, "nearest_any/worst.txt");
}

int main()
{
    std::vector<std::vector<int>> instance = readInstance("TSPA.csv");
    std::vector<std::vector<int>> D = distanceMatrix(instance);
    std::vector<int> costs = getCosts(instance);
    // for (int i = 0; i < D.size(); i++)
    // {
    //     for (int j = 0; j < D.size(); j++)
    //     {
    //         std::cout << D[i][j] << ' ';
    //     }
    //     std::cout << std::endl;
    // }
    randomTSPTask(D, costs);
    nearestNeighborTSPTask(D, costs);
    nearestNeighborAnyTSPTask(D, costs);
    // std::vector<int> sol = nearestNeighborAnyTSP(D, costs, 0);
    // // std::vector<int>
    // //     sol = nearestNeighborTSP(D, costs, 0);
    // printResults(sol);
    // std::cout << evaluateSolution(D, costs, sol) << std::endl;

    return 0;
}
