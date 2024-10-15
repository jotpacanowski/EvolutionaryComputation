#include <algorithm>
#include <cmath>
// #include <filesystem>
#include <fstream>
#include <iostream>
#include <iterator>
#include <numeric>
#include <random>
#include <sstream>  // std::stringstream
#include <string>
#include <vector>

// instead of `using namespace`
using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

// namespace fs = std::filesystem;

vector<vector<int>> readInstance(const char* filename)
{
    string path = "data/instances/";
    path.append(filename);
    std::ifstream data(path);
    vector<vector<int>> instance;
    string line;
    while (std::getline(data, line)) {
        std::stringstream lineStream(line);
        vector<int> row;
        string cell;
        while (std::getline(lineStream, cell, ';')) {
            row.push_back(std::stoi(cell));
        }
        instance.push_back(row);
    }
    return instance;
}

vector<vector<int>> distanceMatrix(const vector<vector<int>>& instance)
{
    const int N = instance.size();

    vector<vector<int>> D(N, vector<int>(N, 0));

    for (int i = 0; i < N; i++) {
        for (int j = 0; j < N; j++) {
            int dist = (int)nearbyint(std::hypot(instance[i][0] - instance[j][0],
                                                 instance[i][1] - instance[j][1]));
            D[i][j] = dist;
        }
    }
    return D;
}

vector<int> getCosts(const vector<vector<int>>& instance)
{
    vector<int> costs(instance.size(), 0);
    for (int i = 0; i < instance.size(); i++) {
        costs[i] = instance[i][2];
    }
    return costs;
}

int evaluateSolution(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
                     const vector<int>& solution)
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
