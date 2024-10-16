#include <cmath>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>

using std::cerr;
using std::cout;
using std::endl;
using std::string;
using std::vector;

struct TSPInstance {
    vector<vector<int>> instance;
    vector<vector<int>> distances;
    vector<int> costs;

    static TSPInstance readFromFile(const std::string& filename)
    {
        string path = "data/instances/" + filename;
        std::ifstream data(path);
        if (!data.is_open()) {
            cerr << "Open " << path << " failed: ";
            perror("");
            exit(1);
        }
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
        if (instance.size() < 2) {
            cerr << "[warn] Read " << instance.size() << " rows\n";
            exit(1);
        }
        cerr << "[ok] Read " << filename << endl;
        return TSPInstance(instance);
    }

    explicit TSPInstance(const vector<vector<int>>& instanceData) : instance(instanceData)
    {
        distances = calculateDistanceMatrix();
        costs = extractCosts();
    }

    vector<vector<int>> calculateDistanceMatrix() const
    {
        const int N = instance.size();
        vector<vector<int>> D(N, vector<int>(N, 0));

        for (int i = 0; i < N; i++) {
            for (int j = 0; j < N; j++) {
                if (i == j) continue;

                double hyp = std::hypot(instance[i][0] - instance[j][0],
                                        instance[i][1] - instance[j][1]);
                D[i][j] = std::lround(hyp);
                // https://en.cppreference.com/w/cpp/numeric/math/round
                // rounding halfway cases away from zero,
                // regardless of the current rounding mode
            }
        }
        return D;
    }

    vector<int> extractCosts() const
    {
        vector<int> costs(instance.size(), 0);
        for (int i = 0; i < instance.size(); i++) {
            costs[i] = instance[i][2];
        }
        return costs;
    }

    int evaluateSolution(const vector<int>& solution) const
    {
        int result = costs[solution[0]];
        for (int i = 1; i < solution.size(); i++) {
            int node = solution[i];
            int prev_node = solution[i - 1];
            result += distances[prev_node][node] + costs[node];
        }
        result += distances[solution.front()][solution.back()];
        return result;
    }
};
