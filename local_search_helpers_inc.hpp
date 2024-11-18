#pragma once

#include <vector>

using namespace std;

// Return the previous node of solution[index], wrapping around
constexpr int cyclePrev(const vector<int> &solution, int index)
{
    if (index == 0) return solution[solution.size() - 1];
    return solution[index - 1];
}

// Return the next node of solution[index], wrapping around
constexpr int cycleNext(const vector<int> &solution, int index)
{
    if (index == solution.size() - 1) return solution[0];
    return solution[index + 1];
}

// Subtract 1 from index, wrapping around solution size
constexpr int cycleIndexBefore(const vector<int> &solution, int index)
{
    if (index == 0) return solution.size() - 1;
    return index - 1;
}

// Add 1 to index, wrapping around solution size
constexpr int cycleIndexAfter(const vector<int> &solution, int index)
{
    if (index == solution.size() - 1) return 0;
    return index + 1;
}

static inline int findIndex(const vector<int> &arr, int item)
{
    auto ret = std::find(arr.begin(), arr.end(), item);
    if (ret == arr.end()) return -1;
    return ret - arr.begin();
}
