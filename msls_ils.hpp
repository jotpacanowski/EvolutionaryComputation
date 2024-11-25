#pragma once

#include "headers.hpp"

vector<int> generate_random_solution_sliding_window(const vector<vector<int>>& D,
                                                    const vector<int>& C, int seed);

vector<int> multiple_start_steepestLS(const vector<vector<int>>& distanceMatrix,
                                      const vector<int>& costs, int seed);

vector<int> iterative_steepest_LS(const vector<vector<int>>& distanceMatrix,
                                  const vector<int>& costs, int seed);

vector<int> iterative_steepest_LS2(const vector<vector<int>>& distanceMatrix,
                                   const vector<int>& costs, int seed);