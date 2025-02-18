#pragma once

#include "headers.hpp"

vector<int> greedyCycleRepair(vector<int> vec, const vector<vector<int>>& D,
                              const vector<int>& C, int desired_length);

vector<int> generate_random_solution_sliding_window(const vector<vector<int>>& D,
                                                    const vector<int>& C, int seed);

vector<int> multiple_start_steepestLS(const vector<vector<int>>& distanceMatrix,
                                      const vector<int>& costs, int seed);

vector<int> iterative_steepest_LS(const vector<vector<int>>& distanceMatrix,
                                  const vector<int>& costs, int seed);

vector<int> large_scale_neighbourhood_LS(const vector<vector<int>>& distanceMatrix,
                                         const vector<int>& costs, int seed);

vector<int> large_scale_neighbourhood_LS2(const vector<vector<int>>& distanceMatrix,
                                          const vector<int>& costs, int seed);

// ILS perturbation
vector<int> perturb(const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
                    vector<int>& solution, int seed);

// LNS perturbation
vector<int> destroy_repair_non_deterministic_multiple_chains(
    const vector<vector<int>>& distanceMatrix, const vector<int>& costs,
    vector<int>& solution, int seed);
