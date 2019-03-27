#pragma once

#include <vector>

// returns a list of cores for which the calling thread has affinity
std::vector<int> get_affinity();

void set_affinity(const std::vector<int>& cores);
