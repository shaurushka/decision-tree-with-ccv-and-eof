#include <vector>
#include <queue>
#include <iostream>
#include <string>
#include <cmath>
#include <cstdlib>
#include "types_definition.h"

using std::vector;
using std::queue;
using std::string;

const ChainMatrixType GenerateChainMatrixFromErrorVector(int L,
                                                         const vector<int>& errors);


struct Chain {
  vector<int> error_vector;
  string type;
};


vector<int> generateClassesFromErrorVector(const vector<int>& chain_error_vector);

