#include "generate_chains.h"

const ChainMatrixType GenerateChainMatrixFromErrorVector(int L,
                                                         const vector<int>& errors) {
  size_t chain_size = errors.size();
  ChainMatrixType chain(chain_size);
  
  chain[0].resize(L);
  queue<int> zero_indices;
  queue<int> one_indices;
  for (int index = 0; index < errors[0]; ++index) {
    one_indices.push(index);
    chain[0][index] = 1;
  }
  for (int index = errors[0]; index < L; ++index) {
    zero_indices.push(index);
    chain[0][index] = 0;
  }
  int prev_algor_error = errors[0];
  for (int index = 1; index < chain_size; ++index) {
    chain[index] = chain[index - 1];
    if (errors[index] < prev_algor_error) {
      int one_index = one_indices.front();
      chain[index][one_index] = 0;
      one_indices.pop();
      zero_indices.push(one_index);
    }
    else {
      int zero_index = zero_indices.front();
      chain[index][zero_index] = 1;
      zero_indices.pop();
      one_indices.push(zero_index);
    }
    prev_algor_error = errors[index];
  }
  return chain;
}


/* проверить */
vector<int> generateClassesFromErrorVector(const vector<int>& chain_error_vector) {
  size_t chain_size = chain_error_vector.size();
  vector<int> classes(chain_size - 1, -1);
  for (int i = 0; i < chain_size - 1; ++i) {
    if (chain_error_vector[i + 1] > chain_error_vector[i]) {
      classes[i] = 1;
    }
  }
  return classes;
}



