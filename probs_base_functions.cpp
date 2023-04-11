#include "probs_base_functions.h"

vector<int> GetGraphEdges(const vector<vector<int> >& chain) {
  size_t chain_size = chain.size();
  if (chain_size < 1) {
    return vector<int>(0);
  }
  
  vector<int> edges(chain_size - 1);
  for (int index = 1; index < chain_size; ++index) {
    int element_index = 0;
    while (chain[index][element_index] == chain[index - 1][element_index]) {
      ++element_index;
    }
    edges[index - 1] = element_index;
  }
  return edges;
}


int LeftSampleErrorCount(const vector<vector<int> >& chain,
                         int d,
                         const vector<int>& graph_edges) {
  int err_count = 0;
  for (int index = 0; index < d; ++index) {
    err_count += chain[d][graph_edges[index]];
  }
  return err_count;
}


bool PessimististicConditionIsHeld(int first_alg_err, int current_alg_err, const CHAIN_DIRECTION& chain_dir) {
  return ((chain_dir == RIGHT && current_alg_err < first_alg_err) ||
          (chain_dir == LEFT && current_alg_err <= first_alg_err));
}

int GetMaxTrainErrorIncreasingPath(int negative_edges_count,
                                   int max_train_edges,
                                   int max_train_error,
                                   const vector<vector<long double> >& last_domain) {
  return min(min(negative_edges_count, min(max_train_error, max_train_edges / 2)),
             static_cast<int>(last_domain.size()));
}

int GetMaxTrainError(int negative_edges_count,
                     int max_train_edges,
                     int max_train_error,
                     const vector<vector<long double> >& last_domain) {
  //  return min(min(negative_edges_count, min(max_train_error, max_train_edges / 2)),
  //             static_cast<int>(last_domain.size()));
  return min(min(negative_edges_count, max_train_edges / 2),
             static_cast<int>(last_domain.size()));
}

int GetMaxDeltaIncreasingPath(const vector<vector<long double> >& last_domain,
                              int train_error,
                              int first_alg_index,
                              int current_alg_index,
                              int max_train_edges,
                              int negative_edges_count) {
//  return min(static_cast<int>(last_domain[train_error].size()),
//             min(abs(current_alg_index - first_alg_index) - negative_edges_count,
//                 max_train_edges - 2 * train_error));
  return min(static_cast<int>(last_domain[train_error].size()),
             min(abs(current_alg_index - first_alg_index) - negative_edges_count,
                 max_train_edges));

}


int GetMaxDeltaDecreasingPath(const vector<vector<long double> >& last_domain,
                              int train_error,
                              int first_alg_index,
                              int current_alg_index,
                              int max_train_edges,
                              int negative_edges_count) {
//  int max_delta = min(abs(current_alg_index - first_alg_index) - negative_edges_count,
//                      max_train_edges - 2 * train_error);
  int max_delta = min(abs(current_alg_index - first_alg_index) - negative_edges_count,
                      max_train_edges);
  if (train_error == 0) {
    max_delta = min(max_delta,
                    static_cast<int>(last_domain[train_error].size()));
  }
  else {
    if (train_error < last_domain.size()) {
      max_delta = min(max_delta,
                      static_cast<int>(max(last_domain[train_error - 1].size() - 1,
                                           last_domain[train_error].size())));
    }
    else {
      max_delta = min(max_delta, static_cast<int>(last_domain[train_error - 1].size() - 1));
    }
  }
  return max_delta;
}


int GetMaxTrainErrorDecreasingPath(int negative_edges_count,
                                   int max_train_edges,
                                   int max_train_error,
                                   const vector<vector<long double> >& last_domain) {
  return min(min(negative_edges_count, min(max_train_error, max_train_edges / 2)),
             static_cast<int>(last_domain.size()));
}


