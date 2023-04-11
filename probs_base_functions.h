#include <vector>
#include <cstdlib>

using std::vector;
using std::min;
using std::max;

enum CHAIN_DIRECTION {LEFT = -1, RIGHT = 1};

vector<int> GetGraphEdges(const vector<vector<int> >& chain);

int LeftSampleErrorCount(const vector<vector<int> >& chain,
                         int d,
                         const vector<int>& graph_edges);

bool PessimististicConditionIsHeld(int first_alg_index,
                                   int current_alg_index,
                                   const CHAIN_DIRECTION& direction);

int GetMaxTrainErrorIncreasingPath(int negative_edges_count,
                                   int max_train_edges,
                                   int max_train_error,
                                   const vector<vector<long double> >& last_domain);

int GetMaxDeltaIncreasingPath(const vector<vector<long double> >& last_domain,
                              int train_error,
                              int first_alg_index,
                              int current_alg_index,
                              int max_train_edges,
                              int negative_edges_count);


int GetMaxDeltaDecreasingPath(const vector<vector<long double> >& last_domain,
                              int train_error,
                              int first_alg_index,
                              int current_alg_index,
                              int max_train_edges,
                              int negative_edges_count);

int GetMaxTrainError(int negative_edges_count,
                     int max_train_edges,
                     int max_train_error,
                     const vector<vector<long double> >& last_domain);


int GetMaxTrainErrorDecreasingPath(int negative_edges_count,
                                   int max_train_edges,
                                   int max_train_error,
                                   const vector<vector<long double> >& last_domain);

