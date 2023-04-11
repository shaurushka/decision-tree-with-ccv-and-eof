#include "chain_with_gaps.h"

#include "base_functions.h"
#include "generate_chains.h"
#include "probs_base_functions.h"
#include "math.h"


class ChainWithCorrectedGaps {
public:
  ChainWithCorrectedGaps() {}
  
  ChainWithCorrectedGaps(const ChainMatrixType& chain_with_gaps) {
    errors_vectors_with_gaps_ = chain_with_gaps;
    vector<int> chain_errors_count = GetChainErrors(chain_with_gaps);
    
    for (int index = 0; index + 1 < chain_errors_count.size(); ++index) {
      chain_with_no_gaps_errors_count_.push_back(chain_errors_count[index]);
      
      if (HammingDistance(chain_with_gaps[index + 1], chain_with_gaps[index]) > 1) {
        int left_classifier_correct_samples_count = 0;
        int left_classifier_incorrect_samples_count = 0;
        
        FindCorrectAndIncorrectSamplesCount(chain_with_gaps[index],
                                            chain_with_gaps[index + 1],
                                            &left_classifier_correct_samples_count,
                                            &left_classifier_incorrect_samples_count);
        
        if (left_classifier_incorrect_samples_count > 0) {
          --left_classifier_incorrect_samples_count;
        }
        else {
          --left_classifier_correct_samples_count;
        }
        
        int added_chain_begin_index = chain_with_no_gaps_errors_count_.size();
        
        AddMonotoneChain(left_classifier_correct_samples_count, INCREASE);
        AddMonotoneChain(left_classifier_incorrect_samples_count, DECREASE);
        
        
        int added_chain_end_index = chain_with_no_gaps_errors_count_.size();
        begin_end_of_added_chains_.push_back(ChainDescription(added_chain_begin_index,
                                                              added_chain_end_index));
      }
    }
    chain_with_no_gaps_errors_count_.push_back(chain_errors_count.back());
    GetIsInInitialChainIndicator();
    MapNewIndicesToOld();
  }
  
  vector<int> GetChainWithNoGapsErrorsCount() const {
    return chain_with_no_gaps_errors_count_;
  }
  
  vector<ChainDescription> GetAddedChainsBorders() const {
    return begin_end_of_added_chains_;
  }
  
  vector<ChainDescription> GetBeginEndOfAddedChains() const {
    return begin_end_of_added_chains_;
  }
  
  int GetOldIndex(int index_in_chain_with_gaps) const {
    return new_index_to_old_map_[index_in_chain_with_gaps];
  }
  
  vector<bool> isInInitialChainIndicator() const {
    return is_in_initial_chain_;
  }
  
private:
  ChainMatrixType errors_vectors_with_gaps_;
  vector<int> chain_with_no_gaps_errors_count_;
  vector<ChainDescription> begin_end_of_added_chains_;
  vector<ChainDescription> begin_end_of_initial_chains_;
  vector<bool> is_in_initial_chain_;
  vector<int> new_index_to_old_map_;
  
  enum MonotoneChainType {DECREASE, INCREASE};
  
  void AddMonotoneChain(size_t chain_size, const MonotoneChainType& type) {
    int errors_count = chain_with_no_gaps_errors_count_.back();
    for (int index = 0; index < chain_size; ++index) {
      if (type == INCREASE) {
        ++errors_count;
      }
      else {
        --errors_count;
      }
      chain_with_no_gaps_errors_count_.push_back(errors_count);
    }
  }
  
  void FindCorrectAndIncorrectSamplesCount(const vector<int>& left_classifier,
                                           const vector<int>& right_classifier,
                                           int* left_classifier_correct_samples_count,
                                           int* left_classifier_incorrect_samples_count) const {
    for (size_t sample_index = 0; sample_index < left_classifier.size(); ++sample_index) {
      if (left_classifier[sample_index] < right_classifier[sample_index]) {
        ++(*left_classifier_correct_samples_count);
      }
      else if (left_classifier[sample_index] > right_classifier[sample_index]) {
        ++(*left_classifier_incorrect_samples_count);
      }
    }
  }
  
  void AddInitialChainEdgeIndicators(int edge_begin_index, int edge_end_index) {
    for (int initial_chain_index = edge_begin_index;
         initial_chain_index < edge_end_index;
         ++initial_chain_index) {
      is_in_initial_chain_[initial_chain_index] = true;
    }
  }
  
  void GetIsInInitialChainIndicator() {
    begin_end_of_initial_chains_.clear();
    int initial_chain_begin = 0;
    int initial_chain_end;
    is_in_initial_chain_.resize(chain_with_no_gaps_errors_count_.size());
    
    
    for (int index = 0; index < begin_end_of_added_chains_.size(); ++index) {
      initial_chain_end = begin_end_of_added_chains_[index].begin_index;
      AddInitialChainEdgeIndicators(initial_chain_begin, initial_chain_end);
      begin_end_of_initial_chains_.push_back(ChainDescription(initial_chain_begin,
                                                              initial_chain_end));
      
      initial_chain_begin = begin_end_of_added_chains_[index].end_index;
    }
    initial_chain_end = chain_with_no_gaps_errors_count_.size();
    AddInitialChainEdgeIndicators(initial_chain_begin, initial_chain_end);
    begin_end_of_initial_chains_.push_back(ChainDescription(initial_chain_begin,
                                                            initial_chain_end));
    
  }
  
  void MapNewIndicesToOld() {
    size_t chain_size = is_in_initial_chain_.size();
    new_index_to_old_map_.assign(chain_size, -1);
    int index_in_chain_with_gaps = 0;
    for (int index = 0; index < chain_size; ++index) {
      if (is_in_initial_chain_[index]) {
        new_index_to_old_map_[index] = index_in_chain_with_gaps;
        ++index_in_chain_with_gaps;
      }
    }
  }
};


long double GetPathsCountFromLastDomain(int delta_from_middle,
                                        int train_error,
                                        const LargeValuesMatrixType& last_domain) {
  if (train_error < last_domain.size()) {
    int last_domain_size = last_domain[train_error].size();
    int last_domain_middle = (last_domain_size - 1) / 2;
    int last_domain_point = last_domain_middle + delta_from_middle;
    if (last_domain_point < last_domain_size && last_domain_point >= 0) {
      return last_domain[train_error][last_domain_point];
    }
  }
  return 0;
}


int GetMinDelta(int max_delta, bool is_in_initial_chain, const MU_TYPE& mu_type) {
  if ((!is_in_initial_chain) || (mu_type == MAXD)){
    return -max_delta;
  }
  return 0;
}

bool deltaIsNotAcceptable(int delta, bool is_in_initial_chain, 
                          int first_alg_err, int current_alg_err, 
                          int l, int L, 
                          const CHAIN_DIRECTION& chain_dir, 
                          const MU_TYPE& mu_type) {
  if (mu_type == ERM) {
    return (delta == 0 &&
            is_in_initial_chain &&
            !PessimististicConditionIsHeld(first_alg_err, current_alg_err, chain_dir));
  }
  // (mu_type == MAXD) 
  else {
    long double lower_bound = (l * 1.0 / L) * (current_alg_err - first_alg_err);
    return (is_in_initial_chain && 
            ((delta < lower_bound) || 
             ((chain_dir == RIGHT) && (fabs(delta - lower_bound) < 1e-10))));
  }
}

enum EDGE_DIR {UP, DOWN};
LargeValuesMatrixType getPathsIntoMiddleDomainCount(int first_alg_err,
                                                    int current_alg_err,
                                                    const LargeValuesMatrixType& last_domain,
                                                    bool is_in_initial_chain,
                                                    int l, 
                                                    int L, 
                                                    const EDGE_DIR& edge_dir,
                                                    const CHAIN_DIRECTION& chain_dir, 
                                                    const MU_TYPE& mu_type, 
                                                    bool* flag) {
  int max_train_err = last_domain.size() - 1;
  LargeValuesMatrixType current_domain(max_train_err + 1);
  for (int train_err = 0; train_err <= max_train_err; ++train_err) {
    current_domain[train_err].assign(2 * max_train_err + 1, 0);
    int min_delta = GetMinDelta(max_train_err, is_in_initial_chain, mu_type);
    for (int delta = min_delta; delta <= max_train_err; ++delta) {
      if (deltaIsNotAcceptable(delta, is_in_initial_chain, first_alg_err, current_alg_err, 
                               l, L, chain_dir, mu_type)) {
        current_domain[train_err][max_train_err + delta] = 0;
      }
      else {
        long double new_term;
        if (edge_dir == UP) {
          new_term = GetPathsCountFromLastDomain(delta - 1, train_err, last_domain) +
            GetPathsCountFromLastDomain(delta, train_err, last_domain);
        }
        else {
          new_term = GetPathsCountFromLastDomain(delta, train_err, last_domain) +
            GetPathsCountFromLastDomain(delta + 1, train_err - 1, last_domain);
        }
        current_domain[train_err][max_train_err + delta] += new_term;
        if (new_term > 1e20) {
          (*flag) = true;
        }
      }
    }
  }
  return current_domain;
}

size_t getSubChainEdgesCount(const CHAIN_DIRECTION& dir, int first_alg_index, int full_size) {
  if (dir == LEFT) {
    return first_alg_index;
  }
  else {
    return full_size - 1 - first_alg_index;
  }
}

LargeValuesMatrixType getPathsIntoFinalDomainCount(const ChainWithCorrectedGaps& chain_with_corrected_gaps,
                                                   int first_alg_index,
                                                   int l, 
                                                   int L,
                                                   const CHAIN_DIRECTION& step_in_dir, 
                                                   const MU_TYPE& mu_type, 
                                                   long double log_denominator) {
  vector<int> chain_errors = chain_with_corrected_gaps.GetChainWithNoGapsErrorsCount();
  /* subchain_size - это количество ребер в цепи! */
  size_t subchain_size = getSubChainEdgesCount(step_in_dir, first_alg_index, chain_errors.size());
  
  LargeValuesMatrixType last_domain(subchain_size + 1, vector<long double>(2 * subchain_size + 1));
  last_domain[0][subchain_size] = 1;

  bool flag_has_large_values = false;
  bool flag_has_been_modified = false;
  for (int shift = 1; shift <= subchain_size; ++shift) {
    int current_alg_index = first_alg_index + step_in_dir * shift;
    bool is_in_initial_chain = chain_with_corrected_gaps.isInInitialChainIndicator()[current_alg_index];
    
    EDGE_DIR edge_dir = DOWN;
    if (chain_errors[current_alg_index] > chain_errors[current_alg_index - step_in_dir]) {
      edge_dir = UP;
    }

    last_domain = getPathsIntoMiddleDomainCount(chain_errors[first_alg_index],
                                                chain_errors[current_alg_index],
                                                last_domain,
                                                is_in_initial_chain,
                                                l, 
                                                L, 
                                                edge_dir,
                                                step_in_dir, 
                                                mu_type, 
                                                &flag_has_large_values);
    if ((flag_has_large_values == true) && (flag_has_been_modified == false)) {
      for (int i = 0; i < last_domain.size(); i++) {
        for (int j = 0; j < last_domain[i].size(); j++) {
          last_domain[i][j] /= std::expl(log_denominator);
        }
      }
      flag_has_been_modified = true;
    }
  }
  if (flag_has_been_modified == false) {
    for (int i = 0; i < last_domain.size(); i++) {
      for (int j = 0; j < last_domain[i].size(); j++) {
        last_domain[i][j] /= std::expl(log_denominator);
      }
    }
  }
  return last_domain;
}

LargeValuesMatrixType getSubChainSplitsCount(const ChainWithCorrectedGaps& chain_with_corrected_gaps,
                                             int first_alg_index,
                                             int l, 
                                             int L, 
                                             const CHAIN_DIRECTION& chain_dir, 
                                             const MU_TYPE& mu_type, 
                                             long double log_denominator) {
  LargeValuesMatrixType pathsToLastDomain(getPathsIntoFinalDomainCount(chain_with_corrected_gaps,
                                                                       first_alg_index,
                                                                       l, 
                                                                       L, 
                                                                       chain_dir, 
                                                                       mu_type, 
                                                                       log_denominator));
  int subchain_size = pathsToLastDomain.size() - 1;
  LargeValuesMatrixType splits_count(subchain_size + 1,
                                     vector<long double>(subchain_size + 1, 0));
  for (int train_edges = 0; train_edges <= subchain_size; ++train_edges) {
    for (int train_err = 0; train_err <= train_edges; ++train_err) {
      splits_count[train_edges][train_err] = 
        pathsToLastDomain[train_err][subchain_size + train_edges - 2 * train_err];
    }
  }
  return splits_count;
}

LargeValuesMatrixType getNeutralSplitsCount(int train_edges_max, int alg_chain_err, int alg_neutral_err,
                                            int L, int l, int chain_edges_count, long double eps) {
  LargeValuesMatrixType splits_count(train_edges_max + 1);
  for (int train_edges = 0; train_edges <= train_edges_max; ++train_edges) {
    int train_err_max = fmin(train_edges, alg_chain_err);
    splits_count[train_edges].resize(train_err_max + 1);
    for (int train_err = 0; train_err <= train_err_max; ++train_err) {
      int s = getIntOfValue(l * (alg_neutral_err + alg_chain_err - eps * (L - l)) / L) - train_err;
      if (s >= 0) {
        splits_count[train_edges][train_err] = HH(L, l, chain_edges_count, alg_neutral_err, train_edges, s);
      }
    }
  }
  return splits_count;
}

LargeValuesMatrixType getNeutralValidError(int train_edges_max, int alg_chain_err, int m,
                                           int L, int l, int chain_edges_count) {
  LargeValuesMatrixType neutral_valid_err(train_edges_max + 1, vector<long double>(train_edges_max + 1));

  for (int train_edges = 0; train_edges <= train_edges_max; ++train_edges)
    for (int train_err = 0; train_err <= fmin(alg_chain_err, train_edges); ++train_err)
      for (int s = 0; s <= fmin(l - train_edges, m); ++s) 
        neutral_valid_err[train_edges][train_err] += combination(s, m) * 
                                                     combination(l - train_edges - s, L - chain_edges_count - m) * 
                                                     (alg_chain_err + m - s - train_err) * 1.0 / (L - l);
  return neutral_valid_err;
}

LargeValuesMatrixType getDiscrepancy(int train_edges_max, int alg_chain_err, int m,
                                     int L, int l, int chain_edges_count) {
  LargeValuesMatrixType discrepancy(train_edges_max + 1, vector<long double>(train_edges_max + 1));

  for (int train_edges = 0; train_edges <= train_edges_max; ++train_edges)
    for (int train_err = 0; train_err <= fmin(alg_chain_err, train_edges); ++train_err)
      for (int s = 0; s <= fmin(l - train_edges, m); ++s) 
        discrepancy[train_edges][train_err] += 
          combination(s, m) * combination(l - train_edges - s, L - chain_edges_count - m) * 
          fmax(alg_chain_err + m - (s + train_err) * 1.0 * L / l, 0) * 1.0 / (L - l);
  return discrepancy;
}

enum Functional {QEPS, CCV, EOFF};
LargeValuesMatrixType getFunctionalValues(int train_edges_max, int alg_chain_err, int alg_neutral_err,
                                          int L, int l, int chain_edges_count, long double eps, 
                                          const Functional& functional) {
  LargeValuesMatrixType functional_values;
  switch (functional) {
    case QEPS : 
      functional_values = getNeutralSplitsCount(train_edges_max, alg_chain_err, alg_neutral_err,
                                                L, l, chain_edges_count, eps);
      break;
    case CCV : 
      functional_values = getNeutralValidError(train_edges_max, alg_chain_err, alg_neutral_err,
                                               L, l, chain_edges_count);
      break;
    case EOFF : 
      functional_values = getDiscrepancy(train_edges_max, alg_chain_err, alg_neutral_err,
                                         L, l, chain_edges_count);
      break;
  }
  return functional_values;
}

void getAlgContribution(int alg_index,
                        const ChainWithCorrectedGaps& corrected_chain,
                        const vector<int> &corrected_chain_errors,
                        const vector<vector<int> > &corrected_chain_matrix,
                        vector<int> &graph_edges,
                        int chain_edges_count,
                        int L,
                        int l,
                        int m,
                        long double eps, 
                        const MU_TYPE& mu_type, 
                        long double log_denominator,
                        AlgContribution* contribution) {
  long double half_log_denom = log_denominator / 2;
  LargeValuesMatrixType left_splits = getSubChainSplitsCount(corrected_chain, alg_index, l, L, LEFT, mu_type, half_log_denom);
  LargeValuesMatrixType right_splits = getSubChainSplitsCount(corrected_chain, alg_index, l, L, RIGHT, mu_type, half_log_denom);
  int alg_chain_err = corrected_chain_errors[alg_index] - m;
  int train_edges_max = fmin(l, chain_edges_count);
  LargeValuesMatrixType neutral_splits(getFunctionalValues(train_edges_max, alg_chain_err, m,
                                                           L, l, chain_edges_count, eps, QEPS));
  LargeValuesMatrixType valid_err(getFunctionalValues(train_edges_max, alg_chain_err, m,
                                                      L, l, chain_edges_count, eps, CCV));
  LargeValuesMatrixType discrepancy(getFunctionalValues(train_edges_max, alg_chain_err, m,
                                                        L, l, chain_edges_count, eps, EOFF));
  int left_chain_err = LeftSampleErrorCount(corrected_chain_matrix, alg_index, graph_edges);
  
  for (int left_train_edges = 0; left_train_edges <= fmin(l, alg_index); ++left_train_edges) {
    for (int right_train_edges = 0; right_train_edges <= fmin(l, chain_edges_count - alg_index); ++right_train_edges) {
      for (int left_train_err = 0; left_train_err <= fmin(left_train_edges, left_chain_err); ++left_train_err) {
        for (int right_train_err = 0; right_train_err <= fmin(right_train_edges, alg_chain_err - left_chain_err);
             ++right_train_err) {
          int train_edges = left_train_edges + right_train_edges;
          int train_err = left_train_err + right_train_err;
          if ((train_edges <= train_edges_max) and (train_err <= fmin(train_edges, alg_chain_err))) {
            long double chain_splits = left_splits[left_train_edges][left_train_err] *
                                       right_splits[right_train_edges][right_train_err];
            // std::cout << "(" << left_train_edges << ";" << 
            //                     right_train_edges << ";" << 
            //                     left_train_err << ";" << 
            //                     right_train_err << ":" << chain_splits << ";" <<
            //                     left_splits[left_train_edges][left_train_err] << ";" <<
            //                     right_splits[right_train_edges][right_train_err] << ";" << 
            //                     // neutral_splits[train_edges][train_err] << ";" <<
            //                     log_denominator << "\n";
            contribution->add(chain_splits * neutral_splits[train_edges][train_err], 
                              chain_splits * valid_err[train_edges][train_err], 
                              chain_splits * discrepancy[train_edges][train_err]); 
          }
        }
      }
    }
  }
}


void getChainWithGapsGenBounds(int L, int l, int m, long double eps,
                               const ChainMatrixType& chain_with_gaps,
                               GenBoundsPack* gen_bounds, 
                               const MU_TYPE& mu_type) {
  if (chain_with_gaps.size() <= 1) {
    throw std::bad_exception();
  }
  ChainWithCorrectedGaps chain_with_corrected_gaps(chain_with_gaps);
  vector<int> corrected_chain_errors(chain_with_corrected_gaps.GetChainWithNoGapsErrorsCount());
  vector<vector<int> > corrected_chain_matrix(GenerateChainMatrixFromErrorVector(L, corrected_chain_errors));
  vector<int> graph_edges = GetGraphEdges(corrected_chain_matrix);
  
  long double log_denominator = logCombination(l, L);
  // cout << "denominator:" << log_denominator << "\n";
  (*gen_bounds) = GenBoundsPack(chain_with_gaps.size());
  
  size_t chain_size = corrected_chain_errors.size();
  for (int alg_index = 0; alg_index < chain_size; ++alg_index) {
    if (chain_with_corrected_gaps.isInInitialChainIndicator()[alg_index]) {
      AlgContribution contribution;
      getAlgContribution(alg_index,
                         chain_with_corrected_gaps,
                         corrected_chain_errors,
                         corrected_chain_matrix,
                         graph_edges,
                         chain_size - 1,
                         L, l, m, eps, 
                         mu_type,
                         log_denominator, 
                         &contribution);
      gen_bounds->addContribution(contribution, chain_with_corrected_gaps.GetOldIndex(alg_index)); //denominator
    }
  }
}

long double getChainWithGapsQeps(int L, int l, int m, long double eps,
                                 const ChainMatrixType& chain_with_gaps,
                                 vector<long double>* contributions, 
                                 const MU_TYPE& mu_type) {
  GenBoundsPack gen_bounds;
  getChainWithGapsGenBounds(L, l, m, eps, chain_with_gaps, &gen_bounds, mu_type); 
  if (contributions) {
    (*contributions) = gen_bounds.getQeps().getContributions();
  }
  return gen_bounds.getQeps().getValue();
}

long double getChainWithGapsCCV(int L, int l, int m,
                                const ChainMatrixType& chain_with_gaps,
                                vector<long double>* contributions, 
                                const MU_TYPE& mu_type) {
  long double eps_dummy = 0.05;
  GenBoundsPack gen_bounds;
  getChainWithGapsGenBounds(L, l, m, eps_dummy, chain_with_gaps, &gen_bounds, mu_type); 
  if (contributions) {
    (*contributions) = gen_bounds.getCCV().getContributions();
  }
  return gen_bounds.getCCV().getValue();
}

long double getChainWithGapsEOF(int L, int l, int m,
                                const ChainMatrixType& chain_with_gaps,
                                vector<long double>* contributions, 
                                const MU_TYPE& mu_type) {
  long double eps_dummy = 0.05;
  GenBoundsPack gen_bounds;
  getChainWithGapsGenBounds(L, l, m, eps_dummy, chain_with_gaps, &gen_bounds, mu_type); 
  if (contributions) {
    (*contributions) = gen_bounds.getEOF().getContributions();
  }
  return gen_bounds.getEOF().getValue();
}

uint32_t getMaxGapsCount(size_t chain_size) {
  return (chain_size - 1) / 2;
}

vector<uint32_t> generateSample(size_t full_size, size_t sample_size) {
  vector<uint32_t> indices(full_size - 1);
  std::iota(indices.begin(), indices.end(), 1);
  std::random_shuffle(indices.begin(), indices.end());
  vector<uint32_t> random_sample(indices.begin(), indices.begin() + sample_size);
  std::sort(random_sample.begin(), random_sample.end());
  return random_sample;
}


vector<ChainDescription> generateGapsDescriptions(size_t chain_size, uint32_t gaps_count) {
  vector<ChainDescription> gaps;
  gaps.reserve(gaps_count);
  vector<uint32_t> sample(generateSample(chain_size, 2 * gaps_count));
  
  for (int index = 0; index < sample.size(); index+= 2) {
    gaps.push_back(ChainDescription(sample[index], sample[index + 1]));
  }
  return gaps;
}


ChainMatrixType makeGapsFromDescription(const ChainMatrixType& chain,
                                        const vector<ChainDescription>& gapsDecriptions) {
  ChainMatrixType chain_with_gaps(chain);
  
  for (vector<ChainDescription>::const_reverse_iterator it = gapsDecriptions.rbegin();
       it != gapsDecriptions.rend();
       ++it) {
    int gap_begin = it->begin_index;
    int gap_end = it->end_index;
    chain_with_gaps.erase(chain_with_gaps.begin() + gap_begin,
                          chain_with_gaps.begin() + gap_end);
  }
  return chain_with_gaps;
}

