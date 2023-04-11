#include "base_functions.h"
#include "math.h"

long double combination(int of_count, int from_count) {
  if ((from_count < of_count) || (of_count < 0) || (from_count < 0)) {
    return 0;
  }
  
  /* use symmetry of combination to optimize calculation */
  if (from_count - of_count < of_count) {
    of_count = from_count - of_count;
  }
  // std:: cout << of_count << " " << from_count << "\n";
  long double current_combination = 1;
  for (size_t of = 1; of <= of_count; ++of) {
    current_combination = current_combination / of;
    current_combination = current_combination * (from_count - of + 1);
    // std:: cout << current_combination << " ";
  }

  return current_combination;
}


long double logCombination(int of_count, int from_count) {
  if ((from_count < of_count) || (of_count < 0) || (from_count < 0)) {
    return -1e307;
  }
  
  /* use symmetry of combination to optimize calculation */
  if (from_count - of_count < of_count) {
    of_count = from_count - of_count;
  }
  // std:: cout << of_count << " " << from_count << "\n";
  long double current_combination = 0;
  for (size_t of = 1; of <= of_count; ++of) {
    current_combination += std::logl(from_count - of + 1) -  std::logl(of);
    // std:: cout << current_combination << " ";
  }
  
  return current_combination;
}


long double H(int L, int l, int m, int s) {
  if ((s < 0) || (L < l) || (l < 0) || (L < 0)) {
    return 0;
  }
  long double sum = 0;
  long double log_denominator = logCombination(l, L);
  for (int i = fmax(0, m - L + l); i <= fmin(fmin(l, m), s); ++i) {
    long double log_term = logCombination(i, m) + logCombination(l - i, L - m) - log_denominator;
    if (log_term > 1e-300) {
      sum += std::expl(log_term);  
    }
    // sum += combination(i, m) * combination(l - i, L - m);
  }
  return sum;
  // return sum / combination(l, L);
}

long double TildeH(int L, int l, int D, int m, int z_0,  int s) {
  if ((s < 0)||(L < l)||(l < 0) || (L < 0)) {
    return 0;
  }
  long double multH = H(L - D, l - z_0, m, s);
  if (multH > 1e-300) {
    multH = std::logl(multH);
    long double log_term = logCombination(l - z_0, L - D) * multH - logCombination(l, L);
    return std::expl(log_term);
  }
  return 0;
  
  
  // return combination(l - z_0, L - D) * H(L - D, l - z_0, m, s) / combination(l, L);

}

long double HH(int L, int l, int D, int m, int z_0,  int s) {
  if ((s < 0)||(L < l)||(l < 0) || (L < 0)) {
    return 0;
  }
  long double multH = H(L - D, l - z_0, m, s);
  if (multH > 1e-300) {
    multH = std::logl(multH);
    long double log_term = logCombination(l - z_0, L - D) * multH;
    return std::expl(log_term);
  }
  return 0;
  // return combination(l - z_0, L - D) * H(L - D, l - z_0, m, s);
}

int ErrorsCount(const vector<vector<int> >& chain, int d) {
  int errors = 0;
  for (int index = 0; index < chain[d].size(); ++index) {
    errors += chain[d][index];
  }
  return errors;
}

vector<int> GetChainErrors(const vector<vector<int> >& chain) {
  vector<int> chain_errors(chain.size());
  for (int index = 0; index < chain.size(); ++index) {
    chain_errors[index] = ErrorsCount(chain, index);
  }
  return chain_errors;
}

int HammingDistance(const vector<int>& first_vector,
                    const vector<int>& second_vector) {
  if (first_vector.size() != second_vector.size()) {
    throw std::bad_exception();
  }
  int distance = 0;
  for (int index = 0; index < first_vector.size(); ++index) {
    distance += std::abs(first_vector[index] - second_vector[index]);
  }
  return distance;
}


void PrintIntVector(const vector<int>& v) {
  for (int index = 0; index < v.size(); ++index) {
    std::cout << v[index] << " ";
  }
  std::cout << "\n";
}


void PrintLongDoubleVector(const vector<long double>& v) {
  for (int i = 0; i < v.size(); ++i) {
    std::cout << v[i] << " ";
  }
  std::cout << "\n";
}


int getIntOfValue(long double value) {
  int ceil_of_value = ceil(value);
  if (ceil_of_value - value < 1e-9) {
    return ceil_of_value;
  }
  else {
    return floor(value);
  }
}
