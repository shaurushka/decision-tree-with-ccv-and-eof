#ifndef TYPES_DEF_H
#define TYPES_DEF_H

#include <vector>
#include <iostream>

using std::vector;


typedef vector<vector<int> > ChainMatrixType;

typedef vector<vector<long double> > LargeValuesMatrixType;

struct ChainDescription {
  int begin_index;
  int end_index;
  ChainDescription();
  
  ChainDescription(int begin_index_, int end_index_);
  
};

class AlgContribution {
  long double contribution_qeps;
  long double contribution_ccv;
  long double contribution_eof;

public:
  AlgContribution();

  AlgContribution(long double contribution_qeps, long double contribution_ccv, long double contribution_eof);

  void add(long double qeps, long double ccv, long double eof);

  long double getContributionQeps() const;

  long double getContributionCCV() const;

  long double getContributionEOF() const;
};

class GenBound {
  long double value;
  vector<long double> contributions;

public:
  GenBound();
  
  GenBound(size_t chain_size);

  void addContribution(long double contribution, int alg_index);

  long double getValue() const;

  vector<long double> getContributions() const;
};

class GenBoundsPack {
  GenBound qeps;
  GenBound ccv;
  GenBound eof;

public:
  GenBoundsPack();

  GenBoundsPack(size_t chain_size);
  
  void addContribution(const AlgContribution& alg, int alg_index); // long double denominator, 

  GenBound getQeps() const;

  GenBound getCCV() const;

  GenBound getEOF() const;
};

enum MU_TYPE {ERM, MAXD};

enum EXPER_TYPE {EXPER_PROB, EXPER_CCV, EXPER_EOF};

#endif