#include "types_definition.h"

ChainDescription::ChainDescription() : begin_index(0), end_index(0) {}

ChainDescription::ChainDescription(int begin_index_, int end_index_)
  : begin_index(begin_index_), end_index(end_index_) {}

std::ostream& operator << (std::ostream& os, const ChainDescription& chain) {
  os << "(" << chain.begin_index << ", " << chain.end_index << ")";
  return os;
}

AlgContribution::AlgContribution() : contribution_qeps(0), contribution_ccv(0), contribution_eof(0) {}

AlgContribution::AlgContribution(long double contribution_qeps, 
                                 long double contribution_ccv, 
                                 long double contribution_eof)
  : contribution_qeps(contribution_qeps), 
    contribution_ccv(contribution_ccv),
    contribution_eof(contribution_eof) {}

void AlgContribution::add(long double qeps, long double ccv, long double eof) {
  contribution_qeps += qeps;
  contribution_ccv += ccv;
  contribution_eof += eof;
}

long double AlgContribution::getContributionQeps() const {
  return contribution_qeps;
}

long double AlgContribution::getContributionCCV() const {
  return contribution_ccv;
}

long double AlgContribution::getContributionEOF() const {
  return contribution_eof;
}

GenBound::GenBound() {}
  
GenBound::GenBound(size_t chain_size) : value(0) {
  contributions.assign(chain_size, 0);
}

void GenBound::addContribution(long double contribution, int alg_index) {
  contributions[alg_index] += contribution; 
  value += contribution;
}

long double GenBound::getValue() const {
  return value;
}

vector<long double> GenBound::getContributions() const {
  return contributions;
}

GenBoundsPack::GenBoundsPack() {}

GenBoundsPack::GenBoundsPack(size_t chain_size) 
  : qeps(GenBound(chain_size)), ccv(GenBound(chain_size)), eof(GenBound(chain_size)) {}

void GenBoundsPack::addContribution(const AlgContribution& alg, int alg_index) { //, long double denominator
  qeps.addContribution(alg.getContributionQeps(), alg_index); //  / denominator
  ccv.addContribution(alg.getContributionCCV(), alg_index); //  / denominator
  eof.addContribution(alg.getContributionEOF(), alg_index); //  / denominator
}

GenBound GenBoundsPack::getQeps() const {
  return qeps;
}

GenBound GenBoundsPack::getCCV() const {
  return ccv;
}

GenBound GenBoundsPack::getEOF() const {
  return eof;
}
