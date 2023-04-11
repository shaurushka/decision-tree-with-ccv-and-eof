#include <vector>
#include <iostream>
#include <fstream>
#include <list>
#include <iterator>
#include <algorithm> //min, max, reverse, random_shuffle
#include <cmath> //floor, fmin, fabs
#include <numeric> //iota
#include <string>
#include <sstream> //
#include "types_definition.h"

using std::vector;
using std::cout;
using std::list;
using std::min;
using std::max;
using std::string;

ChainMatrixType makeGapsFromDescription(const ChainMatrixType& chain,
                                        const vector<ChainDescription>& gapsDecriptions);

long double getChainWithGapsQeps(int L, int l, int m, long double eps,
                                 const ChainMatrixType& chain_with_gaps,
                                 vector<long double>* contributions = NULL, 
                                 const MU_TYPE& mu_type = ERM);

long double getChainWithGapsCCV(int L, int l, int m,
                                const ChainMatrixType& chain_with_gaps,
                                vector<long double>* contributions = NULL, 
                                const MU_TYPE& mu_type = ERM);

long double getChainWithGapsEOF(int L, int l, int m,
                                const ChainMatrixType& chain_with_gaps,
                                vector<long double>* contributions = NULL, 
                                const MU_TYPE& mu_type = ERM);
