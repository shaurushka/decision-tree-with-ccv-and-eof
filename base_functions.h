#include <iostream> //for size_t
#include <cmath> // for fmin
#include <vector>

using std::vector;

long double combination(int of_count, int from_count);

long double logCombination(int of_count, int from_count);

long double H(int L, int l, int m, int s);

long double TildeH(int L, int l, int D, int m, int z_0,  int s);

long double HH(int L, int l, int D, int m, int z_0,  int s);

int ErrorsCount(const vector<vector<int> >& chain, int d);

vector<int> GetChainErrors(const vector<vector<int> >& chain);

void PrintIntVector(const vector<int>& v);

int HammingDistance(const vector<int>& first_vector,
                    const vector<int>& second_vector);

void PrintLongDoubleVector(const vector<long double>& v);

int getIntOfValue(long double value);