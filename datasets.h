
#include <vector>
#include <cmath> //fabs
#include <algorithm> // sort, min
#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <random> //shuffle
#include "types_definition.h"


using std::vector;
using std::ifstream;
using std::ofstream;
using std::string;
using std::cout;

const long double EPS = 1e-9;


// template<typename T>
// class FeatureSample {
// public:
//   FeatureSample();
  
//   FeatureSample(T feature, int class_label);
  

//   FeatureSample(const FeatureSample<T>&);

//   T feature;
//   int class_label;
// };


// template<typename T>
// struct ThresholdClassifier {
//   ThresholdClassifier<T>(T threshold, int errors);

//   T threshold;
//   int errors;
// };


// template<typename T>
// vector<ThresholdClassifier<T> > getThresholdClassifierErrors(
//                                     const vector<FeatureSample<T> >& singleFeatureDataset,
//                                     ChainMatrixType* error_matrix);

// template<typename T>
// bool operator < (const FeatureSample<T>& first, const FeatureSample<T>& second);


class FeatureSample {
public:
  FeatureSample();
  
  FeatureSample(float feature, int class_label);
  
  // FeatureSample(const FeatureSample&); // = default;

  // FeatureSample(const FeatureSample<T>&);

  float feature;
  int class_label;
};


struct ThresholdClassifier {
  ThresholdClassifier(float threshold, int errors);

  float threshold;
  int errors;
};


vector<ThresholdClassifier> getThresholdClassifierErrors(
                                    const vector<FeatureSample>& singleFeatureDataset,
                                    ChainMatrixType* error_matrix);

// template<typename T>
bool operator < (const FeatureSample& first, const FeatureSample& second);
