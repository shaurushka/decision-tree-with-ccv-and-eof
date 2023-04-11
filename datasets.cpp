#include "datasets.h"
#include "base_functions.h"
#include "chain_with_gaps.h"
#include "generate_chains.h"


// template<typename T>
// FeatureSample< T >::FeatureSample()
//   : feature(0), class_label(0) {}
  

// template<typename T>
// FeatureSample< T >:: FeatureSample(T feature_, int class_label_) {
//   feature = feature_;
//   class_label = class_label_;
// }


// template<typename T>
// FeatureSample< T >::FeatureSample(const FeatureSample<T>&f){
//     this->feature = f.feature;
//     this->class_label = f.class_label;
// }
    
  
// template<typename T>
// bool operator < (const FeatureSample<T>& first, const FeatureSample<T>& second) {
//   return ((first.feature < second.feature) ||
//           (fabs(first.feature - second.feature) < 1e-9 && first.class_label < second.class_label));
// }

FeatureSample::FeatureSample()
  : feature(0), class_label(0) {}
  

FeatureSample:: FeatureSample(float feature_, int class_label_) {
  feature = feature_;
  class_label = class_label_;
}
    
// FeatureSample:: FeatureSample(const FeatureSample& f){
//     this->feature = f.feature;
//     this->class_label = f.class_label;
// }  

bool operator < (const FeatureSample& first, const FeatureSample& second) {
  return ((first.feature < second.feature) ||
          (fabs(first.feature - second.feature) < 1e-9 && first.class_label < second.class_label));
}

int getZeroClassSize(const vector<FeatureSample>& dataset) {
  int zeros_count = 0;
  for (int index = 0; index < dataset.size(); ++index) {
    zeros_count += (1 - dataset[index].class_label);
  }
  return zeros_count;
}

vector<FeatureSample> sortDataset(const vector<FeatureSample>& dataset) {
  vector<FeatureSample> sorted_dataset(dataset);
  sort(sorted_dataset.begin(), sorted_dataset.end());
  return sorted_dataset;
}


ThresholdClassifier::ThresholdClassifier(float threshold, int errors): threshold(threshold), errors(errors) {}


bool operator < (const ThresholdClassifier& t1, const ThresholdClassifier& t2) {
  return t1.errors < t2.errors;
}


/* все-таки проверить, как строятся gaps */
// template<typename T>
vector<ThresholdClassifier> getThresholdClassifierErrors(
                                    const vector<FeatureSample>& singleFeatureDataset,
                                    ChainMatrixType* error_matrix) {
  vector<FeatureSample> sorted_dataset(sortDataset(singleFeatureDataset));
  // for (int i = 0; i < sorted_dataset.size(); ++i) {
  //   std::cout << "[" << sorted_dataset[i].feature << ";" << sorted_dataset[i].class_label << "]; ";
  // }
  // std::cout << std::endl;
  
  size_t dataset_size = sorted_dataset.size();
  vector<ThresholdClassifier> classifiers;
  classifiers.reserve(dataset_size + 1);
  
  vector<int> errors;
  errors.reserve(dataset_size + 1);

  vector<ChainDescription> gaps;
  
  float prev_feature_value = sorted_dataset[0].feature;
  int classifier_error = getZeroClassSize(sorted_dataset);

  classifiers.push_back(ThresholdClassifier(prev_feature_value - 1, classifier_error));
  errors.push_back(classifier_error);
  
  int gap_begin = 1, gap_end = 0;
  float feature_value;
  int index = 0;
  while (index < dataset_size) {
    feature_value = sorted_dataset[index].feature;
    if (fabs(feature_value - prev_feature_value) > EPS) {
      classifiers.push_back(ThresholdClassifier((prev_feature_value + feature_value) / 2,
                                                 classifier_error));
      if (gap_end > gap_begin) {
        gaps.push_back(ChainDescription(gap_begin, gap_end));
      }
      gap_begin = gap_end = index + 1;
      prev_feature_value = feature_value;
     }
    else {
      ++gap_end;
    }
    classifier_error += 2 * sorted_dataset[index].class_label - 1;
    errors.push_back(classifier_error);
    ++index;
  }
  

  if (gap_end > gap_begin) {
    gaps.push_back(ChainDescription(gap_begin, gap_end));
  }
  
  *error_matrix = makeGapsFromDescription(GenerateChainMatrixFromErrorVector(dataset_size, errors), gaps);
  classifiers.push_back(ThresholdClassifier(feature_value + 1, classifier_error));
  return classifiers;
}


// template<typename T>
// int getZeroClassSize(const vector<FeatureSample<T> >& dataset) {
//   int zeros_count = 0;
//   for (int index = 0; index < dataset.size(); ++index) {
//     zeros_count += (1 - dataset[index].class_label);
//   }
//   return zeros_count;
// }


// template<typename T>
// vector<FeatureSample<T> > sortDataset(const vector<FeatureSample<T> >& dataset) {
//   vector<FeatureSample<T> > sorted_dataset(dataset);
//   sort(sorted_dataset.begin(), sorted_dataset.end());
//   return sorted_dataset;
// }


// template<typename T>
// ThresholdClassifier<T>::ThresholdClassifier(T threshold, int errors): threshold(threshold), errors(errors) {}


// template<typename T>
// bool operator < (const ThresholdClassifier<T>& t1, const ThresholdClassifier<T>& t2) {
//   return t1.errors < t2.errors;
// }


// /* все-таки проверить, как строятся gaps */
// template<typename T>
// vector<ThresholdClassifier<T> > getThresholdClassifierErrors(
//                                     const vector<FeatureSample<T> >& singleFeatureDataset,
//                                     ChainMatrixType* error_matrix) {
//   vector<FeatureSample<T> > sorted_dataset(sortDataset(singleFeatureDataset));
//   size_t dataset_size = sorted_dataset.size();
//   vector<ThresholdClassifier<T> > classifiers;
//   classifiers.reserve(dataset_size + 1);
  
//   vector<int> errors;
//   errors.reserve(dataset_size + 1);
//   vector<ChainDescription> gaps;
  
//   T prev_feature_value = sorted_dataset[0].feature;
//   int classifier_error = getZeroClassSize(sorted_dataset);
//   classifiers.push_back(ThresholdClassifier<T>(prev_feature_value - 1, classifier_error));
//   errors.push_back(classifier_error);
  
//   int gap_begin = 1, gap_end = 0;
//   T feature_value;
//   int index = 0;
//   while (index < dataset_size) {
//     feature_value = sorted_dataset[index].feature;
//     if (fabs(feature_value - prev_feature_value) > EPS) {
//       classifiers.push_back(ThresholdClassifier<T>((prev_feature_value + feature_value) / 2,
//                                                    classifier_error));
//       if (gap_end > gap_begin) {
//         gaps.push_back(ChainDescription(gap_begin, gap_end));
//       }
//       gap_begin = gap_end = index + 1;
//       prev_feature_value = feature_value;
//      }
//     else {
//       ++gap_end;
//     }
//     classifier_error += 2 * sorted_dataset[index].class_label - 1;
//     errors.push_back(classifier_error);
//     ++index;
//   }
//   if (gap_end > gap_begin) {
//     gaps.push_back(ChainDescription(gap_begin, gap_end));
//   }
//   *error_matrix = makeGapsFromDescription(GenerateChainMatrixFromErrorVector(dataset_size, errors), gaps);
//   classifiers.push_back(ThresholdClassifier<T>(feature_value + 1, classifier_error));
//   return classifiers;
// }



