#include <Rcpp.h>
#include <cmath>

using namespace Rcpp;

namespace {

const double LOG_SQRT_2PI = 0.91893853320467274178;

bool isMissing(double x) {
  return NumericVector::is_na(x) || std::isnan(x);
}

NumericVector listElementAsVector(const List& values, int index, const char* name) {
  if (index < 0 || index >= values.size()) {
    stop("Index %d is outside '%s'.", index + 1, name);
  }

  SEXP element = values[index];
  if (Rf_isNull(element)) {
    stop("'%s' element %d is NULL.", name, index + 1);
  }

  return as<NumericVector>(element);
}

double logNormalDensity(double x, double mean, double sd, bool includeNormalConstant) {
  if (sd <= 0.0 || isMissing(sd) || isMissing(mean) || isMissing(x)) {
    return R_NegInf;
  }

  const double z = (x - mean) / sd;
  double value = -std::log(sd) - 0.5 * z * z;
  if (includeNormalConstant) {
    value -= LOG_SQRT_2PI;
  }
  return value;
}

double logSumExp(const std::vector<double>& values) {
  if (values.empty()) {
    return R_NegInf;
  }

  double maxValue = R_NegInf;
  for (double value : values) {
    if (value > maxValue) {
      maxValue = value;
    }
  }

  if (!R_finite(maxValue)) {
    return maxValue;
  }

  long double total = 0.0;
  for (double value : values) {
    total += std::exp(value - maxValue);
  }

  return maxValue + std::log(static_cast<double>(total));
}

double componentSd(double sdOrVariance, bool valuesAreVariances) {
  if (valuesAreVariances) {
    if (sdOrVariance < 0.0) {
      return R_NaN;
    }
    return std::sqrt(sdOrVariance);
  }
  return sdOrVariance;
}

double lTargetGMMInternal(double newX,
                         NumericVector means,
                         NumericVector sdOrVariances,
                         Nullable<NumericVector> weights,
                         bool valuesAreVariances,
                         bool includeNormalConstant) {
  if (means.size() != sdOrVariances.size()) {
    stop("'means' and 'sds' must have the same length.");
  }

  NumericVector componentWeights;
  const bool hasWeights = weights.isNotNull();
  if (hasWeights) {
    componentWeights = NumericVector(weights);
    if (componentWeights.size() != means.size()) {
      stop("'weights' must have the same length as 'means' and 'sds'.");
    }
  }

  std::vector<double> terms;
  terms.reserve(means.size());

  for (int k = 0; k < means.size(); ++k) {
    const double sd = componentSd(sdOrVariances[k], valuesAreVariances);
    double logWeight = -std::log(static_cast<double>(means.size()));

    if (hasWeights) {
      const double weight = componentWeights[k];
      if (isMissing(weight) || weight <= 0.0) {
        continue;
      }
      logWeight = std::log(weight);
    }

    const double logDensity = logNormalDensity(newX, means[k], sd, includeNormalConstant);
    if (R_finite(logDensity)) {
      terms.push_back(logWeight + logDensity);
    }
  }

  return logSumExp(terms);
}

Nullable<NumericVector> weightVectorForIndex(Nullable<List> weights, int index) {
  if (weights.isNull()) {
    return R_NilValue;
  }

  List weightList(weights);
  return listElementAsVector(weightList, index, "weights");
}

} // namespace

// [[Rcpp::export]]
double lTargetGMMCpp(int d,
                     double newX,
                     NumericVector means,
                     NumericVector sds,
                     Nullable<NumericVector> weights = R_NilValue,
                     bool sdsAreVariances = false,
                     bool includeNormalConstant = true) {
  // 'd' is kept for API symmetry with lTargetKernelCpp.
  (void)d;
  return lTargetGMMInternal(newX, means, sds, weights, sdsAreVariances, includeNormalConstant);
}

// [[Rcpp::export]]
double loglikelihood_uqGMMCpp(NumericVector params,
                              NumericVector t_times,
                              List means,
                              List sds,
                              Nullable<List> weights = R_NilValue,
                              bool sdsAreVariances = false,
                              bool includeNormalConstant = true) {
  // 'params' is kept for API symmetry with loglikelihood_uqCpp.
  (void)params;

  if (means.size() != sds.size()) {
    stop("'means' and 'sds' must contain the same number of samples.");
  }
  if (t_times.size() != means.size()) {
    stop("'t_times' length must match the number of GMM samples.");
  }
  if (weights.isNotNull()) {
    List weightList(weights);
    if (weightList.size() != means.size()) {
      stop("'weights' must contain one vector per GMM sample.");
    }
  }

  long double sumll = 0.0;
  for (int i = 0; i < t_times.size(); ++i) {
    NumericVector currentMeans = listElementAsVector(means, i, "means");
    NumericVector currentSds = listElementAsVector(sds, i, "sds");
    Nullable<NumericVector> currentWeights = weightVectorForIndex(weights, i);

    sumll += lTargetGMMInternal(t_times[i],
                                currentMeans,
                                currentSds,
                                currentWeights,
                                sdsAreVariances,
                                includeNormalConstant);
  }

  return static_cast<double>(sumll);
}
