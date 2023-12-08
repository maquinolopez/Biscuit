#include <Rcpp.h>
#include <cstdio>
#include <math.h>
#include <vector>
#include <cmath> // For pow, log, and sum

using namespace Rcpp;

// Helper function to calculate normal density
  
double NormalDensity(double x, double mean, double std_dev) {
    // Function: NormalDensity
    // Purpose: Computes the density of a normal distribution for a given value.
    //
    // Parameters:
    // - double x: The value for which the density is to be calculated.
    // - double mean: The mean of the normal distribution.
    // - double std_dev: The standard deviation of the normal distribution.
    //
    // Returns:
    // - A double value representing the density of the normal distribution at the point 'x'.
    // _____________________________________________
    // Normalize 'x' by subtracting the mean and dividing by the standard deviation
    double nx = (x - mean)/std_dev;
    // Calculate the logarithm of the standard deviation
    double log_sig = log(std_dev );
    // Square of the normalized value 'nx'
    double sqr_x = nx *nx;
    // Compute the logarithm of the density
    double densi = -log_sig - .5 * sqr_x ;
    // Calculate the density using the exponential of 'densi'
    long double result = exp(densi);
    // Return the calculated density
    return result;
    
}

NumericVector getColumn(NumericMatrix tar,int col) {
  int ncol = tar.ncol();
  NumericVector result(ncol);
  
  for(int i = 0; i < ncol; i++){
    result[i] = tar(col, i); // 1 = second column 
  }
  
  return result;
}

// [[Rcpp::export]]
double lTargetKernelCpp(int d, double newX, double kdeValue, NumericVector tar, double sdConvertor) {
    // Function: lTargetKernelCpp
    // Purpose: Computes a log-transformed value based on kernel density estimation.
    //
    // Parameters:
    // - int d: Unused in the function. Placeholder or for future use.
    // - double newX: The new data point for which the kernel density is estimated.
    // - double kdeValue: The bandwidth of the kernel density estimator, controlling the smoothness of the density estimate.
    // - NumericVector tar: A vector of posterior sample used in the kernel density estimation process.
    // - double sdConvertor: A scaling factor to adjust standard deviation for normalization.
    //
    // Returns:
    // - A double value representing the log-transformed result of the kernel density estimation.
    // _____________________________________________
    // Calculate the scaled bandwidth
    double s_2 = sdConvertor * kdeValue;
    // Initialize sum for densities
    long double sumDnorm = 0;
    // Iterate over each element in the 'tar' vector
    for(int i = 0; i < tar.length(); i++) {
      // Current element from 'tar' vector
      double mu = tar[i];
      // Accumulate the normal density for 'newX' given 'mu' and scaled bandwidth 's_2'
      sumDnorm += NormalDensity(newX, mu, s_2);
    }
    // Divide by the sample size
    //sumDnorm = sumDnorm / tar.length(); // remove because is a constant and it don't affect the overall likelihood
    // Compute the log-transformed result
    // Log of the sum of densities minus the log of the product of 'sdConvertor' and 'kdeValue'
    double result = log(sumDnorm) ;//-  log( sdConvertor * kdeValue);
    return result;
}



// [[Rcpp::export]]
double loglikelihood_uqCpp(NumericVector params, NumericVector t_times, NumericVector kde_bw,NumericMatrix tar, double sdConvertor ){
  // Function: loglikelihood_uqCpp
  // Purpose: Calculates the log-likelihood value using kernel density estimation.
  //
  // Parameters:
  // - NumericVector params: Vector of parameters, not used in this function for extended functionality.
  // - NumericVector t_times: Vector of time points or values for which the kernel density is estimated.
  // - NumericVector kde_bw: Vector of bandwidth values for the kernel density estimator, one for each time point in 't_times'.
  // - NumericMatrix tar: A matrix where each column represents a set of data points for kernel density estimation at a corresponding time point.
  // - double sdConvertor: A scaling factor, possibly to adjust standard deviation for normalization.
  //
  // Returns:
  // - A double value representing the sum of log-likelihood values computed for each time point in 't_times'.
  // _____________________________________________
  // Initialize the sum for log-likelihood values
  long double sumll = 0;

  // Iterate over each time point in 't_times'
  for(int i = 0; i < t_times.length(); i++) {
    double newx = t_times[i];
    // Bandwidth value for the current time point
    double bw_value = kde_bw[i];
    // Extract the column vector corresponding to the current time point from 'tar'
    NumericVector tar_vec = getColumn(tar, i);
    // Compute and accumulate the log-likelihood for the current time point
    // using the lTargetKernelCpp function
    sumll += lTargetKernelCpp(i, newx, bw_value, tar_vec, sdConvertor) ;

  }
  
  // Return the total log-likelihood
  return sumll;
} 




// [[Rcpp::export]]
double tdistroC(const std::vector<double>& X, double Mu, double sigma, double a, double b) {
  sigma = pow(sigma, 2);
  double sum = 0.0;
  
  for (double x : X) {
    if (std::isnan(x)) continue; // Skip NaN values
    double term = ((2 * a + 1) / 2.0) * log(b + pow((x - Mu), 2) / (2 * sigma)) + 0.5 * log(sigma);
    sum -= term; // Subtracting because the original R function multiplies by -1
  }
  
  return sum;
}




