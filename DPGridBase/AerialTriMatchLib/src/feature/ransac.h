#ifndef FEA_RANSAC_HEADER
#define FEA_RANSAC_HEADER

#include "feature/defines.h"

FEA_NAMESPACE_BEGIN

/**
 * The function returns the required number of iterations for a desired
 * RANSAC success rate. If w is the probability of choosing one good sample
 * (the inlier ratio), then w^n is the probability that all n samples are
 * inliers. Then k is the number of iterations required to draw only inliers
 * with a certain probability of success, p:
 *
 *          log(1 - p)
 *     k = ------------
 *         log(1 - w^n)
 *
 * Example: For w = 50%, p = 99%, n = 8: k = log(0.001) / log(0.99609) = 1176.
 * Thus, it requires 1176 iterations for RANSAC to succeed with a 99% chance.
 */
int
compute_ransac_iterations (double inlier_ratio,
    int num_samples,
    double desired_success_rate = 0.99);

FEA_NAMESPACE_END

#endif /* FEA_RANSAC_HEADER */