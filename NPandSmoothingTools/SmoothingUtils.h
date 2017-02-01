
#ifndef SMOOTHING_UTILITIES_H
#define SMOOTHING_UTILITIES_H

#include "TString.h"
#include "TMath.h"

#include "NPandSmoothingTools/DynamicHistogramSmoother.h"

namespace Analysis {

void
eigenDecompCalib (TString fileName,
                  TString unsmoothedFileName = "",
                  TString pathName = "",
                  int     numberOfNPs = -1,
                  bool    preserveMergedError = true);


void
smoothCalibrations (TString  fileName,
                    TString  container = ".*_SF",
                    std::vector<int>  bNPset = {5,3,2},
                    std::vector<int>  cNPset = {4,4,3},
                    std::vector<int>  lightNPset = {12,5,4},
                    unsigned pt_bins = 100,
                    float    pt_smoothing = 0.4,
                    size_t   order = 0,
                    bool     smoothB = true,
                    bool     smoothC = true,
                    bool     smoothLigh = true);

void
smoothContinuousCalibrations (TString  fileName,
                              TString  container = ".*_SF",
                              std::vector<int>  bNPset = {5,3,2},
                              std::vector<int>  cNPset = {4,4,3},
                              std::vector<int>  lightNPset = {12,5,4},
                              unsigned pt_bins = 100,
                              float    b_pt_smoothing = 0.4,
                              float    c_pt_smoothing = 0.4,
                              float    light_pt_smoothing = 0.4,
                              size_t   order = 0,
                              bool     smoothB = true,
                              bool     smoothC = true,
                              bool     smoothLigh = true);

// struct LeaveOneOutCrossValidation {};
//
// template<class Smoother, class OptimizerType>
// Smoother::Numeric_t
// findOptimalGlobalBandwidth (const Smoother &s);


template<class Smoother>
typename Smoother::Numeric_t
optimizeLeaveOneOutCrossValidation(const Smoother                            &s,
                                   const typename Smoother::FullCoordIndex_t &selected_coord = 0,
                                   typename Smoother::Numeric_t        delta = 0.0);
}

// //////////////////////////////////////
#include "BandwidthOptimization.icc"

// //////////////////////////////////////

// template double Analysis::optimizeLeaveOneOutCrossValidation(const Analysis::ROOTHistogramSmoother&, const typename Analysis::ROOTHistogramSmoother::FullCoordIndex_t&, double);

#endif // include gaurds

// double                    GetOptimalChi2NdfWidthX (double delta = 0.0, double target = 1.0, int bias = 0);

// Leave-one-out cross validation with a least-squares metric
// template <class Smoother>
// double  GetOptimalLCVWidthX (double delta = 0.0);
