#ifdef __CINT__

#include "NPandSmoothingTools/LocalPolyKernelEstimator.h"
#include "NPandSmoothingTools/HistogramSmoother.h"
#include "NPandSmoothingTools/BandwidthPolicies.h"
#include "NPandSmoothingTools/ScalingPolicies.h"
#include "NPandSmoothingTools/KernelPolicies.h"
#include "NPandSmoothingTools/DynamicHistogramSmoother.h"
#include "NPandSmoothingTools/SmoothingUtils.h"

#pragma link off all globals;
#pragma link off all classes;
#pragma link off all functions;
#pragma link off all typedef;
// #pragma link C++ nestedclass;

#pragma link C++ namespace Analysis;
// #pragma link C++ defined_in "NPandSmoothingTools/LocalPolyKernelEstimator.h";
// #pragma link C++ defined_in "NPandSmoothingTools/HistogramSmoother.h";
// #pragma link C++ defined_in "NPandSmoothingTools/BandwidthPolicies.h";
// #pragma link C++ defined_in "NPandSmoothingTools/ScalingPolicies.h";
// #pragma link C++ defined_in "NPandSmoothingTools/DynamicHistogramSmoother.h";

#pragma link C++ class Analysis::LocalPolyKernelEstimatorTraits<double>+;
// #pragma link C++ class Analysis::DynamicKernelPolicy<Analysis::HistogramSmoother<double, Analysis::DynamicKernelPolicy, Analysis::GlobalBandwidthPolicy, Analysis::DynamicScalingPolicy>, Analysis::LocalPolyKernelEstimatorTraits<double>>+;
// #pragma link C++ class Analysis::GlobalBandwidthPolicy<Analysis::HistogramSmoother<double, Analysis::DynamicKernelPolicy, Analysis::GlobalBandwidthPolicy, Analysis::DynamicScalingPolicy>, Analysis::LocalPolyKernelEstimatorTraits<double>>+;
// #pragma link C++ class Analysis::DynamicScalingPolicy<Analysis::HistogramSmoother<double, Analysis::DynamicKernelPolicy, Analysis::GlobalBandwidthPolicy, Analysis::DynamicScalingPolicy>, Analysis::LocalPolyKernelEstimatorTraits<double>>+;
#pragma link C++ class Analysis::LocalPolyKernelEstimator<double, Analysis::DynamicKernelPolicy, Analysis::GlobalBandwidthPolicy, Analysis::DynamicScalingPolicy>+;
#pragma link C++ class Analysis::HistogramSmoother<double, Analysis::DynamicKernelPolicy, Analysis::GlobalBandwidthPolicy, Analysis::DynamicScalingPolicy>+;
#pragma link C++ typedef Analysis::DynamicHistogramSmoother+;
#pragma link C++ class Analysis::LocalPolyKernelEstimator<double, Analysis::TF1RadialKernelPolicy, Analysis::GlobalBandwidthPolicy, Analysis::TF1ScalingPolicy>+;
#pragma link C++ class Analysis::HistogramSmoother<double, Analysis::TF1RadialKernelPolicy, Analysis::GlobalBandwidthPolicy, Analysis::TF1ScalingPolicy>+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::NumericalVector_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::NumericalVectorIndex_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::NumericalMatrix_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::NumericalMatrixIndex_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::DataStorage_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::DataStorageIndex_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::ScaledDataStorage_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::ScaledDataStorageIndex_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::FullCoord_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::FullCoordIndex_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::Covariates_t+;
#pragma link C++ typedef Analysis::ROOTHistogramSmoother::CovariatesIndex_t+;

// #include <functional>
// typedef double(scaleFunc_t)(const double&);
// #pragma link C++ class std::function<scaleFunc_t>+;

#pragma link C++ function Analysis::eigenDecompCalib(TString, TString, TString, int, bool);
#pragma link C++ function Analysis::smoothCalibrations(TString, unsigned, float, size_t, bool, bool, bool);
#pragma link C++ function Analysis::optimizeLeaveOneOutCrossValidation(const Analysis::ROOTHistogramSmoother&, const typename Analysis::ROOTHistogramSmoother::FullCoordIndex_t&, double);

#endif
