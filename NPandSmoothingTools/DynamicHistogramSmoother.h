
#ifndef DYNAMIC_HISTOGRAM_SMOOTHER_H
#define DYNAMIC_HISTOGRAM_SMOOTHER_H

#include "NPandSmoothingTools/HistogramSmoother.h"
#include "NPandSmoothingTools/BandwidthPolicies.h"
#include "NPandSmoothingTools/ScalingPolicies.h"
#include "NPandSmoothingTools/KernelPolicies.h"

namespace Analysis {
// dynamic smoother (must use 'typedef' keyword rather than 'using' for CINT)
// using DynamicHistogramSmoother = ::Analysis::HistogramSmoother<double,
//                                                    ::Analysis::DynamicKernelPolicy,
//                                                    ::Analysis::GlobalBandwidthPolicy,
//                                                    ::Analysis::DynamicScalingPolicy>;

typedef HistogramSmoother<double,
        Analysis::DynamicKernelPolicy,
        Analysis::GlobalBandwidthPolicy,
        Analysis::DynamicScalingPolicy> DynamicHistogramSmoother;

typedef HistogramSmoother<double,
        Analysis::TF1RadialKernelPolicy,
        Analysis::GlobalBandwidthPolicy,
        Analysis::TF1ScalingPolicy> ROOTHistogramSmoother;

// Explicitly instantiate for 1, 2, and 3 dimensions
// (as far as CINT is concerned this class is only valid for these dimensions)
template void DynamicHistogramSmoother::LPKE::SetDataPoint (const DynamicHistogramSmoother::LPKE::FullCoordIndex_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&);
template void DynamicHistogramSmoother::LPKE::SetDataPoint (const DynamicHistogramSmoother::LPKE::FullCoordIndex_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&);
template void DynamicHistogramSmoother::LPKE::SetDataPoint (const DynamicHistogramSmoother::LPKE::FullCoordIndex_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&);
template void DynamicHistogramSmoother::LPKE::SetDataPoint (const DynamicHistogramSmoother::LPKE::FullCoordIndex_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&);
template void DynamicHistogramSmoother::LPKE::AddDataPoint (const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&);
template void DynamicHistogramSmoother::LPKE::AddDataPoint (const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&);
template void DynamicHistogramSmoother::LPKE::AddDataPoint (const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&);
template void DynamicHistogramSmoother::LPKE::AddDataPoint (const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::LPKE::Numeric_t&);
template DynamicHistogramSmoother::LPKE::Numeric_t DynamicHistogramSmoother::LPKE::operator() (const DynamicHistogramSmoother::LPKE::Numeric_t&) const;
template DynamicHistogramSmoother::LPKE::Numeric_t DynamicHistogramSmoother::LPKE::operator() (const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::Numeric_t&) const;
template DynamicHistogramSmoother::LPKE::Numeric_t DynamicHistogramSmoother::LPKE::operator() (const DynamicHistogramSmoother::LPKE::Numeric_t&, const DynamicHistogramSmoother::Numeric_t&, const DynamicHistogramSmoother::Numeric_t&) const;

template void ROOTHistogramSmoother::LPKE::SetDataPoint (const ROOTHistogramSmoother::LPKE::FullCoordIndex_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&);
template void ROOTHistogramSmoother::LPKE::SetDataPoint (const ROOTHistogramSmoother::LPKE::FullCoordIndex_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&);
template void ROOTHistogramSmoother::LPKE::SetDataPoint (const ROOTHistogramSmoother::LPKE::FullCoordIndex_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&);
template void ROOTHistogramSmoother::LPKE::SetDataPoint (const ROOTHistogramSmoother::LPKE::FullCoordIndex_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&);
template void ROOTHistogramSmoother::LPKE::AddDataPoint (const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&);
template void ROOTHistogramSmoother::LPKE::AddDataPoint (const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&);
template void ROOTHistogramSmoother::LPKE::AddDataPoint (const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&);
template void ROOTHistogramSmoother::LPKE::AddDataPoint (const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::Numeric_t&);
template ROOTHistogramSmoother::LPKE::Numeric_t ROOTHistogramSmoother::LPKE::operator() (const ROOTHistogramSmoother::LPKE::Numeric_t&) const;
template ROOTHistogramSmoother::LPKE::Numeric_t ROOTHistogramSmoother::LPKE::operator() (const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&) const;
template ROOTHistogramSmoother::LPKE::Numeric_t ROOTHistogramSmoother::LPKE::operator() (const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&, const ROOTHistogramSmoother::LPKE::Numeric_t&) const;
}

#endif
