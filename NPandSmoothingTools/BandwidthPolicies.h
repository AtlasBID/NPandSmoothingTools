
#ifndef BANDWIDTH_POLICIES_H
#define BANDWIDTH_POLICIES_H

#include "TObject.h"

namespace Analysis {

template <typename LPKE, typename LPKETraits>
class GlobalBandwidthPolicy {
  using Numeric_t             = typename LPKETraits::Numeric_t;
  using NumericalVecT         = typename LPKETraits::NumericalVector_t;
  using NumericalVecIndexT    = typename LPKETraits::NumericalVectorIndex_t;
  using ScaledDataStorage_t   = typename LPKETraits::ScaledDataStorage_t;
  using Covariates_t         = typename LPKETraits::Covariates_t;
  using CovariatesIndex_t    = typename LPKETraits::CovariatesIndex_t;

  mutable const ScaledDataStorage_t   *m_data       = nullptr;
  mutable const Covariates_t          *m_coord      = nullptr;
  mutable NumericalVecT                m_bandwidths = NumericalVecT(0);
protected:

  void Update () const
  {
    auto dim = static_cast<const LPKE&>(*this).m_dim;
    if (dim != static_cast<decltype(dim)>(m_bandwidths.GetNoElements()))
      m_bandwidths.ResizeTo(dim);
    m_data = &static_cast<const LPKE&>(*this).m_scaled_coordinates;
    m_coord = &static_cast<const LPKE&>(*this).m_scaled_coordinate;
  }

  Numeric_t GetKernelFactor () const
  {
    Numeric_t bandwidth_product = 1;
    for (NumericalVecIndexT i = 0; i < m_bandwidths.GetNoElements(); ++i)
      bandwidth_product *= m_bandwidths[i];
    return 1.0/bandwidth_product;
  }

  Covariates_t GetKernelInput (const Covariates_t &scaled_diff) const
  {
    Covariates_t scaled_coords(m_bandwidths.GetNoElements());
    for (NumericalVecIndexT dim = 0; dim < m_bandwidths.GetNoElements(); ++dim)
      scaled_coords[dim] = scaled_diff[dim]/m_bandwidths[dim];
    return scaled_coords;
  }

public:
  void SetBandwidth (const NumericalVecIndexT &index, const Numeric_t &val)
  {this->m_bandwidths[index] = val;}

  constexpr Numeric_t GetBandwidth (const NumericalVecIndexT &index) const noexcept
  {return this->m_bandwidths[index];}

  ClassDefNV(GlobalBandwidthPolicy, 1);
};

// template <typename LPKE, typename LPKETraits>
// class BallonBandwidthPolicy {
// };
// 
// template <typename LPKE, typename LPKETraits>
// class PointwiseBandwidthPolicy {
// };

}

templateClassImp(Analysis::GlobalBandwidthPolicy)

#endif // include gaurds
