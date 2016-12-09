
#ifndef KERNEL_POLICIES_H
#define KERNEL_POLICIES_H

#include <functional>

#include "TMath.h"
#include "TF1.h"

namespace Analysis {

template <typename LPKE, typename LPKETraits>
class NormalKernelPolicy {
  using Numeric_t = typename LPKETraits::Numeric_t;
  using Covariates_t = typename LPKETraits::Covariates_t;
protected:
  void Update () const
  {}

  Numeric_t EvaluateKernel (const Covariates_t &scaled_coord) const
  {
    Numeric_t normalInput = 0;
    for (auto coord : scaled_coord)
      normalInput += coord*coord;
    return Numeric_t(TMath::Gaus(TMath::Sqrt(normalInput), 0.0, 1.0, kTRUE));
  }
};

template <typename LPKE, typename LPKETraits>
class DynamicKernelPolicy {
  using Numeric_t     = typename LPKETraits::Numeric_t;
  using Covariates_t = typename LPKETraits::Covariates_t;
  using Func_t        = std::function<Numeric_t(const Covariates_t&)>;
  Func_t  m_kernel;

protected:

  void Update () const
  {}

  Numeric_t EvaluateKernel (const Covariates_t &scaled_coord) const
  {return m_kernel(scaled_coord);}

public:
  void SetKernel (const Func_t &func)
  {m_kernel = func;}

  const Func_t& GetKernel () const
  {return m_kernel;}

  ClassDefNV(DynamicKernelPolicy, 1);
};

template <typename LPKE, typename LPKETraits>
class TF1RadialKernelPolicy {
  using Covariates_t = typename LPKETraits::Covariates_t;
  TF1  m_kernel = TF1("guassian", "TMath::Gaus(x, 0.0, 1.0, kTRUE)", -1.0, 1.0);

protected:

  void Update () const
  {}

  double EvaluateKernel (const Covariates_t &scaled_coord) const
  {
    double reduced = 0.0;
    for (auto coord : scaled_coord)
      reduced += coord*coord;
    reduced = TMath::Sqrt(reduced);
    return m_kernel(reduced);
  }

public:

  double EvaluateKernel (const double &reduced_scaled_coord) const
  {
    return m_kernel(reduced_scaled_coord);
  }

  void SetKernel (const TF1 &func)
  {m_kernel = func;}

  const TF1& GetKernel () const
  {return m_kernel;}

  ClassDefNV(TF1RadialKernelPolicy, 1);
};

}

templateClassImp(Analysis::DynamicKernelPolicy)

#endif // include gaurds
