
#ifndef SCALING_POLICIES_H
#define SCALING_POLICIES_H

#include <type_traits>
#include <vector>
#include <functional>

#include "TMath.h"
#include "TF1.h"

namespace Analysis {

template <typename LPKE, typename LPKETraits>
class NoScalingPolicy {
  using Numeric_t = typename LPKETraits::Numeric_t;

protected:
  void Update () const
  {}

public:

  Numeric_t ScaleCoordinate (const size_t &/*dim_index*/, const Numeric_t &coord) const
  {return coord;}

  Numeric_t InverseScaleCoordinate (const size_t &/*dim_index*/, const Numeric_t &coord) const
  {return coord;}
};

template <typename LPKE, typename LPKETraits>
class LnScalingPolicy {
  using Numeric_t = typename LPKETraits::Numeric_t;
  std::vector<bool> m_lnScaleIndices;

protected:
  void Update () const
  {
    auto dim = static_cast<const LPKE&>(*this).m_dim;
    if (dim != static_cast<decltype(dim)>(m_lnScaleIndices.size()))
      m_lnScaleIndices.resize(dim);
    for (auto &f : m_lnScaleIndices)
      f = false;
  }

public:

  Numeric_t ScaleCoordinate (const size_t &dim_index, const Numeric_t &coord) const
  {
    if (dim_index <= m_lnScaleIndices.size()) Update(); 
    return m_lnScaleIndices[dim_index] ? Numeric_t(TMath::Log(coord)) : coord;
  }

  Numeric_t InverseScaleCoordinate (const size_t &dim_index, const Numeric_t &coord) const
  {
    if (dim_index <= m_lnScaleIndices.size()) Update(); 
    return m_lnScaleIndices[dim_index] ? Numeric_t(TMath::Exp(coord)) : coord;
  }

  void SetLnScale (const size_t &dim_index, bool enable = true)
  {m_lnScaleIndices[dim_index] = enable;}
};


template <typename LPKE, typename LPKETraits>
class DynamicScalingPolicy {
  using Numeric_t           = typename LPKETraits::Numeric_t;
  using Func_t              = std::function<Numeric_t(const Numeric_t&)>;
  // using Func_t              = std::function<double(const double&)>;
  using FuncStorage_t       = std::vector<Func_t>;
  using FuncStorageIndex_t  = typename FuncStorage_t::size_type;

  mutable FuncStorage_t  m_scales, m_invscales;

protected:
  void Update () const
  {
    auto dim = static_cast<const LPKE&>(*this).m_dim;
    if (dim != m_scales.size()) {
      m_scales.resize(dim);
      m_invscales.resize(dim);
    }
  }

  Numeric_t ScaleCoordinate (const size_t &dim_index, const Numeric_t &coord) const
  {return m_scales[dim_index](coord);}

  Numeric_t InverseScaleCoordinate (const size_t &dim_index, const Numeric_t &coord) const
  {return m_invscales[dim_index](coord);}

public:


  void SetScaleFunction (const FuncStorageIndex_t &index, const Func_t &func)
  {if (index <= m_scales.size()) Update(); m_scales[index] = func;}

  void SetInvScaleFunction (const FuncStorageIndex_t &index, const Func_t &func)
  {if (index <= m_invscales.size()) Update(); m_invscales[index] = func;}

  const Func_t& GetScaleFunction (const FuncStorageIndex_t &index) const
  {return m_scales[index];}

  const Func_t& GetInvScaleFunction (const FuncStorageIndex_t &index) const
  {return m_invscales[index];}

  ClassDefNV(DynamicScalingPolicy, 1);
};

template <typename LPKE, typename LPKETraits>
class TF1ScalingPolicy {
  mutable std::vector<TF1>  m_scales = {{TF1("no scale", "x[0]", -1, 1)}}, 
                            m_invscales = {{TF1("no scale", "x[0]", -1, 1)}};

protected:
  void Update () const
  {
    auto dim = static_cast<const LPKE&>(*this).m_dim;
    if (dim != m_scales.size()) {
      m_scales.resize(dim);
      m_invscales.resize(dim);
    }
  }

public:

  double ScaleCoordinate (const size_t &dim_index, const double &coord) const
  {return m_scales[dim_index](coord);}

  double InverseScaleCoordinate (const size_t &dim_index, const double &coord) const
  {return m_invscales[dim_index](coord);}

  void SetScaleFunction (const std::vector<TF1>::size_type &index, const TF1 &func)
  {if (index <= m_scales.size()) Update(); m_scales[index] = func;}

  void SetInvScaleFunction (const std::vector<TF1>::size_type &index, const TF1 &func)
  {if (index <= m_invscales.size()) Update(); m_invscales[index] = func;}

  const TF1& GetScaleFunction (const std::vector<TF1>::size_type &index)
  {return m_scales[index];}

  const TF1& GetInvScaleFunction (const std::vector<TF1>::size_type &index)
  {return m_invscales[index];}

  ClassDefNV(TF1ScalingPolicy, 1);
};

}

templateClassImp(Analysis::DynamicScalingPolicy)
templateClassImp(Analysis::TF1DynamicScalingPolicy)

#endif // include gaurds
