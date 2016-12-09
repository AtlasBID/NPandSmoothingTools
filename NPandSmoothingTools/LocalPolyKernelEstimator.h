#ifndef LOCAL_POLY_KERNEL_ESTIMATOR_H
#define LOCAL_POLY_KERNEL_ESTIMATOR_H

#include <cmath>
#include <iostream>
#include <exception>
#include <limits>
#include <string>
#include <vector>
#include <algorithm>

#include "TMath.h"
#include "TError.h"
#include "TNamed.h"
#include "TMatrixT.h"
#include "TVectorT.h"

#include "NPandSmoothingTools/TraitsHelpers.h"


namespace Analysis {
template<typename NumericType>
struct LocalPolyKernelEstimatorTraits {
  using Numeric_t              = NumericType;
  using NumericalVector_t      = TVectorT<Numeric_t>;
  using NumericalMatrix_t      = TMatrixT<Numeric_t>;
  using NumericalVectorIndex_t = decltype(NumericalVector_t().GetNoElements());
  using NumericalMatrixIndex_t = decltype(NumericalMatrix_t().GetNoElements());

  using Coord_t             = std::vector<Numeric_t>;
  using CoordIndex_t        = typename Coord_t::size_type;
  using CoordStorage_t      = std::vector<Coord_t>;
  using CoordStorageIndex_t = typename CoordStorage_t::size_type;

  using DataStorage_t            = CoordStorage_t;
  using DataStorageIndex_t       = typename CoordStorage_t::size_type;
  using ScaledDataStorage_t      = CoordStorage_t;
  using ScaledDataStorageIndex_t = typename CoordStorage_t::size_type;
  using FullCoord_t              = Coord_t;
  using FullCoordIndex_t         = typename Coord_t::size_type;
  using Covariates_t            = Coord_t;
  using CovariatesIndex_t       = typename Coord_t::size_type;

  //ClassDefNV(LocalPolyKernelEstimatorTraits, 1);
};

template<typename NumericType,
         template<class, class> class KernelPolicy,
         template<class, class> class BandwidthPolicy,
         template<class, class> class ScalingPolicy>
class LocalPolyKernelEstimator :
  public TNamed,
  // Policy Bridges
  public KernelPolicy < LocalPolyKernelEstimator<NumericType, KernelPolicy, BandwidthPolicy, ScalingPolicy>, LocalPolyKernelEstimatorTraits < NumericType >>,
  public BandwidthPolicy < LocalPolyKernelEstimator<NumericType, KernelPolicy, BandwidthPolicy, ScalingPolicy>, LocalPolyKernelEstimatorTraits < NumericType >>,
  public ScalingPolicy < LocalPolyKernelEstimator<NumericType, KernelPolicy, BandwidthPolicy, ScalingPolicy>, LocalPolyKernelEstimatorTraits < NumericType >> {
  using Traits_t = LocalPolyKernelEstimatorTraits<NumericType>;

  // Policies
  // friend class KernelPolicy<LocalPolyKernelEstimator,  Traits_t>;
  // friend class BandwidthPolicy<LocalPolyKernelEstimator,  Traits_t>;
  // friend class ScalingPolicy<LocalPolyKernelEstimator,  Traits_t>;
  friend class KernelPolicy < LocalPolyKernelEstimator<NumericType, KernelPolicy, BandwidthPolicy, ScalingPolicy>, LocalPolyKernelEstimatorTraits < NumericType >>;
  friend class BandwidthPolicy < LocalPolyKernelEstimator<NumericType, KernelPolicy, BandwidthPolicy, ScalingPolicy>, LocalPolyKernelEstimatorTraits < NumericType >>;
  friend class ScalingPolicy < LocalPolyKernelEstimator<NumericType, KernelPolicy, BandwidthPolicy, ScalingPolicy>, LocalPolyKernelEstimatorTraits < NumericType >>;

public:

  using Numeric_t                = typename Traits_t::Numeric_t;
  using NumericalVector_t        = typename Traits_t::NumericalVector_t;
  using NumericalVectorIndex_t   = typename Traits_t::NumericalVectorIndex_t;
  using NumericalMatrix_t        = typename Traits_t::NumericalMatrix_t;
  using NumericalMatrixIndex_t   = typename Traits_t::NumericalMatrixIndex_t;
  using DataStorage_t            = typename Traits_t::DataStorage_t;
  using DataStorageIndex_t       = typename Traits_t::DataStorageIndex_t;
  using ScaledDataStorage_t      = typename Traits_t::ScaledDataStorage_t;
  using ScaledDataStorageIndex_t = typename Traits_t::ScaledDataStorageIndex_t;
  using FullCoord_t              = typename Traits_t::FullCoord_t;
  using FullCoordIndex_t         = typename Traits_t::FullCoordIndex_t;
  using Covariates_t            = typename Traits_t::Covariates_t;
  using CovariatesIndex_t       = typename Traits_t::CovariatesIndex_t;

  using KernelPolicy_t    = KernelPolicy<LocalPolyKernelEstimator,    Traits_t>;
  using BandwidthPolicy_t = BandwidthPolicy<LocalPolyKernelEstimator, Traits_t>;
  using ScalingPolicy_t   = ScalingPolicy<LocalPolyKernelEstimator,   Traits_t>;

protected:

  // Policies
  // void KernelPolicy_t::Update () const
  // void BandwidthPolicy_t::Update () const
  // void ScalingPolicy_t::Update () const

public:

  // Policies
  using KernelPolicy_t::EvaluateKernel;          // Numeric_t EvaluateKernel (const Coord_t<Dim>&) const
  using BandwidthPolicy_t::GetKernelFactor;      // Numeric_t GetKernelFactor () const
  using BandwidthPolicy_t::GetKernelInput;       // Coord_t<Dim> GetKernelInput (const Coord_t<Dim>&) const
  using ScalingPolicy_t::ScaleCoordinate;        // Numeric_t ScaleCoordinate (const size_t&, const Numeric_t&) const
  using ScalingPolicy_t::InverseScaleCoordinate; // Numeric_t InverseScaleCoordinate (const size_t&, const Numeric_t&) const

public:

  // // // default constructor
  // // constexpr LocalPolyKernelEstimator () noexcept;
  // constexpr LocalPolyKernelEstimator () = default;
  // // // copy constructor
  // // LocalPolyKernelEstimator (const LocalPolyKernelEstimator<Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t> &other);
  // LocalPolyKernelEstimator (const LocalPolyKernelEstimator<Order, Dim, CoordTransformFunc_t, KernelFunc_t, Numeric_t>&) = default;
  // // // move constructor
  // // LocalPolyKernelEstimator (LocalPolyKernelEstimator<Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t> &&other);
  // LocalPolyKernelEstimator (LocalPolyKernelEstimator<Order, Dim, CoordTransformFunc_t, KernelFunc_t, Numeric_t> &&other) = default;
  // // // copy assignment operator (copy-swap)
  // // LocalPolyKernelEstimator<Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t>& operator= (LocalPolyKernelEstimator<Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t> other);
  // LocalPolyKernelEstimator<Order, Dim, CoordTransformFunc_t, KernelFunc_t, Numeric_t>& operator= (const LocalPolyKernelEstimator<Order, Dim, CoordTransformFunc_t, KernelFunc_t, Numeric_t> &other) = default;
  // // // move assignment operator
  // // LocalPolyKernelEstimator<Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t>& operator= (LocalPolyKernelEstimator<Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t> &&other);
  // LocalPolyKernelEstimator<Order, Dim, CoordTransformFunc_t, KernelFunc_t, Numeric_t>& operator= (LocalPolyKernelEstimator<Order, Dim, CoordTransformFunc_t, KernelFunc_t, Numeric_t> &&other) = default;

  // swap operation
  template<typename NumericType_, template<class, class> class KernelPolicy_, template<class, class> class BandwidthPolicy_, template<class, class> class ScalingPolicy_>
  friend void swap (LocalPolyKernelEstimator<NumericType_, KernelPolicy_, BandwidthPolicy_, ScalingPolicy_> &lhs,
                    LocalPolyKernelEstimator<NumericType_, KernelPolicy_, BandwidthPolicy_, ScalingPolicy_> &rhs) noexcept;


  // destructor
  virtual ~LocalPolyKernelEstimator () = default;

  // // compatible conversions
  // template <size_t OrderOther, typename NumericOther_t,
  //           // to prevent overwriting default copy constructor
  //           typename std::enable_if<!std::is_same<NumericOther_t, Numeric_t>::value, int>::type = 0>
  // LocalPolyKernelEstimator (const LocalPolyKernelEstimator<Dim, OrderOther, CoordTransformFunc_t, KernelFunc_t, NumericOther_t> &other);


  // Setters
  void UseWeightedData (const bool &enable = true);
  void UseFastInvert (const bool &enable = true) { m_fastInvert = enable; }
  void SetDimension (const size_t &dim);
  void SetOrder (const size_t &order);
  void SetTol (const Numeric_t &tol) { m_tol = tol; }
  void SetCovarianceMatrix (const NumericalMatrix_t &cov) {
    m_covarianceMatrix.ResizeTo(cov.GetNrows(), cov.GetNcols());
    m_covarianceMatrix = cov;
  }
  void ResetCovarianceMatrix () { m_covarianceMatrix.UnitMatrix(); }

  template<typename ... Args_t>
  typename std::enable_if<are_convertible<NumericType, Args_t ...>::value>::type
  SetDataPoint (const FullCoordIndex_t&,
                const Args_t               & ...);

  void SetDataPoint (const FullCoordIndex_t&,
                     const FullCoord_t     &) noexcept;

  template<typename ... Args_t>
  typename std::enable_if<are_convertible<NumericType, Args_t ...>::value>::type
  AddDataPoint (const Args_t& ...);

  void AddDataPoint (const FullCoord_t&);

  void SetNumberOfDataPoints (const DataStorageIndex_t &size) {
    this->m_coordinates.resize(size);
  }

  // Getters
  constexpr bool UsingWeightedData () const noexcept {
    return this->m_useWeightedData;
  }

  constexpr bool UsingFastInvert () const noexcept {
    return this->m_fastInvert;
  }

  constexpr size_t GetDimension () const noexcept {
    return this->m_dim;
  }

  constexpr size_t GetOrder () const noexcept {
    return this->m_order;
  }

  constexpr Numeric_t GetTol () const noexcept {
    return m_tol;
  }

  constexpr const NumericalMatrix_t& GetCovarianceMatrix () const noexcept {
    return m_covarianceMatrix;
  }

  constexpr const FullCoord_t& GetDataPoint (const DataStorageIndex_t &index) const noexcept {
    return this->m_coordinates[index];
  }

  FullCoord_t& GetDataPoint (const DataStorageIndex_t &index) noexcept {
    return this->m_coordinates[index];
  }

  Covariates_t GetCovariates (const DataStorageIndex_t &index) const {
    Covariates_t coords(m_dim);

    for (size_t idim = 0; idim < m_dim; ++idim) coords[idim] = m_coordinates[index][idim];
    return coords;
  } // GetCovariates

  constexpr const NumericType& GetResponse (const DataStorageIndex_t &index) const noexcept {
    return this->m_coordinates[index][m_dim];
  }

  NumericType& GetResponse (const DataStorageIndex_t &index) noexcept {
    return this->m_coordinates[index][m_dim];
  }

  constexpr const NumericType& GetResponseWeight (const DataStorageIndex_t &index) const noexcept {
    return this->m_coordinates[index][m_dim + 1];
  }

  NumericType& GetResponseWeight (const DataStorageIndex_t &index) noexcept {
    return this->m_coordinates[index][m_dim + 1];
  }

  constexpr DataStorageIndex_t GetNumberOfDataPoints () const noexcept {
    return m_coordinates.size();
  }

  // Evaluation
  template<typename ... Args_t>
  typename std::enable_if<are_convertible<NumericType, Args_t ...>::value, NumericType>::type
  operator()(const Args_t& ...) const;

  NumericType operator()(const Covariates_t&) const;

  // Utilities:
  // Compute smoothing matrix (useful in computing degrees of freedom)
  NumericalMatrix_t ComputeHhat () const;

protected:

  void RequestUpdate (const bool update_this = true) const;
  template<typename OtherNumericType_, template<class, class> class OtherKernelPolicy_, template<class, class> class OtherBandwidthPolicy_, template<class, class> class OtherScalingPolicy_>
  void CopyImpl (const LocalPolyKernelEstimator<OtherNumericType_, OtherKernelPolicy_, OtherBandwidthPolicy_, OtherScalingPolicy_> &other);

  // resizes matrices if necessary
  void              UpdateMatrices () const;
  void              CreateXMatrix () const;
  NumericalMatrix_t CreateWeightsVector () const;

  bool m_useWeightedData = false,                   // use unequally weighted data or not
       m_fastInvert      = false;                   // use fast matrix inversion
  mutable bool m_update  = true;                    // internal
  size_t m_dim           = 1,                       // covariate dimension
         m_order         = 0;                       // polynomial order of estimator
  Numeric_t m_tol        = 1e-104;// std::sqrt(std::numeric_limits<Numeric_t>::min());  // numerical zero
  NumericalMatrix_t           m_covarianceMatrix;   // store variance-covariance matrix
  mutable NumericalVector_t   m_weights,            // keep persistent only to avoid recreating each evaluation
                              m_backup_weights,     // keep persistent only to avoid recreating each evaluation
                              m_yVector;            // keep persistent only to avoid recreating each evaluation
  mutable ScaledDataStorage_t m_scaled_coordinates; // keep persistent only to avoid recreating each evaluation
  DataStorage_t               m_coordinates;        // n x m, where n = number of points, m = dimension + 1 + error
  mutable Covariates_t        m_scaled_coordinate;  // keep persistent only to avoid recreating each evaluation
  mutable NumericalMatrix_t   m_xMatrix,            // keep persistent only to avoid recreating each evaluation
                              m_invCovarianceMatrix;// keep persistent only to avoid recreating each evaluation
  mutable std::vector<bool>   m_goodDataPoint;      // keep persistent only to avoid recreating each evaluation

  ClassDef(LocalPolyKernelEstimator, 1);
};
}

// // TVectorT range-based for loop integration
// template <size_t Dim, size_t Order, typename CoordTransformFunc_t, typename KernelFunc_t, typename Numeric_t>
// auto begin (typename Analysis::LocalPolyKernelEstimator<Dim, Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t>::NumericalVector_t &vec)
//   ->decltype(begin(vec.GetMatrixArray()))
// {return begin(vec.GetMatrixArray());}
//
// template <size_t Dim, size_t Order, typename CoordTransformFunc_t, typename KernelFunc_t, typename Numeric_t>
// auto end (typename Analysis::LocalPolyKernelEstimator<Dim, Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t>::NumericalVector_t &vec)
//   ->decltype(end(vec.GetMatrixArray()))
// {return end(vec.GetMatrixArray());}
//
// template <size_t Dim, size_t Order, typename CoordTransformFunc_t, typename KernelFunc_t, typename Numeric_t>
// auto cbegin (const typename Analysis::LocalPolyKernelEstimator<Dim, Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t>::NumericalVector_t &vec)
//   ->decltype(cbegin(vec.GetMatrixArray()))
// {return cbegin(vec.GetMatrixArray());}
//
// template <size_t Dim, size_t Order, typename CoordTransformFunc_t, typename KernelFunc_t, typename Numeric_t>
// auto cend (const typename Analysis::LocalPolyKernelEstimator<Dim, Order, CoordTransformFunc_t, KernelFunc_t, Numeric_t>::NumericalVector_t &vec)
//   ->decltype(cend(vec.GetMatrixArray()))
// {return cend(vec.GetMatrixArray());}


// //////////////////////////////////////
#include "LocalPolyKernelEstimator.icc"
// //////////////////////////////////////

// templateClassImp(Analysis::LocalPolyKernelEstimatorTraits);
templateClassImp(Analysis::LocalPolyKernelEstimator);

#endif // include gaurds
