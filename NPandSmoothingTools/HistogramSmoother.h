// TODO: make this class properly const-correct

#ifndef HISTOGRAM_SMOOTHER_H
#define HISTOGRAM_SMOOTHER_H

#include <array>
#include <vector>
#include <memory>

#include "TF1.h"
#include "TAxis.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph.h"
#include "TGraphErrors.h"
#include "TGraph2D.h"
#include "TGraph2DErrors.h"


#include "NPandSmoothingTools/LocalPolyKernelEstimator.h"

namespace Analysis {

template <typename NumericType, template <class, class> class KernelPolicy, template <class, class> class BandwidthPolicy, template <class, class> class ScalingPolicy>
class HistogramSmoother :
  public LocalPolyKernelEstimator<NumericType, KernelPolicy, BandwidthPolicy, ScalingPolicy>
{
public:
  using LPKE = LocalPolyKernelEstimator<NumericType, KernelPolicy, BandwidthPolicy, ScalingPolicy>;
  using LPKE::GetDimension;
  using BinStorage = std::vector<size_t>;
  using BinStorageIndex = typename BinStorage::size_type;
  // Policies
  using LPKE::ScalingPolicy_t::ScaleCoordinate;
  using LPKE::ScalingPolicy_t::InverseScaleCoordinate;

  // destructor
  virtual ~HistogramSmoother () = default;

  // Setters
  void SetNbins (const BinStorageIndex &index, size_t bins) {if (m_nBins.size() != GetDimension()) m_nBins.resize(GetDimension()); m_nBins[index] = bins;}

  // Getters
  size_t GetNbins (const BinStorageIndex &index) const {return m_nBins[index];}

  void LoadData (const TH1 &h, bool useBinErrors = true);

  TH1* MakeSmoothedTH1 () const;

  // TODO: put these into SmoothingUtilities
  // Simple chi^2/Ndf (least-squares)
  // double  GetChi2overNdof (int bias = 0, bool smooth_error = true);
  // double  GetOptimalChi2NdfWidthX (double delta = 0.0, double target = 1.0, int bias = 0);

protected:
  std::vector<double>       MakeScaledVector (size_t dim, double min, double max) const;
  void                      CopyAxisInfo (TAxis *from, TAxis *to) const;

  BinStorage            m_nBins = {100};
  std::shared_ptr<TH1>  m_axisHisto;

  ClassDef(HistogramSmoother, 1);
};

}

// //////////////////////////////////////
#include "HistogramSmoother.icc"
// //////////////////////////////////////

templateClassImp(Analysis::HistogramSmoother)

#endif // include gaurds
