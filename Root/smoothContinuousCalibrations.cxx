// NOTE: Weighted data is NOT being used anymore. The covariance matrix is now
// being used. Currently, only the diagonal elements are filled (stat. errors).

// File     : smoothContinuousCalibrations.C
// Author   : Jeffrey Wayne Hetherly
// Purpose  : Pre-processing macro to smooth continuous calibrations


#include <tuple>
#include <string>
#include <unordered_map>
#include <map>
#include <vector>
#include <set>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TMath.h"
#include "TDirectoryFile.h"
#include "TFile.h"
#include "TList.h"
#include "TMap.h"
#include "TKey.h"
#include "TF1.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TString.h"
#include "TObjString.h"
#include "TRegexp.h"
#include "TVectorD.h"

#include "NPandSmoothingTools/SmoothingUtils.h"
#include "NPandSmoothingTools/DynamicHistogramSmoother.h"
// #include "NPandSmoothingTools/HistogramSmoother.h"
// #include "NPandSmoothingTools/BandwidthPolicies.h"
// #include "NPandSmoothingTools/ScalingPolicies.h"
// #include "NPandSmoothingTools/KernelPolicies.h"


using namespace ROOT;
using std::string;
using std::unordered_map;
using std::map;
using std::vector;
using std::set;


# include "CalibrationDataInterface/CalibrationDataContainer.h"

using Analysis::CalibrationDataContainer;
using Analysis::CalibrationDataFunctionContainer;
using Analysis::CalibrationDataHistogramContainer;

// locally visible code
namespace {

// // this can change within the smoother later
// static const size_t KERNEL_DIM = 1;

// dynamic
using KernelSmoother = Analysis::DynamicHistogramSmoother;
// using KernelSmoother = Analysis::HistogramSmoother<double,
//                                                    Analysis::DynamicKernelPolicy,
//                                                    Analysis::GlobalBandwidthPolicy,
//                                                    Analysis::DynamicScalingPolicy>;

// static
// using KernelSmoother = Analysis::HistogramSmoother<double,
//       Analysis::NormalKernelPolicy,
//       Analysis::GlobalBandwidthPolicy,
//       Analysis::LnScalingPolicy>;

string kernel = "Normal";
/*
  Uniform
  Triangular
  Epanechnikov
  Biweight
  Triweight
  Tricube
  Normal
  Logistic
  Dummy
 * */

auto kernels = std::make_tuple(
    [](const typename KernelSmoother::Covariates_t &c) // UniformKernel
    {
    typename KernelSmoother::Numeric_t normalInput = 0;

    for (auto coord : c) normalInput += coord * coord;
    auto r = TMath::Sqrt(normalInput);
    return typename KernelSmoother::Numeric_t(r > 1.0 ? 0.0 : 0.5);
    },
    [](const typename KernelSmoother::Covariates_t &c) // TriangularKernel
    {
    typename KernelSmoother::Numeric_t normalInput = 0;

    for (auto coord : c) normalInput += coord * coord;
    auto r = TMath::Sqrt(normalInput);
    return typename KernelSmoother::Numeric_t(r > 1.0 ? 0.0 : 1.0 - TMath::Abs(r));
    },
    [](const typename KernelSmoother::Covariates_t &c) // EpanechnikovKernel
    {
    typename KernelSmoother::Numeric_t normalInput = 0;

    for (auto coord : c) normalInput += coord * coord;
    auto r = TMath::Sqrt(normalInput);
    return typename KernelSmoother::Numeric_t(r > 1.0 ? 0.0 : 0.75*(1.0 - r*r));
    },
    [](const typename KernelSmoother::Covariates_t &c) // BiweightKernel
    {
      typename KernelSmoother::Numeric_t normalInput = 0;

      for (auto coord : c) normalInput += coord * coord;
      auto r = TMath::Sqrt(normalInput);
      auto term = 1 - r*r;
      return typename KernelSmoother::Numeric_t(r > 1.0 ? 0.0 : 0.9375*(term*term));
    },
    [](const typename KernelSmoother::Covariates_t &c) // TriweightKernel
    {
      typename KernelSmoother::Numeric_t normalInput = 0;

      for (auto coord : c) normalInput += coord * coord;
      auto r = TMath::Sqrt(normalInput);
      auto term = 1 - r*r;
      return typename KernelSmoother::Numeric_t(r > 1.0 ? 0.0 : 1.09375*(term*term*term));
    },
    [](const typename KernelSmoother::Covariates_t &c) // TricubeKernel
    {
      typename KernelSmoother::Numeric_t normalInput = 0;

      for (auto coord : c) normalInput += coord * coord;
      auto r = TMath::Sqrt(normalInput);
      auto term = 1 - TMath::Abs(r*r*r);
      return typename KernelSmoother::Numeric_t(r > 1.0 ? 0.0 : 0.8641975308642*(term*term*term));
    },
    [](const typename KernelSmoother::Covariates_t &c) // NormalKernel
    {
      typename KernelSmoother::Numeric_t normalInput = 0;

      for (auto coord : c) normalInput += coord * coord;
      return typename KernelSmoother::Numeric_t(TMath::Gaus(TMath::Sqrt(normalInput), 0.0, 1.0, kTRUE));
    },
    [](const typename KernelSmoother::Covariates_t &c) // LogisticKernel
    {
      typename KernelSmoother::Numeric_t normalInput = 0;

      for (auto coord : c) normalInput += coord * coord;
      auto r = TMath::Sqrt(normalInput);
      return typename KernelSmoother::Numeric_t(1.0/(TMath::Exp(r) + 2.0 + TMath::Exp(-r)));
    },
    [](const typename KernelSmoother::Covariates_t&) // DummyKernel
    {
      return 1.0;
    }
    );

auto LnScale = [](const typename KernelSmoother::Numeric_t &v)
               { return typename KernelSmoother::Numeric_t(TMath::Log(v)); };

auto LnInvScale = [](const typename KernelSmoother::Numeric_t &v)
                  { return typename KernelSmoother::Numeric_t(TMath::Exp(v)); };

auto NoScale = [](const typename KernelSmoother::Numeric_t &v)
               {return v;};


// "directory/name within file" versus TObject* objects
typedef map<TString, TObject*> data_map_t;

struct FlavorInfo {
  string fileTag;
  bool   smooth;
  int    loose, medium, tight;
};

// math
double ComputePValue (double chi2, double ndf)
{
  return 1.0 - TMath::Gamma(ndf*0.5, chi2*0.5)/TMath::Gamma(ndf*0.5);
}

// histogram helpers
auto GetTotalNbins = [](TH1 &h)->decltype(h.GetNbinsX()) {
  auto ndim     = h.GetDimension();
  auto nbins    = h.GetNbinsX() + 2;

  if (ndim > 1) nbins *= (h.GetNbinsY() + 2);
  if (ndim > 2) nbins *= (h.GetNbinsZ() + 2);
  return nbins;
};

auto GetTotalRealNbins = [](TH1 &h)->decltype(h.GetNbinsX()) {
  auto nbins = 0;
  for (Int_t bin = 0; bin < GetTotalNbins(h); ++bin)
    if (!h.IsBinOverflow(bin) && !h.IsBinUnderflow(bin)) ++nbins;
  return nbins;
};

auto GetRealBinIndex = [](TH1 &h, Int_t index)->decltype(h.GetNbinsX()) {
  auto bin_counter = -1;
  for (Int_t bin = 0; bin < GetTotalNbins(h); ++bin) {
    if (!h.IsBinOverflow(bin) && !h.IsBinUnderflow(bin)) {
      if (bin_counter < 0) bin_counter = 0;
      if (bin == index) break;
      ++bin_counter;
    }
  }

  return bin_counter;
};

auto RemoveHistogramError = [](TH1 &h)->void {
  auto nbins    = GetTotalNbins(h);

  for (decltype(nbins) bin = 0; bin < nbins; ++bin)
    h.SetBinError(bin, 0.0);
};

// auto AddErrorQuadToFrom = [](TH1 &q, const TH1 &a)->void {
//   auto nbins    = GetTotalNbins(q);
//
//   for (decltype(nbins) bin = 0; bin < nbins; ++bin) {
//     auto prev = q.GetBinError(bin),
//          add = a.GetBinError(bin),
//          quad = TMath::Sqrt(prev*prev + add*add);
//     q.SetBinError(bin, quad);
//   }
// };

auto AddQuadToFrom = [](TH1 &q, const TH1 &a)->void {
  auto nbins    = GetTotalNbins(q);

  for (decltype(nbins) bin = 0; bin < nbins; ++bin) {
    auto prev = q.GetBinContent(bin),
         add = a.GetBinContent(bin),
         quad = TMath::Sqrt(prev*prev + add*add);
    q.SetBinContent(bin, quad);
  }
};

double ComputeChiSquared (TH1 &original, const KernelSmoother &smoother, const KernelSmoother::NumericalMatrix_t &inv_cov)
{
  std::vector<double> diffs(GetTotalRealNbins(original));
  double result = 0.0;
  auto dim = original.GetDimension();
  Int_t binx, biny, binz;

  for (Int_t bin = 0; bin < GetTotalNbins(original); ++bin) {
    if (original.IsBinOverflow(bin) || original.IsBinUnderflow(bin)) continue;
    auto real_bin = GetRealBinIndex(original, bin);
    diffs[real_bin] = 0.0;
    auto content = original.GetBinContent(bin);
    original.GetBinXYZ(bin, binx, biny, binz);
    if (dim == 1) {
      auto xcenter = original.GetXaxis()->GetBinCenter(binx);
      diffs[real_bin] = content - smoother(xcenter);
    } else if (dim == 2) {
      auto xcenter = original.GetXaxis()->GetBinCenter(binx),
           ycenter = original.GetYaxis()->GetBinCenter(biny);
      diffs[real_bin] = content - smoother(xcenter, ycenter);
    } else if (dim == 3) {
      auto xcenter = original.GetXaxis()->GetBinCenter(binx),
           ycenter = original.GetYaxis()->GetBinCenter(biny),
           zcenter = original.GetZaxis()->GetBinCenter(binz);
      diffs[real_bin] = content - smoother(xcenter, ycenter, zcenter);
    }
  }
  for (KernelSmoother::NumericalMatrixIndex_t row = 0; row < inv_cov.GetNrows(); ++row)
    for (KernelSmoother::NumericalMatrixIndex_t col = row; col < inv_cov.GetNcols(); ++col) {
      auto term = diffs[row]*inv_cov(row, col)*diffs[col];
      if (row == col) result += term;
      else result += 2.0*term;
    }

  return result;
}

void HarvestHistograms (TFile *f, TFile *of, data_map_t &m)
{
  TList *key_list = f->GetListOfKeys();

  // local functor
  struct local_
  {
    TFile *_f;
    TFile *_of;

    void operator()(TString map_index, data_map_t &m, TList *keys)
    {
      for (int i = 0; i < keys->GetSize(); ++i) {
        TKey   *key = dynamic_cast<TKey*>(keys->At(i));
        TString new_map_index(map_index),
        class_name(key->GetClassName());

        if (class_name == "TDirectoryFile") {
          new_map_index += key->GetName();
          new_map_index += "/";
          TDirectoryFile *dir_file = dynamic_cast<TDirectoryFile*>(key->ReadObj());
          (*this)(new_map_index,
                  m,
                  dir_file->GetListOfKeys());
        }
        else if (class_name == "Analysis::CalibrationDataHistogramContainer") {
          _of->mkdir(new_map_index.Data());
          new_map_index   += key->GetName();
          m[new_map_index] = key->ReadObj();
          CalibrationDataHistogramContainer *c = dynamic_cast<CalibrationDataHistogramContainer*>(m[new_map_index]);

          // loop through histograms
          TIterator  *it = c->MakeIterator();
          TObjString *k;

          while ((k = (TObjString*)it->Next())) {
            // if ((*c)(k)->InheritsFrom("TH2")) {
            //   TH2 *h = (TH2*)(*c)(k);
            //   h->SetDirectory(0);
            // }

            /*else*/ if ((*c)(k)->InheritsFrom("TH1")) {
              TH1 *h = (TH1*)(*c)(k);
              h->SetDirectory(0);
            }
          }
        }
        else if (key->InheritsFrom("TH1")) {
          _of->mkdir(new_map_index.Data());
          new_map_index   += key->GetName();
          m[new_map_index] = key->ReadObj();
          TH1 *h = dynamic_cast<TH1*>(key->ReadObj());
          h->SetDirectory(0);
        }
        else {
          _of->mkdir(new_map_index.Data());
          new_map_index   += key->GetName();
          m[new_map_index] = key->ReadObj();
        }
      }
    } // ()
  } Harvester;

  Harvester._f  = f;
  Harvester._of = of;
  Harvester("", m, key_list);
} // HarvestHistograms

bool SmoothHistogram (TH1             &unsmoothed,
                      TH1             &eff,
                      KernelSmoother  &smoother,
                      const bool      &uniaxis = true,
                      const bool      &is_nominal = true,
                      const KernelSmoother::NumericalMatrix_t *cov = nullptr,
                      const bool      &ignore_tw_relationship = false)
{
  // check for dimension here
  auto ndim  = unsmoothed.GetDimension();
  auto nbins = GetTotalNbins(unsmoothed);

  // ///////////////////////////
  // 3D Histogram
  // ///////////////////////////
  if (ndim == 3) {
    auto        *h        = static_cast<TH3*>(&unsmoothed);
    auto        *e        = static_cast<TH3*>(&eff);
    TAxis       *pt_axis  = nullptr,
                *eta_axis = nullptr,
                *tw_axis  = nullptr;
    bool smooth_errors = smoother.UsingWeightedData(),
         eta_is_x      = false,
         eta_is_y      = false,
         tw_is_x      = false,
         tw_is_y      = false;

    // get axis information
    if (TString(h->GetXaxis()->GetTitle()) == "pt") {
      // pt_axis    = h->GetXaxis();
      pt_axis    = e->GetXaxis();
    }
    else if (TString(h->GetYaxis()->GetTitle()) == "pt") {
      // pt_axis    = h->GetYaxis();
      pt_axis    = e->GetYaxis();
    }
    else {
      // pt_axis    = h->GetZaxis();
      pt_axis    = e->GetZaxis();
    }
    if (TString(h->GetXaxis()->GetTitle()) == "eta" ||
        TString(h->GetXaxis()->GetTitle()) == "abseta") {
      eta_axis    = h->GetXaxis();
      eta_is_x    = true;
    }
    else if (TString(h->GetYaxis()->GetTitle()) == "eta" ||
             TString(h->GetYaxis()->GetTitle()) == "abseta") {
      eta_axis    = h->GetYaxis();
      eta_is_y    = true;
    }
    else {
      eta_axis    = h->GetZaxis();
    }
    if (TString(h->GetXaxis()->GetTitle()) == "tagweight") {
      tw_axis    = h->GetXaxis();
      tw_is_x    = true;
    }
    else if (TString(h->GetYaxis()->GetTitle()) == "tagweight") {
      tw_axis    = h->GetYaxis();
      tw_is_y    = true;
    }
    else {
      tw_axis    = h->GetZaxis();
    }

    for (int i = 0; i <= nbins; ++i) {
      if (!h->IsBinOverflow(i) && !h->IsBinUnderflow(i) &&
          (h->GetBinError(i) == 0.0)) {
        smoother.UseWeightedData(false);
        break;
      }
    }

    if (uniaxis) { // smooth only along pT axis
      // NOTE: eta bin, tw bin, pt histo
      map<int, map<int, TH1*>>    smoothed_histos;
      auto                        last_tw_bin = tw_axis->GetNbins();

      // Smooth along pt separately for each eta bin
      for (int i = 1; i <= eta_axis->GetNbins(); ++i) {
        for (int j = 1; j <= last_tw_bin; ++j) {
          TH1 *to_smooth = nullptr;
          if (eta_is_x && tw_is_y) to_smooth = h->ProjectionZ("_pz", i, i, j, j); // assumes pt is z
          else if (eta_is_y && tw_is_x) to_smooth = h->ProjectionZ("_pz", j, j, i, i); // assumes pt is z
          else if (eta_is_x && !tw_is_x && !tw_is_y) to_smooth = h->ProjectionY("_py", i, i, j, j); // assumes pt is y
          else if (!eta_is_x && !eta_is_y && tw_is_x) to_smooth = h->ProjectionY("_py", j, j, i, i); // assumes pt is y
          else if (eta_is_y && !tw_is_x && !tw_is_y) to_smooth = h->ProjectionX("_px", i, i, j, j); // assumes pt is x
          else if (!eta_is_y && !eta_is_x && tw_is_y) to_smooth = h->ProjectionX("_px", j, j, i, i); // assumes pt is x

          to_smooth->SetDirectory(0);

          smoother.LoadData(*to_smooth, false);

          Info("SmoothHistogram", "at %s bin %d and %s bin %d", eta_axis->GetTitle(), i, tw_axis->GetTitle(), j);
          smoothed_histos[i][j] = static_cast<TH1*>(smoother.MakeSmoothedTH1());

          delete to_smooth;
        }
      }

      {
        auto &hist = smoothed_histos[1][1];
        // Rebin histogram
        if (eta_is_x) {
          if (tw_is_y) {
            h->SetBins(h->GetNbinsX(), h->GetXaxis()->GetXbins()->GetArray(),
                       h->GetNbinsY(), h->GetYaxis()->GetXbins()->GetArray(),
                       hist->GetNbinsX(), hist->GetXaxis()->GetXbins()->GetArray());
          }
          else {
            h->SetBins(h->GetNbinsX(), h->GetXaxis()->GetXbins()->GetArray(),
                       hist->GetNbinsX(), hist->GetXaxis()->GetXbins()->GetArray(),
                       h->GetNbinsZ(), h->GetZaxis()->GetXbins()->GetArray());
          }
        }
        else if (eta_is_y) {
          if (tw_is_x) {
            h->SetBins(h->GetNbinsX(), h->GetXaxis()->GetXbins()->GetArray(),
                       h->GetNbinsY(), h->GetYaxis()->GetXbins()->GetArray(),
                       hist->GetNbinsX(), hist->GetXaxis()->GetXbins()->GetArray());
          }
          else {
            h->SetBins(hist->GetNbinsX(), hist->GetXaxis()->GetXbins()->GetArray(),
                       h->GetNbinsY(), h->GetYaxis()->GetXbins()->GetArray(),
                       h->GetNbinsZ(), h->GetZaxis()->GetXbins()->GetArray());
          }
        }
        else {
          if (tw_is_x) {
            h->SetBins(h->GetNbinsX(), h->GetXaxis()->GetXbins()->GetArray(),
                       hist->GetNbinsX(), hist->GetXaxis()->GetXbins()->GetArray(),
                       h->GetNbinsZ(), h->GetZaxis()->GetXbins()->GetArray());
          }
          else {
            h->SetBins(hist->GetNbinsX(), hist->GetXaxis()->GetXbins()->GetArray(),
                       h->GetNbinsY(), h->GetYaxis()->GetXbins()->GetArray(),
                       h->GetNbinsZ(), h->GetZaxis()->GetXbins()->GetArray());
          }
        }
      }

      // std::cout << "JWH: smoothed histogram bins" << std::endl;
      // std::cout << "JWH: # of " << h->GetXaxis()->GetTitle() << " bins = " << h->GetXaxis()->GetNbins() << std::endl;
      // std::cout << "JWH: # of " << h->GetYaxis()->GetTitle() << " bins = " << h->GetYaxis()->GetNbins() << std::endl;
      // std::cout << "JWH: # of " << h->GetZaxis()->GetTitle() << " bins = " << h->GetZaxis()->GetNbins() << std::endl;
      //
      // std::cout << "JWH: efficiency bins" << std::endl;
      // std::cout << "JWH: # of " << eff.GetXaxis()->GetTitle() << " bins = " << eff.GetXaxis()->GetNbins() << std::endl;
      // std::cout << "JWH: # of " << eff.GetYaxis()->GetTitle() << " bins = " << eff.GetYaxis()->GetNbins() << std::endl;
      // std::cout << "JWH: # of " << eff.GetZaxis()->GetTitle() << " bins = " << eff.GetZaxis()->GetNbins() << std::endl;


      {
        auto last_tw_bins = static_cast<TH1*>(smoothed_histos[1][1]->Clone());
        last_tw_bins->SetDirectory(0);

        for (auto &eta_pair : smoothed_histos) {
          int eta_bin = eta_pair.first;
          double eta_value =  eta_is_x ? h->GetXaxis()->GetBinCenter(eta_bin) :
                              eta_is_y ? h->GetYaxis()->GetBinCenter(eta_bin) :
                                         h->GetZaxis()->GetBinCenter(eta_bin);
          last_tw_bins->Reset();

          for (auto &tw_pair : eta_pair.second) {
            int tw_bin = tw_pair.first;
            double tw_value = tw_is_x ? h->GetXaxis()->GetBinCenter(tw_bin) :
                              tw_is_y ? h->GetYaxis()->GetBinCenter(tw_bin) :
                                        h->GetZaxis()->GetBinCenter(tw_bin);
            TH1 *hist = tw_pair.second; // smoothed pt bins
            TH1 *eff_smooth = nullptr;
            {
              auto i = eta_bin,
                   j = tw_bin;
              if (eta_is_x && tw_is_y) eff_smooth = static_cast<TH3&>(eff).ProjectionZ("_pz", i, i, j, j); // assumes pt is z
              else if (eta_is_y && tw_is_x) eff_smooth = static_cast<TH3&>(eff).ProjectionZ("_pz", j, j, i, i); // assumes pt is z
              else if (eta_is_x && !tw_is_x && !tw_is_y) eff_smooth = static_cast<TH3&>(eff).ProjectionY("_py", i, i, j, j); // assumes pt is y
              else if (!eta_is_x && !eta_is_y && tw_is_x) eff_smooth = static_cast<TH3&>(eff).ProjectionY("_py", j, j, i, i); // assumes pt is y
              else if (eta_is_y && !tw_is_x && !tw_is_y) eff_smooth = static_cast<TH3&>(eff).ProjectionX("_px", i, i, j, j); // assumes pt is x
              else if (!eta_is_y && !eta_is_x && tw_is_y) eff_smooth = static_cast<TH3&>(eff).ProjectionX("_px", j, j, i, i); // assumes pt is x

              eff_smooth->SetDirectory(0);
              smoother.LoadData(*eff_smooth, false);
              delete eff_smooth;
              eff_smooth = static_cast<TH1*>(smoother.MakeSmoothedTH1());
              eff_smooth->SetDirectory(0);
            }

            // Actually fill bins
            for (int pt_bin = 1; pt_bin <= hist->GetNbinsX(); ++pt_bin) {
              int global_bin = eta_is_x && !tw_is_x && !tw_is_y ? h->GetBin(eta_bin, pt_bin, tw_bin) :
                               eta_is_x && tw_is_y ? h->GetBin(eta_bin, tw_bin, pt_bin) :
                               eta_is_y && tw_is_x ? h->GetBin(tw_bin, eta_bin, pt_bin) :
                               eta_is_y && !tw_is_x && !tw_is_y ? h->GetBin(pt_bin, eta_bin, tw_bin) :
                               tw_is_x ? h->GetBin(tw_bin, pt_bin, eta_bin) :
                               /*tw_is_y ?*/ h->GetBin(pt_bin, tw_bin, eta_bin);

              // NOTE: do tag weight axis restriction
              double pt_value = hist->GetXaxis()->GetBinCenter(pt_bin);
              int unsmoothed_pt_bin = !eta_is_x && !tw_is_x ? eff.GetXaxis()->FindBin(pt_value) :
                                      !eta_is_y && !tw_is_y ? eff.GetYaxis()->FindBin(pt_value) :
                                      eff.GetZaxis()->FindBin(pt_value);
              int local_global_bin = eta_is_x && !tw_is_x && !tw_is_y ? eff.GetBin(eta_bin, unsmoothed_pt_bin, tw_bin) :
                                     eta_is_x && tw_is_y ? eff.GetBin(eta_bin, tw_bin, unsmoothed_pt_bin) :
                                     eta_is_y && tw_is_x ? eff.GetBin(tw_bin, eta_bin, unsmoothed_pt_bin) :
                                     eta_is_y && !tw_is_x && !tw_is_y ? eff.GetBin(unsmoothed_pt_bin, eta_bin, tw_bin) :
                                     tw_is_x ? eff.GetBin(tw_bin, unsmoothed_pt_bin, eta_bin) :
                                     /*tw_is_y ?*/ eff.GetBin(unsmoothed_pt_bin, tw_bin, eta_bin);
              // double local_eff = eff.GetBinContent(local_global_bin);
              double local_eff = eff_smooth->GetBinContent(pt_bin);
              if (tw_bin == tw_axis->GetNbins() && is_nominal) {
                if (!ignore_tw_relationship) {
                  h->SetBinContent(global_bin, (1.0 - last_tw_bins->GetBinContent(pt_bin))/local_eff);
                  h->SetBinError(global_bin, last_tw_bins->GetBinError(pt_bin)/local_eff);
                }
                else {
                  h->SetBinContent(global_bin, hist->GetBinContent(pt_bin));
                  h->SetBinError(global_bin, hist->GetBinError(pt_bin));
                }
              }
              else {
                last_tw_bins->SetBinContent(pt_bin, last_tw_bins->GetBinContent(pt_bin) + hist->GetBinContent(pt_bin)*local_eff);
                last_tw_bins->SetBinError(pt_bin, TMath::Sqrt(last_tw_bins->GetBinError(pt_bin)*last_tw_bins->GetBinError(pt_bin) + local_eff*local_eff*hist->GetBinError(pt_bin)*hist->GetBinError(pt_bin)));

                h->SetBinContent(global_bin, hist->GetBinContent(pt_bin));
                h->SetBinError(global_bin, hist->GetBinError(pt_bin));
              }
              // if (is_nominal) {
              //   std::cout << "JWH: pt bin " << pt_bin << ", eta bin " << eta_bin << ", tagweight bin " << tw_bin << ": " << h->GetBinContent(global_bin) << " +/- " << h->GetBinError(global_bin) << std::endl;
              //   std::cout << "JWH: eff = " << local_eff << std::endl;
              // }
            }

            delete eff_smooth;
          }
        }

        delete last_tw_bins;
      }


      smoother.UseWeightedData(smooth_errors);
      return true;
    }
    else {
      smoother.LoadData(*h, false);
      if (cov) smoother.SetCovarianceMatrix(*cov);

      auto smoothed = smoother.MakeSmoothedTH1();
      smoothed->SetDirectory(0);

      h->SetBins(smoothed->GetNbinsX(), smoothed->GetXaxis()->GetXbins()->GetArray(),
                 smoothed->GetNbinsY(), smoothed->GetYaxis()->GetXbins()->GetArray(),
                 smoothed->GetNbinsZ(), smoothed->GetZaxis()->GetXbins()->GetArray());

      // Actually fill bins
      for (int xbin = 1; xbin <= smoothed->GetNbinsX(); ++xbin) {
        for (int ybin = 1; ybin <= smoothed->GetNbinsY(); ++ybin) {
          for (int zbin = 1; zbin <= smoothed->GetNbinsZ(); ++zbin) {
            h->SetBinContent(xbin, ybin, zbin, smoothed->GetBinContent(xbin, ybin, zbin));
            h->SetBinError(xbin, ybin, zbin, smoothed->GetBinError(xbin, ybin, zbin));
          }
        }
      }

      delete smoothed;
      // TODO: (maybe not needed)
      // remove "all-zero rows" in other axis bins
    }

    return true;
  }

  // ///////////////////////////
  // 2D Histogram
  // ///////////////////////////
  else if (ndim == 2) {
    auto h             = static_cast<TH2*>(&unsmoothed);
    bool x_is_pt       = false,
         smooth_errors = smoother.UsingWeightedData();
    TAxis *pt_axis     = nullptr,
          *other_axis  = nullptr;

    // get axis information
    if (TString(h->GetXaxis()->GetTitle()) == "pt") {
      pt_axis    = h->GetXaxis();
      other_axis = h->GetYaxis();
      x_is_pt    = true;
    }
    else if (TString(h->GetYaxis()->GetTitle()) == "pt") {
      pt_axis    = h->GetYaxis();
      other_axis = h->GetXaxis();
    }
    else if (uniaxis) {
      Error("SmoothHistogram", "unable to find \"pt\" axis in %s - required for smoothing procedure", h->GetName());
      return false;
    }

    for (int i = 0; i <= nbins; ++i) {
      if (!h->IsBinOverflow(i) && !h->IsBinUnderflow(i) &&
          (h->GetBinError(i) == 0.0)) smoother.UseWeightedData(false);
    }

    if (uniaxis) { // smooth only along pT axis
      vector<TH1*> smoothed_histos;
      KernelSmoother::NumericalMatrix_t similarity(other_axis->GetNbins()*pt_axis->GetNbins(), pt_axis->GetNbins());
      KernelSmoother::NumericalMatrix_t reduced_cov(pt_axis->GetNbins(), pt_axis->GetNbins());

      // vector<TH1*>  smoothed_histos(other_axis->GetNbins());

      // Smooth along pt separately for each "other" bin
      for (int i = 1; i <= other_axis->GetNbins(); ++i) {
        TH1D to_smooth("to smooth", Form("%s_%s_bin_%d", h->GetName(), other_axis->GetTitle(), i),
                       pt_axis->GetNbins(), pt_axis->GetXbins()->GetArray());

        // TH1D to_smooth(x_is_pt ? *h->ProjectionY(Form("%s_%s_bin_%d", h->GetName(), other_axis->GetTitle(), i), i, i) : *h->ProjectionX(Form("%s_%s_bin_%d", h->GetName(), other_axis->GetTitle(), i), i, i));
        to_smooth.SetDirectory(0);

        // Fill smoothed histogram
        auto column_counter = 0;
        for (int j = 1; j <= pt_axis->GetNbins(); ++j) {
          int global_bin = x_is_pt ? h->GetBin(j, i) : h->GetBin(i, j);
          auto real_bin = GetRealBinIndex(unsmoothed, global_bin);

          for (Int_t row = 0; row < similarity.GetNrows(); ++row) {
            if (row == real_bin) similarity(row, column_counter) = 1;
            else similarity(row, column_counter) = 0;
          }

          to_smooth.SetBinContent(j, h->GetBinContent(global_bin));
          to_smooth.SetBinError(j, h->GetBinError(global_bin));

          // TODO: until CDI nan bug is fixed
          if (h->GetBinError(global_bin) != h->GetBinError(global_bin)) {
            int global_bin_previous = x_is_pt ? h->GetBin(j - 1, i) : h->GetBin(i, j - 1);
            to_smooth.SetBinError(j, h->GetBinError(global_bin_previous));
          }
          ++column_counter;
        }

        if (cov) {
          for (Int_t row = 0; row < reduced_cov.GetNrows(); ++row)
            for (Int_t col = 0; col < reduced_cov.GetNcols(); ++col) {
              reduced_cov(row, col) = 0;
              for (Int_t sim_row1 = 0; sim_row1 < similarity.GetNrows(); ++sim_row1)
                for (Int_t sim_row2 = 0; sim_row2 < similarity.GetNrows(); ++sim_row2)
                  reduced_cov(row, col) += similarity(sim_row2, row)*(*cov)(sim_row2, sim_row1)*similarity(sim_row1, col);
            }
        }

        smoother.LoadData(to_smooth, false);
        if (cov) smoother.SetCovarianceMatrix(reduced_cov);
        // TODO: testing: bandwidth selection
        if (cov) {
          reduced_cov.Invert();
          auto h_prev = smoother.GetBandwidth(0),
               h_optimal = 0.0,
               chired_optimal = 10.0;
          int npoints = 100;
          double h_initial = 0.05,
                 h_delta = 0.05;
          double *hs = new double[npoints],
                 *chi2 = new double[npoints],
                 *effDoF = new double[npoints],
                 *chi2red = new double[npoints],
                 *prob = new double[npoints];
          for (int ibandwidth = 0; ibandwidth < npoints; ++ibandwidth) {
            hs[ibandwidth] = h_initial + (npoints - ibandwidth)*h_delta;
            smoother.SetBandwidth(0, hs[ibandwidth]);
            chi2[ibandwidth] = ComputeChiSquared(to_smooth, smoother, reduced_cov);

            // NOTE: retrieve regression eff. DoF
            auto H = smoother.ComputeHhat();
            // NOTE: convert h to residual eff. DoF (H -> I - H)
            for (Int_t ii = 0; ii < H.GetNrows(); ++ii) {
              for (Int_t jj = ii; jj < H.GetNcols(); ++jj) {
                if (ii == jj) H(ii, ii) = 1.0 - H(ii, ii);
                else {
                  H(ii, jj) = -H(ii, jj);
                  H(jj, ii) = -H(jj, ii);
                }
              }
            }
            auto HH = H*H.T();
            double trace_HH = 0.0,
                   trace_H  = 0.0,
                   trace_2H_HH  = 0.0;

            for (Int_t ii = 0; ii < HH.GetNrows(); ++ii) {
              trace_HH += HH(ii, ii);
              trace_H += H(ii, ii);
            }
            trace_2H_HH = 2.0*trace_H - trace_HH;

            // effDoF[ibandwidth] = trace_H;
            effDoF[ibandwidth] = trace_HH;
            // effDoF[ibandwidth] = trace_2H_HH;
            chi2red[ibandwidth] = chi2[ibandwidth]/effDoF[ibandwidth];
            prob[ibandwidth] = 1.0 - TMath::Prob(chi2[ibandwidth], effDoF[ibandwidth]);
            if (chired_optimal > TMath::Abs(1.0 - chi2red[ibandwidth])) {
              h_optimal = hs[ibandwidth];
              chired_optimal = TMath::Abs(1.0 - chi2red[ibandwidth]);
            }
          }

          // static int counter = 0;
          // auto f = TFile::Open(Form("bandwidth_selection_%d.root", counter), "RECREATE");
          // f->cd();
          // auto c1 = new TCanvas(Form("c1_%d", counter), "#chi^{2} and Effective Residual Degrees of Freedom (effDoF)", 600, 400);
          // ++counter;
          // c1->cd();
          //
          // // auto gr1 = new TGraph(npoints, hs, chi2);
          // auto gr1 = new TGraph(npoints, hs, chi2red);
          // // auto gr1 = new TGraph(npoints, hs, prob);
          // gr1->SetTitle("#chi^{2} and Effective Residual Degrees of Freedom (effDoF)");
          // gr1->SetMarkerColor(kBlue);
          // gr1->SetMarkerStyle(21);
          // gr1->GetXaxis()->SetTitle("bandwidth");
          // // gr1->GetYaxis()->SetTitle("1 - p-value");
          // gr1->GetYaxis()->SetTitle("#chi^{2}_{red}");
          // // gr1->GetYaxis()->SetTitle("#chi^{2}");
          // gr1->GetYaxis()->SetTitleSize(0.07);
          // gr1->Draw();
          // c1->Update();
          //
          // // auto gr2 = new TGraph(npoints, hs, effDoF);
          // // //scale gr2 to the pad coordinates
          // // auto rightmax = 1.1*gr2->GetYaxis()->GetXmax(),
          // //      rightmin = 0.9*gr2->GetYaxis()->GetXmin();
          // // auto fa = new TF1("fa", "([3] - [2])*x/([1] - [0]) + [2]", -5, 1000);
          // // fa->SetParameter(0, rightmin);
          // // fa->SetParameter(1, rightmax);
          // // fa->SetParameter(2, gr1->GetYaxis()->GetXmin());
          // // fa->SetParameter(3, gr1->GetYaxis()->GetXmax());
          // // gr2->Apply(fa);
          // // // auto scale = gPad->GetUymax()/rightmax;
          // // gr2->SetMarkerColor(kRed);
          // // gr2->SetMarkerStyle(20);
          // // // gr2->GetYaxis()->SetTitle("effDoF");
          // // // gr2->Scale(scale);
          // // gr2->Draw("same");
          // // c1->Update();
          // //
          // // // draw an axis on the right side
          // // TGaxis *axis = new TGaxis(gr1->GetXaxis()->GetXmax(), gr1->GetYaxis()->GetXmin(),
          // //                           gr1->GetXaxis()->GetXmax(), gr1->GetYaxis()->GetXmax(), rightmin, rightmax, 510, "+L");
          // // axis->SetTitle("effDoF");
          // // axis->SetLineColor(kRed);
          // // axis->SetLabelColor(kRed);
          // // axis->Draw();
          // // c1->Update();
          //
          // c1->Write();
          // f->Close();
          // // delete fa;

          delete []hs;
          delete []chi2;
          delete []effDoF;
          delete []chi2red;
          delete []prob;

          if (h_prev > 0.0) smoother.SetBandwidth(0, h_prev);
          else {
            Info("SmoothHistogram", "setting \"optimal\" bandwidth at %f", h_optimal);
            smoother.SetBandwidth(0, h_optimal);
          }
        }

        Info("SmoothHistogram", "at %s bin %d", other_axis->GetTitle(), i);
        smoothed_histos.push_back(smoother.MakeSmoothedTH1());
      }

      // Rebin and fill original 2D histogram
      int other_bin = 1;

      for (vector<TH1*>::iterator hist = smoothed_histos.begin();
           hist != smoothed_histos.end(); ++hist) {
        // Rebin histogram
        if (x_is_pt && (other_bin == 1)) {
          h->SetBins((*hist)->GetNbinsX(), (*hist)->GetXaxis()->GetXbins()->GetArray(),
                     h->GetNbinsY(), h->GetYaxis()->GetXbins()->GetArray());
        }
        else if (other_bin == 1) {
          h->SetBins(h->GetNbinsX(), h->GetXaxis()->GetXbins()->GetArray(),
                     (*hist)->GetNbinsX(), (*hist)->GetXaxis()->GetXbins()->GetArray());
        }

        // Actually fill bins
        for (int pt_bin = 1; pt_bin <= (*hist)->GetNbinsX(); ++pt_bin) {
          int global_bin = x_is_pt ? h->GetBin(pt_bin, other_bin) :
                           h->GetBin(other_bin, pt_bin);
          h->SetBinContent(global_bin, (*hist)->GetBinContent(pt_bin));
          h->SetBinError(global_bin, (*hist)->GetBinError(pt_bin));
        }

        ++other_bin;
        (*hist)->SetDirectory(0);
        delete *hist;
      }
      smoother.UseWeightedData(smooth_errors);
      return true;
    }
    else { // full 2D smoothing
      smoother.LoadData(*h, false);
      if (cov) smoother.SetCovarianceMatrix(*cov);

      auto smoothed = smoother.MakeSmoothedTH1();
      smoothed->SetDirectory(0);

      h->SetBins(smoothed->GetNbinsX(), smoothed->GetXaxis()->GetXbins()->GetArray(),
          smoothed->GetNbinsY(), smoothed->GetYaxis()->GetXbins()->GetArray());

      // Actually fill bins
      for (int xbin = 1; xbin <= smoothed->GetNbinsX(); ++xbin) {
        for (int ybin = 1; ybin <= smoothed->GetNbinsY(); ++ybin) {
          h->SetBinContent(xbin, ybin, smoothed->GetBinContent(xbin, ybin));
          h->SetBinError(xbin, ybin, smoothed->GetBinError(xbin, ybin));
        }
      }

      delete smoothed;
      // TODO: (maybe not needed)
      // remove "all-zero rows" in other axis bins
      return true;
    }
  }

  // ///////////////////////////
  // 1D Histogram
  // ///////////////////////////
  else if (ndim == 1) {
    // TH1   *h          = &unsmoothed;
    // TAxis *pt_axis    = h->GetXaxis();
  }

  return false;
} // SmoothHistogram

TH1* GetNuisanceVariation (KernelSmoother &smoother,

                           // const TH1 &result, const TH1 &unsmoothed_result,
                           // TH1 &variation, const int &bin = -1,
                           TH1 &variation, TH1 &eff, const int &bin = -1,
                           const bool &uniaxis = true,
                           const bool &ignore_tw_relationship = false)
{
  auto var = bin >= 0 ? static_cast<TH1*>(variation.Clone()) : &variation;

  if (bin >= 0) {
    var->Reset();
    auto ndim  = var->GetDimension();
    auto nbins = var->GetNbinsX() + 2;

    if (ndim > 1) nbins *= (var->GetNbinsY() + 2);

    if (ndim > 2) nbins *= (var->GetNbinsZ() + 2);

    for (auto i = 0; i <= nbins; ++i) {
      if (var->IsBinOverflow(i) || var->IsBinUnderflow(i)) continue;

      if (i == bin) var->SetBinContent(bin, variation.GetBinContent(bin));
      var->SetBinError(i, variation.GetBinError(i));
    }
  }

  // NOTE: in general this the correct way
  //       however, the addition/subtraction introduces numerical instability
  //       and the functional form of a kernel smoother makes this unnecessary
  // var->Add(&unsmoothed_result);
  if (!SmoothHistogram(*var, eff, smoother, uniaxis, false, nullptr, ignore_tw_relationship)) {
    Error("GetNuisanceVariation", "unable to smooth nuisance variation %s", var->GetName());
    return nullptr;
  }

  // var->Add(&result, -1.0);

  return var;
} // GetNuisanceVariation

void SetSmootherFor3D (KernelSmoother &smoother, Int_t order, Int_t bins0, Int_t bins1, Int_t bins2,
                       double bandwidth0, double bandwidth1, double bandwidth2)
{
  smoother.SetDimension(3);
  smoother.SetOrder(order);
  smoother.SetNbins(0, bins0);
  smoother.SetNbins(1, bins1);
  smoother.SetNbins(2, bins2);
  smoother.SetBandwidth(0, bandwidth0);
  smoother.SetBandwidth(1, bandwidth1);
  smoother.SetBandwidth(2, bandwidth2);
  if (kernel == "Uniform")      smoother.SetKernel(std::get<0>(kernels));
  if (kernel == "Triangular")   smoother.SetKernel(std::get<1>(kernels));
  if (kernel == "Epanechnikov") smoother.SetKernel(std::get<2>(kernels));
  if (kernel == "Biweight")     smoother.SetKernel(std::get<3>(kernels));
  if (kernel == "Triweight")    smoother.SetKernel(std::get<4>(kernels));
  if (kernel == "Tricube")      smoother.SetKernel(std::get<5>(kernels));
  if (kernel == "Normal")       smoother.SetKernel(std::get<6>(kernels));
  if (kernel == "Logistic")     smoother.SetKernel(std::get<7>(kernels));
  if (kernel == "Dummy")        smoother.SetKernel(std::get<8>(kernels));

  // assumes pt, eta, tagweight
  smoother.SetScaleFunction(0, LnScale);
  smoother.SetInvScaleFunction(0, LnInvScale);
  smoother.SetScaleFunction(1, NoScale);
  smoother.SetInvScaleFunction(1, NoScale);
  smoother.SetScaleFunction(2, NoScale);
  smoother.SetInvScaleFunction(2, NoScale);
}

void SetSmootherFor2D (KernelSmoother &smoother, Int_t order, Int_t bins0, Int_t bins1, double bandwidth0, double bandwidth1)
{
  smoother.SetDimension(2);
  smoother.SetOrder(order);
  smoother.SetNbins(0, bins0);
  smoother.SetNbins(1, bins1);
  smoother.SetBandwidth(0, bandwidth0);
  smoother.SetBandwidth(1, bandwidth1);
  if (kernel == "Uniform")      smoother.SetKernel(std::get<0>(kernels));
  if (kernel == "Triangular")   smoother.SetKernel(std::get<1>(kernels));
  if (kernel == "Epanechnikov") smoother.SetKernel(std::get<2>(kernels));
  if (kernel == "Biweight")     smoother.SetKernel(std::get<3>(kernels));
  if (kernel == "Triweight")    smoother.SetKernel(std::get<4>(kernels));
  if (kernel == "Tricube")      smoother.SetKernel(std::get<5>(kernels));
  if (kernel == "Normal")       smoother.SetKernel(std::get<6>(kernels));
  if (kernel == "Logistic")     smoother.SetKernel(std::get<7>(kernels));
  if (kernel == "Dummy")        smoother.SetKernel(std::get<8>(kernels));

  // assumes eta, pt or tagweight, pt (Project3D("xz"))
  smoother.SetScaleFunction(0, NoScale);
  smoother.SetInvScaleFunction(0, NoScale);
  smoother.SetScaleFunction(1, LnScale);
  smoother.SetInvScaleFunction(1, LnInvScale);
}

void SetSmootherFor1D (KernelSmoother &smoother, Int_t order, Int_t bins, double bandwidth)
{
  smoother.SetDimension(1);
  smoother.SetOrder(order);
  smoother.SetNbins(0, bins);
  smoother.SetBandwidth(0, bandwidth);
  if (kernel == "Uniform")      smoother.SetKernel(std::get<0>(kernels));
  if (kernel == "Triangular")   smoother.SetKernel(std::get<1>(kernels));
  if (kernel == "Epanechnikov") smoother.SetKernel(std::get<2>(kernels));
  if (kernel == "Biweight")     smoother.SetKernel(std::get<3>(kernels));
  if (kernel == "Triweight")    smoother.SetKernel(std::get<4>(kernels));
  if (kernel == "Tricube")      smoother.SetKernel(std::get<5>(kernels));
  if (kernel == "Normal")       smoother.SetKernel(std::get<6>(kernels));
  if (kernel == "Logistic")     smoother.SetKernel(std::get<7>(kernels));
  if (kernel == "Dummy")        smoother.SetKernel(std::get<8>(kernels));

  smoother.SetScaleFunction(0, LnScale);
  smoother.SetInvScaleFunction(0, LnInvScale);
}


}

void Analysis::smoothContinuousCalibrations (TString  fileName,
                                             TString  container,
                                             std::vector<int> bNPset,
                                             std::vector<int> cNPset,
                                             std::vector<int> lightNPset,
                                             unsigned pt_bins,
                                             float    b_pt_smoothing,
                                             float    c_pt_smoothing,
                                             float    light_pt_smoothing,
                                             size_t   order,
                                             bool     smoothB,
                                             bool     smoothC,
                                             bool     smoothLight)
{
  TString output_tag = TString::Format("_order_%zu_smoothing_%g_%g_%g_ptbins_%u", order, b_pt_smoothing, c_pt_smoothing, light_pt_smoothing, pt_bins),
  new_file_name(fileName);
  KernelSmoother smoother;
  set<string>    toSmooth;
  map<string, FlavorInfo> flavorInfo;

  bool smoothTau = smoothC;

  flavorInfo["B"] = {
    "B",
    smoothB,
    bNPset[0], bNPset[1], bNPset[2]
  };
  flavorInfo["C"] = {
    "C",
    smoothC,
    cNPset[0], cNPset[1], cNPset[2]
  };
  flavorInfo["Tau"] = {
    "T",
    smoothTau,
    cNPset[0], cNPset[1], cNPset[2]
  };
  flavorInfo["Light"] = {
    "Light",
    smoothLight,
    lightNPset[0], lightNPset[1], lightNPset[2]
  };

  // toSmooth.insert("default_SF");
  // toSmooth.insert("IntegratedWP_5TagBins_SF");
  // toSmooth.insert("IntegratedWP_7TagBins_SF");
  toSmooth.insert(container.Data());

  // "directory/name within file" versus TObject objects
  data_map_t data_containers;

  TFile *f = TFile::Open(fileName.Data(), "READ");

  if (!f) {
    Fatal("smoothContinuousCalibrations", "couldn't open CDI file %s - aborting...", fileName.Data());
    return;
  }

  // create output file
  TRegexp reg("[.]root*");

  if (!fileName.Contains(output_tag))
    Info("smoothContinuousCalibrations", "storing smoothed calibration data in %s", fileName.Insert(fileName.Index(reg), output_tag).Data());
  else
    Info("smoothContinuousCalibrations", "storing smoothed calibration data in %s", fileName.Data());
  TFile *of = TFile::Open(fileName.Data(), "RECREATE");

  // spider through file and harvest histograms
  HarvestHistograms(f, of, data_containers);

  f->Close();


  // TODO:
  // loop through histograms and find optimal width here


  // loop through all calibration data containers
  for (data_map_t::iterator it = data_containers.begin(); it != data_containers.end(); ++it) {
    TObject *o = it->second;
    TString  full_name(it->first),
             directory_name(TString(full_name).Remove(full_name.Last('/'))),
             class_name(o->ClassName()),
             key_name(TString(full_name).Remove(0, full_name.Last('/') + 1));

    auto shouldSmooth = [&full_name, &class_name, &toSmooth, &flavorInfo](string tag) -> bool {
                          bool result = false;

                          for (auto key : toSmooth) {
                            result = result || (class_name == "Analysis::CalibrationDataHistogramContainer" &&
                                                full_name.Contains(TRegexp(("/" + flavorInfo[tag].fileTag + "/" + key).c_str())) &&
                                                // full_name.EndsWith(("/" + flavorInfo[tag].fileTag + "/" + key).c_str()) &&
                                                flavorInfo[tag].smooth);
                          }
                          return result;
                        };

    // smooth histograms if necessary
    if (shouldSmooth("B") || shouldSmooth("C") || shouldSmooth("Tau") || shouldSmooth("Light")) {
      auto eff_full_name                    = full_name;
      auto uniaxis                          = false;
      auto ignore_tw_relationship           = false;
      auto *c = static_cast<CalibrationDataHistogramContainer*>(o);
      auto *eff = static_cast<CalibrationDataHistogramContainer*>(data_containers[eff_full_name.ReplaceAll("_SF", "_Eff").Data()]);
      TList                              syst_nuis_names;
      TIterator                         *itt = c->MakeIterator();
      TObjString                        *k;
      TH1                               &result               = *static_cast<TH1*>((*c)("result")),
                                        &total_systematics    = *static_cast<TH1*>((*c)("systematics")),
                                        &efficiency           = *static_cast<TH1*>((*eff)("result"));

      // if (directory_name.Contains("continuous")) {
      if (directory_name.Contains("Continuous") || directory_name.Contains("continuous")) {
        if (shouldSmooth("Light")) {
          // pt, eta, tagweight
          // uniaxis = false;
          // SetSmootherFor3D(smoother, order,
          //     pt_bins, 5, 20,
          //     0.4, 0.5, 0.25);
          uniaxis = true;
          // ignore_tw_relationship = true;
          SetSmootherFor1D(smoother, order,
              pt_bins,
              light_pt_smoothing == -1.0 ? 0.4 : light_pt_smoothing);
        } else if (shouldSmooth("C") || shouldSmooth("Tau")) {
          // pt, eta, tagweight
          // uniaxis = false;
          // SetSmootherFor3D(smoother, order,
          //     pt_bins, 5, 20,
          //     0.4, 0.5, 0.25);
          uniaxis = true;
          // ignore_tw_relationship = true;
          SetSmootherFor1D(smoother, order,
              pt_bins,
              c_pt_smoothing == -1.0 ? 0.4 : c_pt_smoothing);
        } else if (shouldSmooth("B")) {
          // pt vs. tagweight
          uniaxis = true;
          SetSmootherFor1D(smoother, order,
              pt_bins,
              b_pt_smoothing == -1.0 ? 0.4 : b_pt_smoothing);
        }
      }

      if ((TString(result.GetXaxis()->GetTitle()) == "pt" && result.GetXaxis()->GetNbins() <= order && uniaxis) ||
          (TString(result.GetYaxis()->GetTitle()) == "pt" && result.GetYaxis()->GetNbins() <= order && uniaxis) ||
          (TString(result.GetZaxis()->GetTitle()) == "pt" && result.GetZaxis()->GetNbins() <= order && uniaxis)) {

        int nbins = (result.GetXaxis()->GetNbins() > 0 ? result.GetXaxis()->GetNbins() : 1)*
                    (result.GetYaxis()->GetNbins() > 0 ? result.GetYaxis()->GetNbins() : 1)*
                    (result.GetZaxis()->GetNbins() > 0 ? result.GetZaxis()->GetNbins() : 1);
        medium_too_few_pt_points = nbins;
        tight_too_few_pt_points = nbins;

        // set reduction list
        if (full_name.Contains(("/" + flavorInfo["B"].fileTag + "/").c_str())) {
          TVectorD *reductionSets = new TVectorD(2);
          (*reductionSets)[0] = tight_too_few_pt_points;
          (*reductionSets)[1] = medium_too_few_pt_points;
          TObjString *str = new TObjString("ReducedSets");
          c->Add(str, reductionSets);
        }
        else if (full_name.Contains(("/" + flavorInfo["C"].fileTag + "/").c_str())) {
          TVectorD *reductionSets = new TVectorD(2);
          (*reductionSets)[0] = tight_too_few_pt_points;
          (*reductionSets)[1] = medium_too_few_pt_points;
          TObjString *str = new TObjString("ReducedSets");
          c->Add(str, reductionSets);
        }
        else if (full_name.Contains(("/" + flavorInfo["Tau"].fileTag + "/").c_str())) {
          TVectorD *reductionSets = new TVectorD(2);
          (*reductionSets)[0] = tight_too_few_pt_points;
          (*reductionSets)[1] = medium_too_few_pt_points;
          TObjString *str = new TObjString("ReducedSets");
          c->Add(str, reductionSets);
        }
        else if (full_name.Contains(("/" + flavorInfo["Light"].fileTag + "/").c_str())) {
          TVectorD *reductionSets = new TVectorD(2);
          (*reductionSets)[0] = tight_too_few_pt_points;
          (*reductionSets)[1] = medium_too_few_pt_points;
          TObjString *str = new TObjString("ReducedSets");
          c->Add(str, reductionSets);
        }

        Warning("smoothContinuousCalibrations", "calibration has too few bins to smooth!");

        continue;
      }

      TH1                               &result_no_error      = *static_cast<TH1*>((*c)("result")->Clone()),
                                        &unsmoothed_result    = *static_cast<TH1*>((*c)("result")->Clone());

      result_no_error.SetDirectory(0);
      unsmoothed_result.SetDirectory(0);

      smoother.UseWeightedData(false);
      if (!SmoothHistogram(total_systematics, efficiency, smoother, uniaxis, false, nullptr, ignore_tw_relationship)) {
        Fatal("smoothContinuousCalibrations", "unable to smooth \"systematics\" histogram - aborting...");
        return;
      }
      total_systematics.Reset();

      // rebin extrapolation in regions that kinematically overlap with all other systematics
      if (c->Contains("extrapolation")) {
        Info("smoothContinuousCalibrations", "modifying extrapolation uncertainties to match number of bins in systematic \"overlap\" region");
        auto &extrapolation     = *static_cast<TH1*>((*c)("extrapolation")),
             &ref_extrapolation = *static_cast<TH1*>((*c)("extrapolation")->Clone());
        ref_extrapolation.SetDirectory(0);
        auto ndim = ref_extrapolation.GetDimension();
        auto xmin_new = 0.0,
             xmax_new = 0.0,
             ymin_new = 0.0,
             ymax_new = 0.0,
             zmin_new = 0.0,
             zmax_new = 0.0;

        for (auto dim = 1; dim <= ndim; ++dim) {
          double new_xmin = 0.0, new_xmax = 0.0;
          vector<Double_t> new_bins;
          TAxis *perm_axis = nullptr,
                *new_axis = nullptr;

          if (dim == 1) {
            perm_axis = extrapolation.GetXaxis();
            new_axis = total_systematics.GetXaxis();
            xmin_new = new_axis->GetXmin();
            xmax_new = new_axis->GetXmax();
          }
          if (dim == 2) {
            perm_axis = extrapolation.GetYaxis();
            new_axis = total_systematics.GetYaxis();
            ymin_new = new_axis->GetXmin();
            ymax_new = new_axis->GetXmax();
          }
          if (dim == 3) {
            perm_axis = extrapolation.GetZaxis();
            new_axis = total_systematics.GetZaxis();
            zmin_new = new_axis->GetXmin();
            zmax_new = new_axis->GetXmax();
          }

          new_xmin = new_axis->GetXmin();
          new_xmax = new_axis->GetXmax();

          if (perm_axis->GetNbins() == new_axis->GetNbins()) continue;

          bool filled_with_new = false;
          for (auto perm_bin = 0; perm_bin <= perm_axis->GetNbins(); ++perm_bin) {
            auto perm_value = (*perm_axis->GetXbins())[perm_bin];
            // the additional tolerance protection is due to small numerical differences in bin ranges
            if (perm_value - new_xmin < -smoother.GetTol() || perm_value - new_xmax > smoother.GetTol())
              new_bins.push_back(perm_value);
            else if (!filled_with_new) {
              for (auto new_bin = 0; new_bin <= new_axis->GetNbins(); ++new_bin)
                new_bins.push_back((*new_axis->GetXbins())[new_bin]);
              filled_with_new = true;
            }
          }

          Double_t *raw_new_bins = new Double_t[new_bins.size()];
          auto i = 0;
          for (auto val : new_bins)
            raw_new_bins[i++] = val;
          perm_axis->Set(new_bins.size() - 1, raw_new_bins);

        }

        if (ndim == 1) {
          extrapolation.SetBins(extrapolation.GetXaxis()->GetNbins(), extrapolation.GetXaxis()->GetXbins()->GetArray());
          for (auto xbin = 1; xbin <= extrapolation.GetNbinsX(); ++xbin) {
            auto gbin = extrapolation.GetBin(xbin);
            auto xcenter = extrapolation.GetXaxis()->GetBinCenter(xbin);
            bool x_is_zero = ((total_systematics.GetXaxis()->FindFixBin(xcenter) <= total_systematics.GetNbinsX()) &&
                              (total_systematics.GetXaxis()->FindFixBin(xcenter) >= 1));
            if (x_is_zero) extrapolation.SetBinContent(gbin, 0.0);
            else extrapolation.SetBinContent(gbin, ref_extrapolation.GetBinContent(ref_extrapolation.FindFixBin(xcenter)));
            extrapolation.SetBinError(gbin, extrapolation.GetBinContent(gbin));
          }
        }
        else if (ndim == 2) {
          extrapolation.SetBins(extrapolation.GetXaxis()->GetNbins(), extrapolation.GetXaxis()->GetXbins()->GetArray(),
                                extrapolation.GetYaxis()->GetNbins(), extrapolation.GetYaxis()->GetXbins()->GetArray());
          for (auto xbin = 1; xbin <= extrapolation.GetNbinsX(); ++xbin) {
            for (auto ybin = 1; ybin <= extrapolation.GetNbinsY(); ++ybin) {
              auto gbin = extrapolation.GetBin(xbin, ybin);
              auto xcenter = extrapolation.GetXaxis()->GetBinCenter(xbin),
                   ycenter = extrapolation.GetYaxis()->GetBinCenter(ybin);
              bool x_is_zero = ((total_systematics.GetXaxis()->FindFixBin(xcenter) <= total_systematics.GetNbinsX()) &&
                                (total_systematics.GetXaxis()->FindFixBin(xcenter) >= 1)),
                   y_is_zero = ((total_systematics.GetYaxis()->FindFixBin(ycenter) <= total_systematics.GetNbinsY()) &&
                                (total_systematics.GetYaxis()->FindFixBin(ycenter) >= 1));
              if (x_is_zero && y_is_zero) extrapolation.SetBinContent(gbin, 0.0);
              else extrapolation.SetBinContent(gbin, ref_extrapolation.GetBinContent(ref_extrapolation.FindFixBin(xcenter, ycenter)));
              extrapolation.SetBinError(gbin, extrapolation.GetBinContent(gbin));
            }
          }
        }
        else if (ndim == 3) {
          extrapolation.SetBins(extrapolation.GetXaxis()->GetNbins(), extrapolation.GetXaxis()->GetXbins()->GetArray(),
                                extrapolation.GetYaxis()->GetNbins(), extrapolation.GetYaxis()->GetXbins()->GetArray(),
                                extrapolation.GetZaxis()->GetNbins(), extrapolation.GetZaxis()->GetXbins()->GetArray());
          for (auto xbin = 1; xbin <= extrapolation.GetNbinsX(); ++xbin) {
            for (auto ybin = 1; ybin <= extrapolation.GetNbinsY(); ++ybin) {
              for (auto zbin = 1; zbin <= extrapolation.GetNbinsZ(); ++zbin) {
                auto gbin = extrapolation.GetBin(xbin, ybin, zbin);
                auto xcenter = extrapolation.GetXaxis()->GetBinCenter(xbin),
                     ycenter = extrapolation.GetYaxis()->GetBinCenter(ybin),
                     zcenter = extrapolation.GetZaxis()->GetBinCenter(zbin);
                bool x_is_zero = ((total_systematics.GetXaxis()->FindFixBin(xcenter) <= total_systematics.GetNbinsX()) &&
                                  (total_systematics.GetXaxis()->FindFixBin(xcenter) >= 1)),
                     y_is_zero = ((total_systematics.GetYaxis()->FindFixBin(ycenter) <= total_systematics.GetNbinsY()) &&
                                  (total_systematics.GetYaxis()->FindFixBin(ycenter) >= 1)),
                     z_is_zero = ((total_systematics.GetZaxis()->FindFixBin(zcenter) <= total_systematics.GetNbinsZ()) &&
                                  (total_systematics.GetZaxis()->FindFixBin(zcenter) >= 1));
                if (x_is_zero && y_is_zero && z_is_zero) extrapolation.SetBinContent(gbin, 0.0);
                else extrapolation.SetBinContent(gbin, ref_extrapolation.GetBinContent(ref_extrapolation.FindFixBin(xcenter, ycenter, zcenter)));
                extrapolation.SetBinError(gbin, extrapolation.GetBinContent(gbin));
              }
            }
          }
        }

        delete &ref_extrapolation;
      }
      else
        Warning("smoothContinuousCalibrations", "couldn't find extrapolation uncertainty");

      // TODO: need to replace these with standalone methods
      // auto optimal_bandwidth = GetOptimalChi2NdfWidthX(0.05, 1, 1);
      // auto L1CV_optimal = GetOptimalLCVWidthX(0.05);
      // cout << "Optimal Bandwidth: " << optimal_bandwidth << "\t\t" << L1CV_optimal << endl;
      // smoother.SetBandwidth(0, optimal_bandwidth < 0.4 ? 0.4 : optimal_bandwidth);

      // calculate varaince-covariance matrix
      KernelSmoother::NumericalMatrix_t cov(GetTotalRealNbins(result), GetTotalRealNbins(result));
      {
        Info("smoothContinuousCalibrations", "calculating variance-covariance matrix");
        KernelSmoother::NumericalMatrix_t cor_cov(GetTotalRealNbins(result), GetTotalRealNbins(result));
        KernelSmoother::NumericalMatrix_t uncor_cov(GetTotalRealNbins(result), GetTotalRealNbins(result));
        uncor_cov.Zero();

        // loop through all systematic histograms to find bin-correlated and -uncorrelated ones
        itt = c->MakeIterator();
        while ((k = (TObjString*)itt->Next())) {
          if ((*c)(k)->InheritsFrom("TH1")) {
            TH1 &h = *static_cast<TH1*>((*c)(k));

            if ((k->GetString() == "result") ||
                (k->GetString() == "comment") ||
                (k->GetString() == "combined") ||
                (k->GetString() == "MCreference") ||
                (k->GetString() == "MChadronisation") ||
                (k->GetString() == "ReducedSets") ||
                (k->GetString() == "extrapolation") ||
                (k->GetString() == "statistics") ||
                (k->GetString() == "systematics")) continue;

            auto isCorrelated = c->isBinCorrelated(k->GetString().Data());
            auto covIndex = 0;
            for (Int_t bin = 0; bin < GetTotalNbins(h); ++bin) {
              if (h.IsBinOverflow(bin) || h.IsBinUnderflow(bin)) continue;
              auto average_error = (h.GetBinContent(bin) + h.GetBinError(bin))*.5; // average of up/down variation
              if (isCorrelated) {
                auto corcovIndex = 0;
                for (Int_t corbin = 0; corbin < GetTotalNbins(h); ++corbin) {
                  if (h.IsBinOverflow(corbin) || h.IsBinUnderflow(corbin)) continue;
                  auto coraverage_error = (h.GetBinContent(corbin) + h.GetBinError(corbin))*.5; // average of up/down variation
                  cor_cov(covIndex, corcovIndex) += average_error*coraverage_error;
                  ++corcovIndex;
                }
              }
              else
                uncor_cov(covIndex, covIndex) += average_error*average_error;
              ++covIndex;
            }
          }
        }
        {
          auto covIndex = 0;
          for (Int_t bin = 0; bin < GetTotalNbins(result); ++bin) {
            if (result.IsBinOverflow(bin) || result.IsBinUnderflow(bin)) continue;
            auto average_error = result.GetBinError(bin);
            uncor_cov(covIndex, covIndex) += average_error*average_error;
            ++covIndex;
          }
        }
        for (KernelSmoother::NumericalMatrixIndex_t row = 0; row < cov.GetNrows(); ++row)
          for (KernelSmoother::NumericalMatrixIndex_t col = row; col < cov.GetNcols(); ++col) {
            cov(row, col) = cor_cov(row, col) + uncor_cov(row, col);
            // TODO: remove whenever new way is developed
            // cov(row, col) = uncor_cov(row, col); // OLD WAY
            if (row != col) cov(col, row) = cov(row, col);
          }
      }
      for (decltype(smoother.GetDimension()) dim = 0; dim < smoother.GetDimension(); ++dim)
        Info("smoothContinuousCalibrations", "smoothing using a global bandwidth %f", smoother.GetBandwidth(dim));

      // smooth central value "result"
      smoother.UseWeightedData(false);
      Info("smoothContinuousCalibrations", "smoothing %s/%s", directory_name.Data(), key_name.Data());
      Info("smoothContinuousCalibrations", "\t\tusing bin errors - %d", smoother.UsingWeightedData());
      Info("smoothContinuousCalibrations", "\t\thistogram name, title - %s, %s", result.GetName(), result.GetTitle());
      // for (decltype(smoother.GetDimension()) dim = 0; dim < smoother.GetDimension(); ++dim) {
      for (decltype(result.GetDimension()) dim = 0; dim < result.GetDimension(); ++dim) {
        // if (smoother.GetDimension() == 1) ++dim; // Hack
        const char *title;
        auto naxis_bins = 0;
        if (dim == 0) {
          title = result.GetXaxis()->GetTitle();
          naxis_bins = result.GetXaxis()->GetNbins();
        }
        else if (dim == 1) {
          title = result.GetYaxis()->GetTitle();
          naxis_bins = result.GetYaxis()->GetNbins();
        }
        else if (dim == 2) {
          title = result.GetZaxis()->GetTitle();
          naxis_bins = result.GetZaxis()->GetNbins();
        }
        // if (smoother.GetDimension() == 1) --dim; // Hack
        // Info("smoothContinuousCalibrations", "\t\tnumber of smoothed %s bins - %lu", title, smoother.GetNbins(dim));
        Info("smoothContinuousCalibrations", "\t\tnumber of %s bins - %d", title, naxis_bins);
      }

      if (!SmoothHistogram(result, efficiency, smoother, uniaxis, true, &cov, ignore_tw_relationship)) {
        Fatal("smoothContinuousCalibrations", "unable to smooth \"result\" histogram - aborting...");
        return;
      }

      smoother.UseWeightedData(false);
      Info("smoothContinuousCalibrations", "smoothing %s/%s", directory_name.Data(), key_name.Data());
      Info("smoothContinuousCalibrations", "\t\tusing bin errors - %d", smoother.UsingWeightedData());
      Info("smoothContinuousCalibrations", "\t\thistogram name, title - %s, %s", result.GetName(), result.GetTitle());
      // for (decltype(smoother.GetDimension()) dim = 0; dim < smoother.GetDimension(); ++dim) {
      for (decltype(result.GetDimension()) dim = 0; dim < result.GetDimension(); ++dim) {
        // // Hack
        // if (smoother.GetDimension() == 1) ++dim;
        const char *title;
        auto naxis_bins = 0;
        if (dim == 0) {
          title = result.GetXaxis()->GetTitle();
          naxis_bins = result.GetXaxis()->GetNbins();
        }
        else if (dim == 1) {
          title = result.GetYaxis()->GetTitle();
          naxis_bins = result.GetYaxis()->GetNbins();
        }
        else if (dim == 2) {
          title = result.GetZaxis()->GetTitle();
          naxis_bins = result.GetZaxis()->GetNbins();
        }
        // // Hack
        // if (smoother.GetDimension() == 1) --dim;
        Info("smoothContinuousCalibrations", "\t\tnumber of smoothed %s bins - %d", title, naxis_bins);
      }

      if (!SmoothHistogram(result_no_error, efficiency, smoother, uniaxis, true, nullptr, ignore_tw_relationship)) {
        Fatal("smoothContinuousCalibrations", "unable to smooth \"result\" (no bin errors) histogram - aborting...");
        return;
      }

      // this is needed for CDI interface
      RemoveHistogramError(result);
      RemoveHistogramError(result_no_error);

      // perform stat. nuisance treatment of errors
      {
        TH1 *stat_var = static_cast<TH1*>(unsmoothed_result.Clone());
        int  nbins    = GetTotalNbins(*stat_var);

        for (int bin = 0; bin < nbins; ++bin)
          stat_var->SetBinContent(bin, stat_var->GetBinError(bin));

        stat_var->SetDirectory(0);

        // smoother.UseWeightedData(true);
        smoother.UseWeightedData(false);
        auto count       = 0,
             total_count = 0;

        for (int bin = 0; bin <= nbins; ++bin) {
          if (stat_var->IsBinOverflow(bin) || stat_var->IsBinUnderflow(bin)) continue;
          ++total_count;
        }

        for (int bin = 0; bin <= nbins; ++bin) {
          if (stat_var->IsBinOverflow(bin) || stat_var->IsBinUnderflow(bin)) continue;
          ++count;
          TObjString *new_key = new TObjString(Form("FT_EFF_Stat Nuis %d of %d", count, total_count));

          //Info("smoothContinuousCalibrations", "creating %s/%s", directory_name.Data(), Form("stat_nuis_%d_of_%d", count, total_count));
          Info("smoothContinuousCalibrations", "creating %s/%s", directory_name.Data(), new_key->GetString().Data());
          Info("smoothContinuousCalibrations", "\t\tusing bin errors - %d", smoother.UsingWeightedData());
          Info("smoothContinuousCalibrations", "\t\thistogram name, title - %s, %s", new_key->GetString().Data(), Form("stat_nuis_%d_of_%d", count, total_count));
          // for (decltype(smoother.GetDimension()) dim = 0; dim < smoother.GetDimension(); ++dim) {
          //   // Hack
          //   if (smoother.GetDimension() == 1) ++dim;
          //   const char *title;
          //   if (dim == 0)
          //     title = result.GetXaxis()->GetTitle();
          //   else if (dim == 1)
          //     title = result.GetYaxis()->GetTitle();
          //   else if (dim == 2)
          //     title = result.GetZaxis()->GetTitle();
          //   // Hack
          //   if (smoother.GetDimension() == 1) --dim;
          //   Info("smoothContinuousCalibrations", "\t\tnumber of smoothed %s bins - %lu", title, smoother.GetNbins(dim));
          // }
          TH1 *stat_nuis = GetNuisanceVariation(smoother, /*result_no_error, unsmoothed_result,*/ *stat_var, efficiency, bin, uniaxis, ignore_tw_relationship);
          stat_nuis->SetName(new_key->GetString().Data());
          stat_nuis->SetTitle(Form("stat_nuis_%d_of_%d", count, total_count));

          {
            auto nbins = GetTotalNbins(*stat_nuis);

            AddQuadToFrom(total_systematics, *stat_nuis);
            for (auto i = 0; i <= nbins; ++i) {
              // make errors conform to CDI convention
              stat_nuis->SetBinError(i, stat_nuis->GetBinContent(i));
            }

          }

          stat_nuis->SetDirectory(0);
          c->Add(new_key, stat_nuis);
        }

        delete stat_var;
      }

      // loop through all systematic histograms
      itt = c->MakeIterator();

      while ((k = (TObjString*)itt->Next())) {
        if ((*c)(k)->InheritsFrom("TH1")) {
          auto &h               = *static_cast<TH1*>((*c)(k));
          bool  bins_correlated = c->isBinCorrelated(k->GetString().Data());

          if ((k->GetString() == "result") ||
              k->GetString().Contains("FT_EFF_Stat Nuis") ||
              (k->GetString() == "comment") ||
              (k->GetString() == "combined") ||
              (k->GetString() == "MCreference") ||
              (k->GetString() == "MChadronisation") ||
              (k->GetString() == "ReducedSets") ||
              (k->GetString() == "extrapolation") ||
              (k->GetString() == "statistics") ||
              (k->GetString() == "systematics")) continue;

          if (bins_correlated) {
            auto &h_error         = *static_cast<TH1*>(h.Clone());
            h_error.SetDirectory(0);

            for (auto bin = 0; bin <= GetTotalNbins(h_error); ++bin)
              h_error.SetBinContent(bin, h_error.GetBinError(bin));

            smoother.UseWeightedData(false);
            Info("smoothContinuousCalibrations", "smoothing %s/%s", directory_name.Data(), key_name.Data());
            Info("smoothContinuousCalibrations", "\t\tusing bin errors - %d", smoother.UsingWeightedData());
            Info("smoothContinuousCalibrations", "\t\thistogram name - %s", k->GetString().Data());
            // for (decltype(smoother.GetDimension()) dim = 0; dim < smoother.GetDimension(); ++dim) {
            //   // Hack
            //   if (smoother.GetDimension() == 1) ++dim;
            //   const char *title;
            //   if (dim == 0)
            //     title = result.GetXaxis()->GetTitle();
            //   else if (dim == 1)
            //     title = result.GetYaxis()->GetTitle();
            //   else if (dim == 2)
            //     title = result.GetZaxis()->GetTitle();
            //   // Hack
            //   if (smoother.GetDimension() == 1) --dim;
            //   Info("smoothContinuousCalibrations", "\t\tnumber of smoothed %s bins - %lu", title, smoother.GetNbins(dim));
            // }

            GetNuisanceVariation(smoother, /*result_no_error, unsmoothed_result,*/ h, efficiency, -1, uniaxis, ignore_tw_relationship);
            GetNuisanceVariation(smoother, /*result_no_error, unsmoothed_result,*/ h_error, efficiency, -1, uniaxis, ignore_tw_relationship);

            {
              auto nbins = GetTotalNbins(h);

              AddQuadToFrom(total_systematics, h);
              for (auto i = 0; i <= nbins; ++i) h.SetBinError(i, h_error.GetBinContent(i));
            }
            delete &h_error;
          }
          else syst_nuis_names.Add(k);
        }
      }

      itt = syst_nuis_names.MakeIterator();

      while ((k = (TObjString*)itt->Next())) {
        TH1 &h        = *static_cast<TH1*>((*c)(k));
        TH1 *syst_var = static_cast<TH1*>(unsmoothed_result.Clone());
        auto nbins    = GetTotalNbins(h);

        auto &h_error = *static_cast<TH1*>(h.Clone());
        h_error.SetDirectory(0);

        for (auto bin = 0; bin <= GetTotalNbins(h_error); ++bin)
          h_error.SetBinContent(bin, h_error.GetBinError(bin));

        for (int bin = 0; bin < nbins; ++bin) {
          syst_var->SetBinContent(bin, h.GetBinContent(bin));
          syst_var->SetBinError(bin, h.GetBinError(bin));
        }

        syst_var->SetDirectory(0);

        auto count       = 0,
             total_count = 0;

        for (int bin = 0; bin <= nbins; ++bin) {
          if (syst_var->IsBinOverflow(bin) || syst_var->IsBinUnderflow(bin)) continue;
          ++total_count;
        }

        for (int bin = 0; bin <= nbins; ++bin) {
          if (syst_var->IsBinOverflow(bin) || syst_var->IsBinUnderflow(bin)) continue;
          ++count;
          TObjString *new_key = new TObjString(Form("%s Nuis %d of %d", k->GetString().Data(), count, total_count));

          smoother.UseWeightedData(false);
          Info("smoothContinuousCalibrations", "creating %s/%s", directory_name.Data(), new_key->GetString().Data());
          Info("smoothContinuousCalibrations", "\t\tusing bin errors - %d", smoother.UsingWeightedData());
          Info("smoothContinuousCalibrations", "\t\thistogram name, title - %s, %s", new_key->GetString().Data(), Form("%s nuis %d_of_%d", k->GetString().Data(), count, total_count));
          // for (decltype(smoother.GetDimension()) dim = 0; dim < smoother.GetDimension(); ++dim) {
          //   // Hack
          //   if (smoother.GetDimension() == 1) ++dim;
          //   const char *title;
          //   if (dim == 0)
          //     title = result.GetXaxis()->GetTitle();
          //   else if (dim == 1)
          //     title = result.GetYaxis()->GetTitle();
          //   else if (dim == 2)
          //     title = result.GetZaxis()->GetTitle();
          //   // Hack
          //   if (smoother.GetDimension() == 1) --dim;
          //   Info("smoothContinuousCalibrations", "\t\tnumber of smoothed %s bins - %lu", title, smoother.GetNbins(dim));
          // }

          TH1 *syst_nuis_error = GetNuisanceVariation(smoother, /*result_no_error, unsmoothed_result,*/ h_error, efficiency, bin, uniaxis, ignore_tw_relationship);
          TH1 *syst_nuis = GetNuisanceVariation(smoother, /*result_no_error, unsmoothed_result,*/ *syst_var, efficiency, bin, uniaxis, ignore_tw_relationship);
          syst_nuis->SetName(new_key->GetString().Data());
          syst_nuis->SetTitle(Form("%s nuis %d_of_%d", k->GetString().Data(), count, total_count));

          {
            auto nbins = GetTotalNbins(*syst_nuis);

            AddQuadToFrom(total_systematics, *syst_nuis);
            for (auto i = 0; i <= nbins; ++i) syst_nuis->SetBinError(i, syst_nuis_error->GetBinContent(i));
          }

          syst_nuis_error->SetDirectory(0);
          delete syst_nuis_error;
          syst_nuis->SetDirectory(0);
          c->Add(new_key, syst_nuis);
        }

        delete &h_error;
        delete syst_var;
        c->DeleteEntry(k);
      }

      // set reduction list
      if (full_name.Contains(("/" + flavorInfo["B"].fileTag + "/").c_str())) {
        TVectorD *reductionSets = new TVectorD(2);
        (*reductionSets)[0] = flavorInfo["B"].tight;
        (*reductionSets)[1] = flavorInfo["B"].medium;
        TObjString *str = new TObjString("ReducedSets");
        c->Add(str, reductionSets);
      }
      else if (full_name.Contains(("/" + flavorInfo["C"].fileTag + "/").c_str())) {
        TVectorD *reductionSets = new TVectorD(2);
        (*reductionSets)[0] = flavorInfo["C"].tight;
        (*reductionSets)[1] = flavorInfo["C"].medium;
        TObjString *str = new TObjString("ReducedSets");
        c->Add(str, reductionSets);
      }
      else if (full_name.Contains(("/" + flavorInfo["Tau"].fileTag + "/").c_str())) {
        TVectorD *reductionSets = new TVectorD(2);
        (*reductionSets)[0] = flavorInfo["Tau"].tight;
        (*reductionSets)[1] = flavorInfo["Tau"].medium;
        TObjString *str = new TObjString("ReducedSets");
        c->Add(str, reductionSets);
      }
      else if (full_name.Contains(("/" + flavorInfo["Light"].fileTag + "/").c_str())) {
        TVectorD *reductionSets = new TVectorD(2);
        (*reductionSets)[0] = flavorInfo["Light"].tight;
        (*reductionSets)[1] = flavorInfo["Light"].medium;
        TObjString *str = new TObjString("ReducedSets");
        c->Add(str, reductionSets);
      }

      delete &result_no_error;
      delete &unsmoothed_result;
    }
    else if (class_name == "Analysis::CalibrationDataHistogramContainer") {
      CalibrationDataHistogramContainer *c = static_cast<CalibrationDataHistogramContainer*>(o);

      // set reduction list
      if (full_name.Contains(("/" + flavorInfo["B"].fileTag + "/").c_str())) {
        TVectorD *reductionSets = new TVectorD(2);
        (*reductionSets)[0] = flavorInfo["B"].tight;
        (*reductionSets)[1] = flavorInfo["B"].medium;
        TObjString *str = new TObjString("ReducedSets");
        c->Add(str, reductionSets);
      }
      else if (full_name.Contains(("/" + flavorInfo["C"].fileTag + "/").c_str())) {
        TVectorD *reductionSets = new TVectorD(2);
        (*reductionSets)[0] = flavorInfo["C"].tight;
        (*reductionSets)[1] = flavorInfo["C"].medium;
        TObjString *str = new TObjString("ReducedSets");
        c->Add(str, reductionSets);
      }
      else if (full_name.Contains(("/" + flavorInfo["Tau"].fileTag + "/").c_str())) {
        TVectorD *reductionSets = new TVectorD(2);
        (*reductionSets)[0] = flavorInfo["Tau"].tight;
        (*reductionSets)[1] = flavorInfo["Tau"].medium;
        TObjString *str = new TObjString("ReducedSets");
        c->Add(str, reductionSets);
      }
      else if (full_name.Contains(("/" + flavorInfo["Light"].fileTag + "/").c_str())) {
        TVectorD *reductionSets = new TVectorD(2);
        (*reductionSets)[0] = flavorInfo["Light"].tight;
        (*reductionSets)[1] = flavorInfo["Light"].medium;
        TObjString *str = new TObjString("ReducedSets");
        c->Add(str, reductionSets);
      }
    }

    // cd to correct directory and write calibration data
    of->cd(directory_name);
    o->Write(key_name.Data(), TObject::kSingleKey);
    of->cd();
  }

  of->Write();
  of->Close();
} // smoothContinuousCalibrations
