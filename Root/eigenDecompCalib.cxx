
// File     : eigenDecompCalib.C
// Author   : Jeff Hetherly
// Purpose  : Perform simple tests of the eigenvector decomposition technique

#include <type_traits>
#include <iostream>
#include <fstream>
#include <iomanip>
#include <functional>
#include <vector>
#include <set>
#include <map>
#include <string>
#include <cmath>

#include "TROOT.h"
#include "TError.h"
#include "TFile.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TGraph2DErrors.h"
#include "TCanvas.h"
#include "TLegend.h"
#include "TLine.h"
#include "TLatex.h"
#include "TString.h"
#include "TVectorD.h"
#include "TMatrixDSym.h"
#include "TDecompSVD.h"
#include "Math/Minimizer.h"
#include "Math/Factory.h"
#include "Math/Functor.h"
// #include "Math/IFunction.h"
// #include "Minuit2/Minuit2Minimizer.h"
// #include "TMinuit.h"

#include "CalibrationDataInterface/CalibrationDataContainer.h"
#include "CalibrationDataInterface/CalibrationDataEigenVariations.h"

#include "NPandSmoothingTools/AtlasStyle.h"
#include "NPandSmoothingTools/AtlasLabels.h"
#include "NPandSmoothingTools/SmoothingUtils.h"

using namespace ROOT;
using std::vector;
using std::set;
using std::map;
using std::string;

using Analysis::CalibrationDataContainer;
using Analysis::CalibrationDataFunctionContainer;
using Analysis::CalibrationDataHistogramContainer;
using Analysis::CalibrationDataEigenVariations;


#include "TVectorT.h"
#include "TMatrixD.h"

// multidimensional minimization by the downhill simplex method (Nelder and Mead)
class Amoeba {
public:
  Amoeba (const double &ftoll = 1.0E-8) : ftol(ftoll), maxiter(500000) {}
  // minimize funciton (func here) starting from simplex Pi = P0 + delta*ei, P0 = point, delta = del
  template <class T>
  TVectorT<double> minimize (const TVectorT<double> &point, const double &del, const T &func)
  // std::valarray<double> minimize (const std::valarray<double> &point, const double &del, const T &func)
  {
    // std::valarray<double> dels(del, point.size());
    // return this->minimize(point, dels, func);
    TVectorT<double> dels(point.GetNrows());
    for (int row = 0; row < point.GetNrows(); ++row)
    dels(row) = del;
    return this->minimize(point, dels, func);
  }
  template <class T>
  TVectorT<double> minimize (const TVectorT<double> &point, const TVectorT<double> &dels, const T &func)
  // std::valarray<double> minimize (const std::valarray<double> &point, const std::valarray<double> &dels, const T &func)
  {
    int ndim = point.GetNrows();
    TMatrixD pp(ndim + 1, ndim);
    for (int row = 0; row <= ndim; ++row) {
      for (int col = 0; col < ndim; ++col) {
        pp(row, col) = point(col);
      }
      if (row != 0) pp(row, row - 1) += dels(row - 1);
    }
    return this->minimize(pp, func);
  }
  template <class T>
  TVectorT<double> minimize (const TMatrixD &pp, const T &func)
  {
    //const int NMAX = 500000;
    const double TINY = 1.0E-10;
    int ihi, ilo, inhi;
    mpts = pp.GetNrows();
    ndim = pp.GetNcols();
    TVectorT<double> psum(ndim),
    pmin(ndim),
    x(ndim);
    p.ResizeTo(pp);
    p = pp;
    y.ResizeTo(mpts);
    for (int row = 0; row < mpts; ++row) {
      for (int col = 0; col < ndim; ++col) x(col) = p(row, col);
      y(row) = func(x);
    }
    nfunc = 0;
    get_psum(p, psum);
    while (true) {
      ilo = 0;
      ihi = y(0) > y(1) ? (inhi = 1, 0) : (inhi = 0, 1);
      for (int row = 0; row < mpts; ++row) {
        if (y(row) <= y(ilo)) ilo = row;
        if (y(row) > y(ihi)) {
          inhi = ihi;
          ihi = row;
        }
        else if (y(row) > y(inhi) && row != ihi)
        inhi = row;
      }
      double rtol = 2.0*fabs(y(ihi) - y(ilo))/(fabs(y(ihi)) + fabs(y(ilo)) + TINY);
      if (rtol < ftol) {
        std::swap(y(0), y(ilo));
        for (int col = 0; col < ndim; ++col) {
          std::swap(p(0, col), p(ilo, col));
          pmin(col) = p(0, col);
        }
        fmin = y(0);
        return pmin;
      }
      if (nfunc >= maxiter) {
        //if (nfunc >= NMAX) {
        TVectorT<double> temp(pp.GetNcols());
        for (int col = 0; col < pp.GetNcols(); ++col)
        temp(col) = pp(0, col);
        Warning("minimize", "maximum iteration count(%ld) was exceeded!", maxiter);
        //  std::cerr << "WARNING: In Ameoba::minimize maximum iteration count("
        //            << maxiter << ") was exceeded!" << std::endl;
        //std::cerr << "WARNING: In minimize NMAX("
        //  << NMAX << ") was exceeded!" << std::endl;
        return temp;
      }
      nfunc += 2;
      double ytry = amotry(p, y, psum, ihi, -1.0, func);
      if (ytry <= y(ilo))
      ytry = amotry(p, y, psum, ihi, 2.0, func);
      else if (ytry >= y(inhi)) {
        double ysave = y(ihi);
        ytry = amotry(p, y, psum, ihi, 0.5, func);
        if (ytry >= ysave) {
          for (int row = 0; row < mpts; ++row) {
            if (row != ilo) {
              for (int col = 0; col < ndim; ++col)
              p(row, col) = psum(col) = 0.5*(p(row, col) + p(ilo, col));
              y(row) = func(psum);
            }
          }
          nfunc += ndim;
          get_psum(p, psum);
        }
      }
      else
      --nfunc;
    }
  }
  Amoeba& SetMaximumIteration (const long &max) {maxiter = max; return *this;}
protected:
  inline void get_psum (const TMatrixD &p, TVectorT<double> &psum)
  {
    for (int col = 0; col < ndim; ++col) {
      double sum = 0.0;
      for (int row = 0; row < mpts; ++row)
      sum += p(row, col);
      psum(col) = sum;
    }
  }
  template <class T>
  double amotry(TMatrixD &p, TVectorT<double> &y, TVectorT<double> &psum,
  const int &ihi, const double &fac, const T &func)
  {
    TVectorT<double> ptry(ndim);
    double fac1 = (1.0 - fac)/ndim;
    double fac2 = fac1 - fac;
    for (int col = 0; col < ndim; ++col)
    ptry(col) = psum(col)*fac1 - p(ihi, col)*fac2;
    double ytry = func(ptry);
    if (ytry < y(ihi)) {
      y(ihi) = ytry;
      for (int col = 0; col < ndim; ++col) {
        psum(col) += ptry(col) - p(ihi, col);
        p(ihi, col) = ptry(col);
      }
    }
    return ytry;
  }
  const double ftol;
  int mpts,
  ndim;
  long maxiter,
  nfunc;
  double fmin;
  TVectorT<double> y;
  TMatrixD p;
};


class R : public ROOT::Math::IBaseFunctionMultiDim {

  friend class RMinimizer;
  friend class RAmoeba;

  unsigned int m_ncomp = 0,
               m_nbins = 0;
  const TMatrixDSym &m_cov;

  double compute_covariance_element (const unsigned int &i, const unsigned int &j, const double *x) const
  {
    double comp_sum = 0.0;
    for (unsigned int k = 0; k < m_ncomp; ++k)
      comp_sum += x[i + k*m_nbins]*x[j + k*m_nbins];
    return comp_sum;
  }

  virtual double DoEval(const double *x) const
  {
    double result = 0.0;
    // double max_rel = 0.0;
    for (unsigned int i = 0; i < m_nbins; ++i) {
      for (unsigned int j = i; j < m_nbins; ++j) {
        double rel = TMath::Abs((compute_covariance_element(i, j, x) - m_cov(i, j))/m_cov(i, j));
        result += (i != j ? 2.0 : 1.0)*rel;
      }
    }
    return result;
  }

  // virtual double DoDerivative(const double* x, unsigned int icoord) const
  // {
  //   double sum = 0.0;
  //   unsigned int i = icoord % m_nbins,
  //                k_prime = icoord / m_nbins;
  //   for (unsigned int j = 0; j < m_nbins; ++j) {
  //     sum += (compute_covariance_element(i, j, x) < 0 ? -1.0 : 1.0)*x[j + k_prime*m_nbins];
  //   }
  //   return x[icoord] < 0 ? -1 : 1;
  // }

public:

  R (const TMatrixDSym &cov, const unsigned int &ncomp) : m_ncomp(ncomp), m_nbins(cov.GetNrows()), m_cov(cov) {}
  virtual unsigned int NDim() const {return m_ncomp*m_cov.GetNrows();}
  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const {return new R(m_cov, m_ncomp);}

  TMatrixDSym ComputeCovarianceMatrix (const double *x) const
  {
    TMatrixDSym result(m_nbins);
    for (unsigned int i = 0; i < m_nbins; ++i) {
      for (unsigned int j = i; j < m_nbins; ++j) {
        result(i, j) = compute_covariance_element(i, j, x);
        result(j, i) = result(i, j);
      }
    }
    return result;
  }

  TMatrixDSym ComputeRelativeDifferenceCovarianceMatrix (const double *x) const
  {
    TMatrixDSym result(m_nbins);
    for (unsigned int i = 0; i < m_nbins; ++i) {
      for (unsigned int j = i; j < m_nbins; ++j) {
        result(i, j) = (compute_covariance_element(i, j, x) - m_cov(i, j))/m_cov(i, j);
        result(j, i) = result(i, j);
      }
    }
    return result;
  }
};

class RMinimizer : public ROOT::Math::IBaseFunctionMultiDim  {

  mutable bool m_minimizerSuccess = true;
  R m_func;
  // R func(covMat, 5);
  // create minimizer giving a name and a name (optionally) for the specific
  // algorithm
  // possible choices are:
  //     minName                  algoName
  // Minuit /Minuit2             Migrad, Simplex,Combined,Scan  (default is Migrad)
  //  Minuit2                     Fumili2
  //  Fumili
  //  GSLMultiMin                ConjugateFR, ConjugatePR, BFGS,
  //                              BFGS2, SteepestDescent
  //  GSLMultiFit
  //   GSLSimAn
  //   Genetic
  ROOT::Math::Minimizer *m_minimizer = nullptr;

  virtual double DoEval(const double *x) const
  {
    for (unsigned int i = 0; i < m_func.NDim(); ++i)
      m_minimizer->SetVariable(i, std::to_string(i), 0, x[0]);
    m_minimizerSuccess = m_minimizer->Minimize();
    auto rel_cov = m_func.ComputeRelativeDifferenceCovarianceMatrix(m_minimizer->X());
    double max_rel_diff = 0.0;
    for (int i = 0; i < rel_cov.GetNrows(); ++i) {
      for (int j = i; j < rel_cov.GetNcols(); ++j) {
        auto abs_rel = TMath::Abs(rel_cov(i, j));
        max_rel_diff = max_rel_diff < abs_rel ? abs_rel : max_rel_diff;
      }
    }
    return max_rel_diff;
  }

public:
  RMinimizer (const TMatrixDSym &cov, const unsigned int &ncomp, const char *minName = "Minuit2", const char *algoName = "Simplex") :
    m_func(cov, ncomp), m_minimizer(ROOT::Math::Factory::CreateMinimizer(minName, algoName))
  {
    m_minimizer->SetPrecision(0.001);
    m_minimizer->SetMaxIterations(100000000); // for GSL
    m_minimizer->SetMaxFunctionCalls(100000000); // for Minuit/Minuit2
    m_minimizer->SetTolerance(0.001);
    m_minimizer->SetPrintLevel(0);

    m_minimizer->SetFunction(m_func);
  }

  ~RMinimizer ()
  {
    if (m_minimizer) delete m_minimizer;
  }

  virtual unsigned int NDim() const {return 1;}
  virtual ROOT::Math::IBaseFunctionMultiDim* Clone() const
  {
    return new RMinimizer(m_func.m_cov, m_func.m_ncomp, m_minimizer->Options().MinimizerType().c_str(), m_minimizer->Options().MinimizerAlgorithm().c_str());
  }

  R& GetFunction () {return m_func;}
  const R& GetFunction () const {return m_func;}
  ROOT::Math::Minimizer& GetMinimizer () {return *m_minimizer;}
  const ROOT::Math::Minimizer& GetMinimizer () const {return *m_minimizer;}

  bool Successful () const {return m_minimizerSuccess;}
  bool Scan () const
  {
    return false;
  }
};


class RAmoeba {
  Amoeba m_minimizer;
  R m_func;
  TVectorT<double> m_initials,
                   m_result;

public:
  RAmoeba (const TMatrixDSym &cov, const unsigned int &ncomp) :
    m_minimizer(1.0e-6), m_func(cov, ncomp)
  {
    m_result.ResizeTo(m_func.NDim());
    m_initials.ResizeTo(m_func.NDim());
    for (int i = 0; i < m_initials.GetNoElements(); ++i) m_initials[i] = 0.0;
    m_minimizer.SetMaximumIteration(1000000);
  }

  R& GetFunction () {return m_func;}
  const R& GetFunction () const {return m_func;}
  Amoeba& GetMinimizer () {return m_minimizer;}
  const Amoeba& GetMinimizer () const {return m_minimizer;}

  double* operator () (const double &x)
  {
    m_result = m_minimizer.minimize(m_initials, x, *this);
    return &m_result[0];
  }

  double operator () (const TVectorT<double> &x) const
  {
    return m_func(&x[0]);
  }
};


// Helper functions
namespace {

struct AnalysisData
{
  TString &info_tag,
          &uinfo_tag,
          &InfoTag,
          &uInfoTag;
  TH1 &result,
      &uresult;
  vector<TH1*> &stat_set,
               &ustat_set,
               &syst_set,
               &usyst_set,
               &eigen_set,
               &ueigen_set;
  TMatrixDSym &original_cov,
              &uoriginal_cov,
              &eigen_cov,
              &ueigen_cov;
};

template<class T1, class T2>
T2 round_to (T1 num, T2 precision)
{
  float correction = (num >= 0) ? 0.5 : -0.5;

  return int(num / precision + correction) * precision;
}

void DoEigenAnalysis(TFile *of, TString pathName, TString tag, TCanvas *c1,
    const vector<unsigned long> &colors, const AnalysisData &data)
{
  bool setLog = (TString(data.result.GetXaxis()->GetTitle()) == "pt");

  // ////////////////////////////////////////////////////////
  // Statistical Nuisances (plotting)
  // ////////////////////////////////////////////////////////
  auto *total_stat_error    = static_cast<TH1*>(data.result.Clone(Form("total_stat_error%s", data.info_tag.Data()))),
       *utotal_stat_error   = static_cast<TH1*>(data.uresult.Clone(Form("total_ustat_error%s", data.info_tag.Data())));
  total_stat_error->Reset();
  utotal_stat_error->Reset();
  {
    Info("eigenDecompCalib", "plotting statistical nuisances: %s", data.InfoTag.Data());
    unsigned color = 1;
    utotal_stat_error->SetName(Form("ustat_error%s", data.info_tag.Data()));
    TLegend *leg = new TLegend(0.67, 0.75, 0.91, 0.90);

    leg->SetBorderSize(0);

    for (auto histo : data.ustat_set) {
      for (auto bin = 1; bin <= histo->GetNbinsX(); ++bin) {
        auto current_error    = utotal_stat_error->GetBinContent(bin),
             additional_error = histo->GetBinContent(bin);
        utotal_stat_error->SetBinContent(bin, sqrt(current_error * current_error + additional_error * additional_error));
      }
    }

    double low_bound  = 0.0,
           high_bound = 1.5 * (utotal_stat_error->GetBinContent(utotal_stat_error->GetMaximumBin()) - low_bound) + low_bound;

    // utotal_stat_error->SetTitle("Statistical \"Nuisance\" Variations");
    utotal_stat_error->SetTitle("");
    utotal_stat_error->SetLineColor(kBlack);
    utotal_stat_error->SetStats(kFALSE);
    utotal_stat_error->GetXaxis()->SetTitleOffset(1.4);
    utotal_stat_error->GetYaxis()->SetTitleOffset(1.5);
    utotal_stat_error->GetYaxis()->SetTitle("Statistical SF Error");
    if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt")
      utotal_stat_error->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta")
      utotal_stat_error->GetXaxis()->SetTitle("#eta");
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta")
      utotal_stat_error->GetXaxis()->SetTitle("#left|#eta#right|");
    utotal_stat_error->SetAxisRange(low_bound, high_bound, "Y");
    utotal_stat_error->GetXaxis()->SetLabelSize(utotal_stat_error->GetYaxis()->GetLabelSize());
    utotal_stat_error->GetXaxis()->SetTitleSize(utotal_stat_error->GetYaxis()->GetTitleSize());
    utotal_stat_error->Draw("hist");
    leg->AddEntry(utotal_stat_error, "Total Stat. Error", "l");

    color = 0;
    map<string, size_t> color_map;

    for (auto histo : data.stat_set) {
      string key("Nuis"), title(histo->GetName());
      if (title.find(key) != string::npos)
        title.erase(title.rfind(key), title.size() - title.rfind(key) + 1);
      if (title.find("FT_EFF_") != string::npos)
        title.erase(title.find("FT_EFF_"), string("FT_EFF_").size());
      auto add_legend = false;
      if (color_map.count(title) == 0) {
        color_map[title] = colors[color % colors.size()];
        ++color;
        add_legend = true;
      }

      for (int bin = 1; bin <= total_stat_error->GetNbinsX(); ++bin)
        total_stat_error->SetBinContent(bin,
            total_stat_error->GetBinContent(bin) +
            histo->GetBinContent(bin) * histo->GetBinContent(bin));
      for (int bin = 1; bin <= histo->GetNbinsX(); ++bin)
        histo->SetBinContent(bin, TMath::Abs(histo->GetBinContent(bin)));

      histo->SetLineColor(color_map[title]);
      histo->Draw("hist same");

      if (add_legend) leg->AddEntry(histo, title.c_str(), "l");
    }

    for (int bin = 1; bin <= total_stat_error->GetNbinsX(); ++bin)
      total_stat_error->SetBinContent(bin,
          sqrt(total_stat_error->GetBinContent(bin)));

    total_stat_error->SetLineColor(kBlack);
    total_stat_error->Draw("hist same");
    leg->Draw();

    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    myText(0.175, 0.9, 0.035, kBlue, pathName.Data());
  }

  c1->SetLogx(setLog);
  c1->SetName(Form("stat_nuis%s", data.info_tag.Data()));
  of->WriteTObject(c1);
  c1->Clear();

  // ////////////////////////////////////////////////////////
  // Eigen-Variations (plotting)
  // ////////////////////////////////////////////////////////
  TH1D *total_error = static_cast<TH1D*>(data.result.Clone(Form("total_error%s", data.info_tag.Data()))),
       *utotal_error     = static_cast<TH1D*>(data.uresult.Clone(Form("total_uerror%s", data.info_tag.Data())));
  total_error->Reset();
  utotal_error->Reset();
  {
    vector<TH1D*> cumulative_eigens, single_eigens,
      ucumulative_eigens, usingle_eigens;
    unsigned int color                    = 1,
                 number_of_eigens_to_plot = 3;
    double low_bound, high_bound;

    auto max_neigen_vars = std::max(data.eigen_set.size(), data.ueigen_set.size());

    for (unsigned int i = 0; i < max_neigen_vars; ++i) {
      auto   drawn = false;
      double y_max = 0.0,
             y_min = 0.0;
      {
        if (i < data.ueigen_set.size()) {
          y_max = data.ueigen_set[i]->GetBinContent(data.ueigen_set[i]->GetMaximumBin());
          y_min = data.ueigen_set[i]->GetBinContent(data.ueigen_set[i]->GetMinimumBin());
        }

        if (i < data.eigen_set.size()) {
          y_max = std::max(y_max, data.eigen_set[i]->GetBinContent(data.eigen_set[i]->GetMaximumBin()));
          y_min = std::min(y_min, data.eigen_set[i]->GetBinContent(data.eigen_set[i]->GetMinimumBin()));
        }
        y_min = y_min < 0.0 ? 1.2 * y_min : 0.8 * y_min;
        y_max = 1.2 * (y_max - y_min) + y_min;
      }

      if (i < data.ueigen_set.size()) {
        TH1 *up0 = data.ueigen_set[i];

        for (int bin = 1; bin <= utotal_error->GetNbinsX(); ++bin)
          utotal_error->SetBinContent(bin,
              utotal_error->GetBinContent(bin) +
              up0->GetBinContent(bin) * up0->GetBinContent(bin));

        if (i < number_of_eigens_to_plot) {
          TH1D *total_temp  = static_cast<TH1D*>(utotal_error->Clone());
          TH1D *single_temp = static_cast<TH1D*>(up0->Clone("|single|"));

          for (int bin = 1; bin <= single_temp->GetNbinsX(); ++bin) {
            total_temp->SetBinContent(bin, sqrt(total_temp->GetBinContent(bin)));
            single_temp->SetBinContent(bin, fabs(single_temp->GetBinContent(bin)));
          }
          ucumulative_eigens.push_back(total_temp);
          usingle_eigens.push_back(single_temp);
        }

        up0->SetLineColor(kBlack);

        // up0->SetTitle(TString::Format("Eigen-Variation #%d", i + 1).Data());
        up0->SetTitle("");
        up0->SetAxisRange(y_min, y_max, "Y");
        up0->SetStats(kFALSE);
        up0->GetXaxis()->SetLabelSize(up0->GetYaxis()->GetLabelSize());
        up0->GetXaxis()->SetTitleSize(up0->GetYaxis()->GetTitleSize());

        up0->GetYaxis()->SetTitle(Form("Eigen-Variation %d", i));
        if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt")
          up0->GetXaxis()->SetTitle("p_{T} [GeV]");
        else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta")
          up0->GetXaxis()->SetTitle("#eta");
        else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta")
          up0->GetXaxis()->SetTitle("#left|#eta#right|");
        up0->Draw("hist");

        myText(0.35, 0.9, 0.04, kBlue, pathName.Data());
        drawn = true;
      }

      if (i < data.eigen_set.size()) {
        TH1 *up0 = data.eigen_set[i];

        for (int bin = 1; bin <= total_error->GetNbinsX(); ++bin) {
          total_error->SetBinContent(bin,
              total_error->GetBinContent(bin) +
              up0->GetBinContent(bin) * up0->GetBinContent(bin));
        }

        if (i < number_of_eigens_to_plot) {
          TH1D *total_temp  = static_cast<TH1D*>(total_error->Clone());
          TH1D *single_temp = static_cast<TH1D*>(up0->Clone("|single|"));

          for (int bin = 1; bin <= single_temp->GetNbinsX(); ++bin) {
            total_temp->SetBinContent(bin, sqrt(total_temp->GetBinContent(bin)));
            single_temp->SetBinContent(bin, fabs(single_temp->GetBinContent(bin)));
          }
          cumulative_eigens.push_back(total_temp);
          single_eigens.push_back(single_temp);
        }

        up0->SetMarkerColor(colors[i % colors.size()]);
        up0->SetLineColor(colors[i % colors.size()]);

        // up0->SetTitle(TString::Format("Eigen-Variation #%d", i + 1).Data());
        up0->SetTitle("");
        up0->SetStats(kFALSE);
        up0->GetXaxis()->SetTitleOffset(0.9);
        up0->GetYaxis()->SetTitleOffset(1.2);
        up0->GetYaxis()->SetTitle(Form("Eigen-Variation %d", i));
        if (TString(data.result.GetXaxis()->GetTitle()) == "pt")
          up0->GetXaxis()->SetTitle("p_{T} [GeV]");
        else if (TString(data.result.GetXaxis()->GetTitle()) == "eta")
          up0->GetXaxis()->SetTitle("#eta");
        else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta")
          up0->GetXaxis()->SetTitle("#left|#eta#right|");

        if (drawn) up0->Draw("hist same");
        else up0->Draw("hist");
      }

      myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
      c1->SetName(Form("eigen-variation_%d%s", i + 1, data.info_tag.Data()));
      of->WriteTObject(c1);
      c1->Clear();
    }

    for (int bin = 1; bin <= total_error->GetNbinsX(); ++bin)
      total_error->SetBinContent(bin,
          sqrt(total_error->GetBinContent(bin)));

    for (int bin = 1; bin <= utotal_error->GetNbinsX(); ++bin)
      utotal_error->SetBinContent(bin,
          sqrt(utotal_error->GetBinContent(bin)));

    // ////////////////////////////////////////////////////////
    // Total Error (smoothed versus unsmoothed)
    // ////////////////////////////////////////////////////////
    // utotal_error->SetTitle("Total Error");
    utotal_error->SetTitle("");
    utotal_error->SetLineColor(kBlack);
    utotal_error->SetStats(kFALSE);
    utotal_error->GetXaxis()->SetLabelSize(utotal_error->GetYaxis()->GetLabelSize());
    utotal_error->GetXaxis()->SetTitleSize(utotal_error->GetYaxis()->GetTitleSize());
    utotal_error->GetXaxis()->SetTitleOffset(0.9);
    utotal_error->GetYaxis()->SetTitleOffset(1.2);
    utotal_error->GetYaxis()->SetTitle("Total SF Error");
    if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt")
      utotal_error->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta")
      utotal_error->GetXaxis()->SetTitle("#eta");
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta")
      utotal_error->GetXaxis()->SetTitle("#left|#eta#right|");
    utotal_error->Draw("hist");
    total_error->SetLineColor(kBlack);
    total_error->Draw("hist same");

    myText(0.25, 0.85, 0.04, kBlue, pathName.Data());

    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    c1->SetName(Form("total_error%s", data.info_tag.Data()));
    of->WriteTObject(c1);
    c1->Clear();


    // ////////////////////////////////////////////////////////
    // Total Error (cumulative)
    // ////////////////////////////////////////////////////////
    {
      high_bound = total_error->GetBinContent(total_error->GetMaximumBin());
      high_bound = std::max(high_bound, utotal_error->GetBinContent(utotal_error->GetMaximumBin()));
      low_bound  = high_bound;

      for (auto hist : cumulative_eigens) {
        low_bound = std::min(low_bound, hist->GetBinContent(hist->GetMinimumBin()));
      }

      for (auto hist : ucumulative_eigens) {
        low_bound = std::min(low_bound, hist->GetBinContent(hist->GetMinimumBin()));
      }
      low_bound  = low_bound < 0.0 ? 1.1 * low_bound : 0.8 * low_bound;
      high_bound = 1.2 * (high_bound - low_bound) + low_bound;
    }

    TLegend *leg = new TLegend(0.2, 0.65, 0.65, 0.95);

    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    // total_error->SetTitle(TString::Format("Total Error Band: Cumulative Leading %d EV",
    //                                    number_of_eigens_to_plot).Data());
    total_error->SetTitle("");
    total_error->SetLineWidth(3);
    total_error->SetLineColor(kBlack);
    total_error->SetStats(kFALSE);
    total_error->SetAxisRange(low_bound, high_bound, "Y");
    total_error->GetXaxis()->SetTitleOffset(0.9);
    total_error->GetYaxis()->SetTitleOffset(1.2);
    total_error->GetYaxis()->SetTitle("SF Error");
    if (TString(data.result.GetXaxis()->GetTitle()) == "pt")
      total_error->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (TString(data.result.GetXaxis()->GetTitle()) == "eta")
      total_error->GetXaxis()->SetTitle("#eta");
    else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta")
      total_error->GetXaxis()->SetTitle("#left|#eta#right|");
    total_error->GetXaxis()->SetLabelSize(total_error->GetYaxis()->GetLabelSize());
    total_error->GetXaxis()->SetTitleSize(total_error->GetYaxis()->GetTitleSize());
    total_error->Draw("hist");
    leg->AddEntry(total_error, "Total Error", "l");
    color = 0;
    string leg_text("EV #");

    for (auto hist : cumulative_eigens) {
      hist->SetLineColor(colors[color % colors.size()]);
      hist->Draw("hist same");
      ++color;
      leg_text += std::to_string(color);
      leg->AddEntry(hist, leg_text.c_str(), "l");
      leg_text += " + EV #";
    }

    utotal_error->SetLineWidth(3);
    utotal_error->SetLineColor(kBlack);
    utotal_error->Draw("hist same");
    color = 0;

    for (auto hist : ucumulative_eigens) {
      hist->SetLineColor(colors[color % colors.size()]);
      hist->Draw("hist same");
      ++color;
    }

    // lt->SetTextColor(kBlue);
    // lt->DrawLatexNDC( 0.15, 0.85, pathName.Data() );
    myText(0.25, 0.6, 0.04, kBlue, pathName.Data());
    ATLASLabel(0.5, 0.9, "Internal");

    leg->Draw();
    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    c1->SetName(Form("cumulative_total_error%s", data.info_tag.Data()));
    of->WriteTObject(c1);
    c1->Clear();
    leg->Clear();


    // ////////////////////////////////////////////////////////
    // Total Error (single)
    // ////////////////////////////////////////////////////////
    // total_error->SetTitle(TString::Format("Total Error Band: Individual Leading %d EV (Absolute Value)",
    //                                    number_of_eigens_to_plot).Data());
    total_error->SetTitle("");
    leg = new TLegend(0.2, 0.7, 0.4, 0.95);

    leg->SetFillStyle(0);
    leg->SetBorderSize(0);

    total_error->SetLineWidth(3);
    total_error->SetLineColor(kBlack);
    total_error->SetStats(kFALSE);
    total_error->SetAxisRange(0, high_bound, "Y");
    total_error->GetXaxis()->SetTitleOffset(0.9);
    total_error->GetYaxis()->SetTitleOffset(1.2);
    total_error->GetYaxis()->SetTitle("SF Error");
    if (TString(data.result.GetXaxis()->GetTitle()) == "pt")
      total_error->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (TString(data.result.GetXaxis()->GetTitle()) == "eta")
      total_error->GetXaxis()->SetTitle("#eta");
    else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta")
      total_error->GetXaxis()->SetTitle("#left|#eta#right|");
    total_error->GetXaxis()->SetLabelSize(total_error->GetYaxis()->GetLabelSize());
    total_error->GetXaxis()->SetTitleSize(total_error->GetYaxis()->GetTitleSize());
    total_error->Draw("hist");
    leg->AddEntry(total_error, "Total Error", "l");
    color = 0;

    for (auto hist : single_eigens) {
      // hist->SetLineStyle(kDashDotted);
      // hist->SetLineStyle(kDashed);
      hist->SetLineColor(colors[color % colors.size()]);
      hist->Draw("hist same");
      ++color;
      leg->AddEntry(hist, ("EV #" + std::to_string(color)).c_str(), "l");
    }

    utotal_error->SetLineWidth(3);
    utotal_error->SetLineColor(kBlack);
    utotal_error->Draw("hist same");
    color = 0;

    for (auto hist : usingle_eigens) {
      hist->SetLineColor(colors[color % colors.size()]);
      hist->Draw("hist same");
      ++color;
    }

    myText(0.35, 0.8, 0.04, kBlue, pathName.Data());
    ATLASLabel(0.5, 0.85, "Internal");

    leg->Draw();

    c1->SetName(Form("single_total_error%s", data.info_tag.Data()));
    of->WriteTObject(c1);
    c1->Clear();
  }


  {
    // ////////////////////////////////////////////////////////
    // Error Bands
    // ////////////////////////////////////////////////////////

    Info("eigenDecompCalib", "plotting error bands: %s", data.InfoTag.Data());
    auto *stat_var_up      = static_cast<TH1*>(data.result.Clone()),
         *stat_var_down    = static_cast<TH1*>(data.result.Clone()),
         *total_var_up     = static_cast<TH1*>(data.result.Clone()),
         *total_var_down   = static_cast<TH1*>(data.result.Clone()),
         *ustat_result     = static_cast<TH1*>(data.uresult.Clone()),
         *utotal_result    = static_cast<TH1*>(data.uresult.Clone());
    TLegend *leg      = new TLegend(0.65, 0.75, 0.9, 0.90);

    stat_var_up->Reset();
    stat_var_down->Reset();
    total_var_up->Reset();
    total_var_down->Reset();
    leg->SetBorderSize(0);

    // leg->SetFillStyle(0);

    for (int bin = 1; bin <= stat_var_up->GetNbinsX(); ++bin) {
      stat_var_up->SetBinContent(bin,
          data.result.GetBinContent(bin) +
          total_stat_error->GetBinContent(bin));
      stat_var_down->SetBinContent(bin,
          data.result.GetBinContent(bin) -
          total_stat_error->GetBinContent(bin));
      total_var_up->SetBinContent(bin,
          data.result.GetBinContent(bin) +
          total_error->GetBinContent(bin));
      total_var_down->SetBinContent(bin,
          data.result.GetBinContent(bin) -
          total_error->GetBinContent(bin));
    }

    for (auto bin = 1; bin <= ustat_result->GetNbinsX(); ++bin) ustat_result->SetBinError(bin, utotal_stat_error->GetBinContent(bin));

    for (int bin = 1; bin <= data.uresult.GetNbinsX(); ++bin) utotal_result->SetBinError(bin, utotal_error->GetBinContent(bin));

    unsigned int total_color = 8,
                 stat_color  = 9;

    // utotal_result->SetTitle("Error Bands");
    utotal_result->SetTitle("");
    utotal_result->SetLineColor(total_color);
    utotal_result->SetMarkerColor(kBlack);
    utotal_result->SetMarkerStyle(20);
    utotal_result->SetStats(kFALSE);
    utotal_result->GetXaxis()->SetTitleOffset(1.4);
    utotal_result->GetYaxis()->SetTitleOffset(1.5);
    utotal_result->GetYaxis()->SetTitle("scale factor");
    if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt")
      utotal_result->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta")
      utotal_result->GetXaxis()->SetTitle("#eta");
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta")
      utotal_result->GetXaxis()->SetTitle("#left|#eta#right|");
    utotal_result->GetXaxis()->SetLabelSize(utotal_result->GetYaxis()->GetLabelSize());
    utotal_result->GetXaxis()->SetTitleSize(utotal_result->GetYaxis()->GetTitleSize());
    utotal_result->Draw("E1 X0");
    ustat_result->SetLineColor(stat_color);
    ustat_result->SetMarkerColor(kBlack);
    ustat_result->SetMarkerStyle(20);

    // ustat_result->SetTitle("");
    // ustat_result->GetYaxis()->SetTitle("scale factor");
    // ustat_result->GetXaxis()->SetLabelSize(ustat_result->GetYaxis()->GetLabelSize());
    // ustat_result->GetXaxis()->SetTitleSize(ustat_result->GetYaxis()->GetTitleSize());
    // ustat_result->Draw("E1 X0");
    ustat_result->Draw("E1 X0 same");
    data.result.SetLineColor(kBlack);
    data.result.Draw("hist same");
    stat_var_up->SetLineColor(stat_color);
    stat_var_up->Draw("hist same");
    stat_var_down->SetLineColor(stat_color);
    stat_var_down->Draw("hist same");
    total_var_up->SetLineColor(total_color);
    total_var_up->Draw("hist same");
    total_var_down->SetLineColor(total_color);
    total_var_down->Draw("hist same");
    leg->AddEntry(&data.result, "Central", "l");
    leg->AddEntry(stat_var_up, "Stat. Error", "l");
    leg->AddEntry(total_var_up, "Stat. + Syst. Error", "l");
    leg->Draw();

    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    myText(0.2, 0.9, 0.035, kBlue, pathName.Data());

    c1->SetLogx(setLog);
    c1->SetName(Form("central%s", data.info_tag.Data()));
    of->WriteTObject(c1);
    c1->Clear();
  }

  {
    // ////////////////////////////////////////////////////////
    // Covariance and Correlations Matrices (construction)
    // ////////////////////////////////////////////////////////
    // original versus eigen-decomp (Var) covariance matrices (smoothed)
    TMatrixDSym covMat    = data.original_cov;
    TMatrixDSym covMatVar = data.eigen_cov;
    TMatrixDSym corrMat(covMat.GetNrows());
    TMatrixDSym corrMatVar(covMatVar.GetNrows());
    TMatrixDSym covarianceDiff(covMatVar.GetNrows());
    TMatrixDSym correlationDiff(corrMatVar.GetNrows());
    // original versus eigen-decomp (Var) covariance matrices (unsmoothed)
    TMatrixDSym ucovMat    = data.uoriginal_cov;
    TMatrixDSym ucovMatVar = data.ueigen_cov;
    TMatrixDSym ucorrMat(ucovMat.GetNrows());
    TMatrixDSym ucorrMatVar(ucovMatVar.GetNrows());
    TMatrixDSym ucovarianceDiff(ucovMatVar.GetNrows());
    TMatrixDSym ucorrelationDiff(ucorrMatVar.GetNrows());

    // TDecompSVD svd(covMat, 1.0e-2);
    // auto sigs = svd.GetSig();
    // for (int row = 0; row < sigs.GetNrows(); ++row)
    //   cout << "(" << row << "): " << sigs(row) << endl;

    corrMat.Zero();
    corrMatVar.Zero();
    correlationDiff.Zero();
    covarianceDiff.Zero();
    ucorrMat.Zero();
    ucorrMatVar.Zero();
    ucorrelationDiff.Zero();
    ucovarianceDiff.Zero();

    for (int u = 0; u < covarianceDiff.GetNrows(); ++u) {
      for (int v = 0; v < covarianceDiff.GetNcols(); ++v) {
        if (covMat(u, v) != 0.0) covarianceDiff(u, v) = covMatVar(u, v) / covMat(u, v) - 1.0;
      }
    }

    for (int u = 0; u < ucovarianceDiff.GetNrows(); ++u) {
      for (int v = 0; v < ucovarianceDiff.GetNcols(); ++v) {
        if (ucovMat(u, v) != 0.0) ucovarianceDiff(u, v) = ucovMatVar(u, v) / ucovMat(u, v) - 1.0;
      }
    }

    for (int u = 0; u < correlationDiff.GetNrows(); ++u) {
      for (int v = 0; v < correlationDiff.GetNcols(); ++v) {
        corrMat(u, v)    = covMat(u, v) / sqrt(covMat(u, u) * covMat(v, v));
        corrMatVar(u, v) = covMatVar(u, v) / sqrt(covMatVar(u, u) * covMatVar(v, v));

        if (corrMat(u, v) != 0.0) correlationDiff(u, v) = corrMatVar(u, v) / corrMat(u, v) - 1.0;
      }
    }

    for (int u = 0; u < ucorrelationDiff.GetNrows(); ++u) {
      for (int v = 0; v < ucorrelationDiff.GetNcols(); ++v) {
        ucorrMat(u, v)    = ucovMat(u, v) / sqrt(ucovMat(u, u) * ucovMat(v, v));
        ucorrMatVar(u, v) = ucovMatVar(u, v) / sqrt(ucovMatVar(u, u) * ucovMatVar(v, v));

        if (ucorrMat(u, v) != 0.0) ucorrelationDiff(u, v) = ucorrMatVar(u, v) / ucorrMat(u, v) - 1.0;
      }
    }

    // ////////////////////////////////////////////////////////
    // Covariance and Correlations Matrices (plotting)
    // ////////////////////////////////////////////////////////
    // c1 = new TCanvas( "cor_eigen_histo",
    //                      "Eigen-Variation Correlation Matrix Relative Difference",
    //                      200, 10, 800, 800 );
    auto old_left_margin   = c1->GetLeftMargin(),
         old_right_margin  = c1->GetRightMargin(),
         old_bottom_margin = c1->GetBottomMargin();

    c1->SetRightMargin(old_left_margin);

    // c1->SetLeftMargin(0.25*old_left_margin);
    // c1->SetBottomMargin(0.25*old_bottom_margin);

    TList lines;

    TH2D *cor_histo = new TH2D(Form("cor_eigen_histo%s", data.info_tag.Data()), "",
        correlationDiff.GetNrows(), 1, correlationDiff.GetNrows(),
        correlationDiff.GetNrows(), 1, correlationDiff.GetNrows());


    if (data.result.GetNbinsX() == correlationDiff.GetNrows())
      cor_histo->SetBins(data.result.GetNbinsX(),
          data.result.GetXaxis()->GetXbins()->GetArray(),
          data.result.GetNbinsX(),
          data.result.GetXaxis()->GetXbins()->GetArray());

    for (int i = 1; i <= cor_histo->GetNbinsX(); ++i) {
      for (int j = i; j <= cor_histo->GetNbinsY(); ++j) {
        cor_histo->SetBinContent(i, j, 100.0 * correlationDiff(i - 1, j - 1));
        cor_histo->SetBinContent(j, i, cor_histo->GetBinContent(i, j));
      }
    }

    cor_histo->GetXaxis()->SetTitleOffset(0.9);
    cor_histo->GetYaxis()->SetTitleOffset(0.9);

    if (data.result.GetNbinsX() == correlationDiff.GetNrows()) {
      if (TString(data.result.GetXaxis()->GetTitle()) == "pt") {
        cor_histo->SetXTitle("p_{T} [GeV]");
        cor_histo->SetYTitle("p_{T} [GeV]");
      }
      else if (TString(data.result.GetXaxis()->GetTitle()) == "eta") {
        cor_histo->SetXTitle("#eta");
        cor_histo->SetYTitle("#eta");
      }
      else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta") {
        cor_histo->SetXTitle("#left|#eta#right|");
        cor_histo->SetYTitle("#left|#eta#right|");
      }
    }
    cor_histo->GetXaxis()->SetRange(1, correlationDiff.GetNrows());
    cor_histo->GetYaxis()->SetRange(1, correlationDiff.GetNrows());
    cor_histo->SetStats(kFALSE);
    cor_histo->Draw("COLZ");

    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    myText(0.2, 0.965, 0.035, kBlack, "Eigen-Variation Correlation Matrix % Relative Difference");

    c1->SetLogx(kFALSE);
    c1->SetLogy(kFALSE);

    if (data.result.GetNbinsX() == correlationDiff.GetNrows()) {
      c1->SetLogx(setLog);
      c1->SetLogy(setLog);
    }
    c1->SetName(Form("correlation_difference%s", data.info_tag.Data()));

    if (data.result.GetNbinsX() == correlationDiff.GetNrows()) {
      for (int ii = 0; ii < data.uresult.GetNbinsX() - 1; ++ii) {
        double x    = data.uresult.GetBinCenter(ii + 1) + data.uresult.GetBinWidth(ii + 1) / 2.0,
               xmin = data.uresult.GetXaxis()->GetXmin(),
               xmax = data.uresult.GetXaxis()->GetXmax();
        TLine *vl   = new TLine(x, xmin, x, xmax);
        vl->Draw();
        lines.Add(vl);
        TLine *hl = new TLine(xmin, x, xmax, x);
        hl->Draw();
        lines.Add(hl);
      }
    }

    of->WriteTObject(c1);
    c1->Clear();

    TH2D *ucor_histo = new TH2D(Form("ucor_eigen_histo%s", data.info_tag.Data()), "",
        ucorrelationDiff.GetNrows(), 1, ucorrelationDiff.GetNrows(),
        ucorrelationDiff.GetNrows(), 1, ucorrelationDiff.GetNrows());

    ucor_histo->SetBins(data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray(),
        data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray());

    for (int i = 1; i <= ucor_histo->GetNbinsX(); ++i) {
      for (int j = i; j <= ucor_histo->GetNbinsY(); ++j) {
        ucor_histo->SetBinContent(i, j, 100.0 * ucorrelationDiff(i - 1, j - 1));
        ucor_histo->SetBinContent(j, i, ucor_histo->GetBinContent(i, j));
      }
    }

    ucor_histo->GetXaxis()->SetTitleOffset(0.9);
    ucor_histo->GetYaxis()->SetTitleOffset(0.9);
    if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt") {
      ucor_histo->SetXTitle("p_{T} [GeV]");
      ucor_histo->SetYTitle("p_{T} [GeV]");
    }
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta") {
      ucor_histo->SetXTitle("#eta");
      ucor_histo->SetYTitle("#eta");
    }
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta") {
      ucor_histo->SetXTitle("#left|#eta#right|");
      ucor_histo->SetYTitle("#left|#eta#right|");
    }
    ucor_histo->GetXaxis()->SetRange(1, ucorrelationDiff.GetNrows());
    ucor_histo->GetYaxis()->SetRange(1, ucorrelationDiff.GetNrows());
    ucor_histo->SetStats(kFALSE);
    ucor_histo->Draw("COLZ");

    myText(0.2, 0.05, 0.05, kBlue, data.uInfoTag.Data());
    myText(0.2, 0.965, 0.035, kBlack, "Eigen-Variation Correlation Matrix % Relative Difference (unsmoothed)");

    c1->SetLogx(setLog);
    c1->SetLogy(setLog);
    c1->SetName(Form("correlation_difference_unsmoothed%s", data.uinfo_tag.Data()));

    if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    c1->Clear();

    TH2D *cov_histo = new TH2D(Form("cov_eigen_histo%s", data.info_tag.Data()), "",
        covarianceDiff.GetNrows(), 1, covarianceDiff.GetNrows(),
        covarianceDiff.GetNrows(), 1, covarianceDiff.GetNrows());

    if (data.result.GetNbinsX() == covarianceDiff.GetNrows())
      cov_histo->SetBins(data.result.GetNbinsX(),
          data.result.GetXaxis()->GetXbins()->GetArray(),
          data.result.GetNbinsX(),
          data.result.GetXaxis()->GetXbins()->GetArray());

    for (int i = 1; i <= cov_histo->GetNbinsX(); ++i) {
      for (int j = i; j <= cov_histo->GetNbinsY(); ++j) {
        cov_histo->SetBinContent(i, j, 100.0 * covarianceDiff(i - 1, j - 1));
        cov_histo->SetBinContent(j, i, cov_histo->GetBinContent(i, j));
      }
    }

    cov_histo->GetXaxis()->SetTitleOffset(0.9);
    cov_histo->GetYaxis()->SetTitleOffset(0.9);

    if (data.result.GetNbinsX() == covarianceDiff.GetNrows()) {
      if (TString(data.result.GetXaxis()->GetTitle()) == "pt") {
        cov_histo->SetXTitle("p_{T} [GeV]");
        cov_histo->SetYTitle("p_{T} [GeV]");
      }
      else if (TString(data.result.GetXaxis()->GetTitle()) == "eta") {
        cov_histo->SetXTitle("#eta");
        cov_histo->SetYTitle("#eta");
      }
      else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta") {
        cov_histo->SetXTitle("#left|#eta#right|");
        cov_histo->SetYTitle("#left|#eta#right|");
      }
    }
    cov_histo->GetXaxis()->SetRange(1, covarianceDiff.GetNrows());
    cov_histo->GetYaxis()->SetRange(1, covarianceDiff.GetNrows());
    cov_histo->SetStats(kFALSE);
    cov_histo->Draw("COLZ");

    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    myText(0.2, 0.965, 0.035, kBlack, "Eigen-Variation Covariance Matrix % Relative Difference");

    c1->SetLogx(kFALSE);
    c1->SetLogy(kFALSE);

    if (data.result.GetNbinsX() == covarianceDiff.GetNrows()) {
      c1->SetLogx(setLog);
      c1->SetLogy(setLog);
    }
    c1->SetName(Form("covaraince_difference%s", data.info_tag.Data()));

    if (data.result.GetNbinsX() == covarianceDiff.GetNrows()) {
      for (int ii = 0; ii < data.uresult.GetNbinsX() - 1; ++ii) {
        double x    = data.uresult.GetBinCenter(ii + 1) + data.uresult.GetBinWidth(ii + 1) / 2.0,
               xmin = data.uresult.GetXaxis()->GetXmin(),
               xmax = data.uresult.GetXaxis()->GetXmax();
        TLine *vl   = new TLine(x, xmin, x, xmax);
        vl->Draw();
        lines.Add(vl);
        TLine *hl = new TLine(xmin, x, xmax, x);
        hl->Draw();
        lines.Add(hl);
      }
    }

    of->WriteTObject(c1);
    c1->Clear();

    TH2D *ucov_histo = new TH2D(Form("ucov_eigen_histo%s", data.info_tag.Data()), "",
        ucovarianceDiff.GetNrows(), 1, ucovarianceDiff.GetNrows(),
        ucovarianceDiff.GetNrows(), 1, ucovarianceDiff.GetNrows());

    ucov_histo->SetBins(data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray(),
        data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray());

    for (int i = 1; i <= ucov_histo->GetNbinsX(); ++i) {
      for (int j = i; j <= ucov_histo->GetNbinsY(); ++j) {
        ucov_histo->SetBinContent(i, j, 100.0 * ucovarianceDiff(i - 1, j - 1));
        ucov_histo->SetBinContent(j, i, ucov_histo->GetBinContent(i, j));
      }
    }

    ucov_histo->GetXaxis()->SetTitleOffset(0.9);
    ucov_histo->GetYaxis()->SetTitleOffset(0.9);
    if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt") {
      ucov_histo->SetXTitle("p_{T} [GeV]");
      ucov_histo->SetYTitle("p_{T} [GeV]");
    }
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta") {
      ucov_histo->SetXTitle("#eta");
      ucov_histo->SetYTitle("#eta");
    }
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta") {
      ucov_histo->SetXTitle("#left|#eta#right|");
      ucov_histo->SetYTitle("#left|#eta#right|");
    }
    ucov_histo->GetXaxis()->SetRange(1, ucovarianceDiff.GetNrows());
    ucov_histo->GetYaxis()->SetRange(1, ucovarianceDiff.GetNrows());
    ucov_histo->SetStats(kFALSE);
    ucov_histo->Draw("COLZ");

    myText(0.2, 0.05, 0.05, kBlue, data.uInfoTag.Data());
    myText(0.2, 0.965, 0.035, kBlack, "Eigen-Variation Covariance Matrix % Relative Difference (unsmoothed)");

    c1->SetLogx(setLog);
    c1->SetLogy(setLog);
    c1->SetName(Form("covaraince_difference_unsmoothed%s", data.uinfo_tag.Data()));

    if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    c1->Clear();

    TH2D *orig_cov_histo = new TH2D(Form("orig_cov_histo%s", data.info_tag.Data()), "",
        covMat.GetNrows(), 1, covMat.GetNrows(),
        covMat.GetNrows(), 1, covMat.GetNrows());


    if (data.result.GetNbinsX() == covMat.GetNrows())
      orig_cov_histo->SetBins(data.result.GetNbinsX(),
          data.result.GetXaxis()->GetXbins()->GetArray(),
          data.result.GetNbinsX(),
          data.result.GetXaxis()->GetXbins()->GetArray());

    for (int i = 1; i <= orig_cov_histo->GetNbinsX(); ++i) {
      for (int j = i; j <= orig_cov_histo->GetNbinsY(); ++j) {
        orig_cov_histo->SetBinContent(i, j, covMat(i - 1, j - 1));
        orig_cov_histo->SetBinContent(j, i, orig_cov_histo->GetBinContent(i, j));
      }
    }

    orig_cov_histo->GetXaxis()->SetTitleOffset(0.9);
    orig_cov_histo->GetYaxis()->SetTitleOffset(0.9);

    if (data.result.GetNbinsX() == covMat.GetNrows()) {
      if (TString(data.result.GetXaxis()->GetTitle()) == "pt") {
        orig_cov_histo->SetXTitle("p_{T} [GeV]");
        orig_cov_histo->SetYTitle("p_{T} [GeV]");
      }
      else if (TString(data.result.GetXaxis()->GetTitle()) == "eta") {
        orig_cov_histo->SetXTitle("#eta");
        orig_cov_histo->SetYTitle("#eta");
      }
      else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta") {
        orig_cov_histo->SetXTitle("#left|#eta#right|");
        orig_cov_histo->SetYTitle("#left|#eta#right|");
      }
    }
    orig_cov_histo->GetXaxis()->SetRange(1, covMat.GetNrows());
    orig_cov_histo->GetYaxis()->SetRange(1, covMat.GetNrows());
    orig_cov_histo->SetStats(kFALSE);
    orig_cov_histo->Draw("COLZ");

    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    myText(0.2, 0.965, 0.035, kBlack, "Original Covariance Matrix");

    c1->SetLogx(kFALSE);
    c1->SetLogy(kFALSE);

    if (data.result.GetNbinsX() == covMat.GetNrows()) {
      c1->SetLogx(setLog);
      c1->SetLogy(setLog);
    }
    c1->SetName(Form("original_covariance%s", data.info_tag.Data()));

    if (data.result.GetNbinsX() == covMat.GetNrows()) {
      for (int ii = 0; ii < data.uresult.GetNbinsX() - 1; ++ii) {
        double x    = data.uresult.GetBinCenter(ii + 1) + data.uresult.GetBinWidth(ii + 1) / 2.0,
               xmin = data.uresult.GetXaxis()->GetXmin(),
               xmax = data.uresult.GetXaxis()->GetXmax();
        TLine *vl   = new TLine(x, xmin, x, xmax);
        vl->Draw();
        lines.Add(vl);
        TLine *hl = new TLine(xmin, x, xmax, x);
        hl->Draw();
        lines.Add(hl);
      }
    }

    if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    c1->Clear();

    TH2D *orig_ucov_histo = new TH2D(Form("orig_ucov_histo%s", data.info_tag.Data()), "",
        ucovMat.GetNrows(), 1, ucovMat.GetNrows(),
        ucovMat.GetNrows(), 1, ucovMat.GetNrows());

    orig_ucov_histo->SetBins(data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray(),
        data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray());

    for (int i = 1; i <= orig_ucov_histo->GetNbinsX(); ++i) {
      for (int j = i; j <= orig_ucov_histo->GetNbinsY(); ++j) {
        orig_ucov_histo->SetBinContent(i, j, ucovMat(i - 1, j - 1));
        orig_ucov_histo->SetBinContent(j, i, orig_ucov_histo->GetBinContent(i, j));
      }
    }

    orig_ucov_histo->GetXaxis()->SetTitleOffset(0.9);
    orig_ucov_histo->GetYaxis()->SetTitleOffset(0.9);
    if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt") {
      orig_ucov_histo->SetXTitle("p_{T} [GeV]");
      orig_ucov_histo->SetYTitle("p_{T} [GeV]");
    }
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta") {
      orig_ucov_histo->SetXTitle("#eta");
      orig_ucov_histo->SetYTitle("#eta");
    }
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta") {
      orig_ucov_histo->SetXTitle("#left|#eta#right|");
      orig_ucov_histo->SetYTitle("#left|#eta#right|");
    }
    orig_ucov_histo->GetXaxis()->SetRange(1, ucovMat.GetNrows());
    orig_ucov_histo->GetYaxis()->SetRange(1, ucovMat.GetNrows());
    orig_ucov_histo->SetStats(kFALSE);
    orig_ucov_histo->Draw("COLZ");

    myText(0.2, 0.05, 0.05, kBlue, data.uInfoTag.Data());
    myText(0.2, 0.965, 0.035, kBlack, "Original Covariance Matrix (unsmoothed)");

    c1->SetLogx(setLog);
    c1->SetLogy(setLog);
    c1->SetName(Form("original_covariance_unsmoothed%s", data.uinfo_tag.Data()));

    if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    c1->Clear();

    TH2D *eigen_cov_histo = new TH2D(Form("eigen_cov_histo%s", data.info_tag.Data()), "",
        covMatVar.GetNrows(), 1, covMatVar.GetNrows(),
        covMatVar.GetNrows(), 1, covMatVar.GetNrows());

    if (data.result.GetNbinsX() == covMatVar.GetNrows())
      eigen_cov_histo->SetBins(data.result.GetNbinsX(),
          data.result.GetXaxis()->GetXbins()->GetArray(),
          data.result.GetNbinsX(),
          data.result.GetXaxis()->GetXbins()->GetArray());

    for (int i = 1; i <= eigen_cov_histo->GetNbinsX(); ++i) {
      for (int j = i; j <= eigen_cov_histo->GetNbinsY(); ++j) {
        eigen_cov_histo->SetBinContent(i, j, covMatVar(i - 1, j - 1));
        eigen_cov_histo->SetBinContent(j, i, eigen_cov_histo->GetBinContent(i, j));
      }
    }

    eigen_cov_histo->GetXaxis()->SetTitleOffset(0.9);
    eigen_cov_histo->GetYaxis()->SetTitleOffset(0.9);

    if (data.result.GetNbinsX() == covMatVar.GetNrows()) {
      if (TString(data.result.GetXaxis()->GetTitle()) == "pt") {
        eigen_cov_histo->SetXTitle("p_{T} [GeV]");
        eigen_cov_histo->SetYTitle("p_{T} [GeV]");
      }
      else if (TString(data.result.GetXaxis()->GetTitle()) == "eta") {
        eigen_cov_histo->SetXTitle("#eta");
        eigen_cov_histo->SetYTitle("#eta");
      }
      else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta") {
        eigen_cov_histo->SetXTitle("#left|#eta#right|");
        eigen_cov_histo->SetYTitle("#left|#eta#right|");
      }
    }
    eigen_cov_histo->GetXaxis()->SetRange(1, covMatVar.GetNrows());
    eigen_cov_histo->GetYaxis()->SetRange(1, covMatVar.GetNrows());
    eigen_cov_histo->SetStats(kFALSE);
    eigen_cov_histo->Draw("COLZ");

    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    myText(0.2, 0.965, 0.035, kBlack, "Eigen-Variation Covariance Matrix");

    c1->SetLogx(kFALSE);
    c1->SetLogy(kFALSE);

    if (data.result.GetNbinsX() == covMatVar.GetNrows()) {
      c1->SetLogx(setLog);
      c1->SetLogy(setLog);
    }
    c1->SetName(Form("eigen_covariance%s", data.info_tag.Data()));

    if (data.result.GetNbinsX() == covMatVar.GetNrows()) {
      for (int ii = 0; ii < data.uresult.GetNbinsX() - 1; ++ii) {
        double x    = data.uresult.GetBinCenter(ii + 1) + data.uresult.GetBinWidth(ii + 1) / 2.0,
               xmin = data.uresult.GetXaxis()->GetXmin(),
               xmax = data.uresult.GetXaxis()->GetXmax();
        TLine *vl   = new TLine(x, xmin, x, xmax);
        vl->Draw();
        lines.Add(vl);
        TLine *hl = new TLine(xmin, x, xmax, x);
        hl->Draw();
        lines.Add(hl);
      }
    }

    if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    c1->Clear();

    TH2D *eigen_ucov_histo = new TH2D(Form("eigen_ucov_histo%s", data.info_tag.Data()), "",
        ucovMatVar.GetNrows(), 1, ucovMatVar.GetNrows(),
        ucovMatVar.GetNrows(), 1, ucovMatVar.GetNrows());

    eigen_ucov_histo->SetBins(data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray(),
        data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray());

    for (int i = 1; i <= eigen_ucov_histo->GetNbinsX(); ++i) {
      for (int j = i; j <= eigen_ucov_histo->GetNbinsY(); ++j) {
        eigen_ucov_histo->SetBinContent(i, j, ucovMatVar(i - 1, j - 1));
        eigen_ucov_histo->SetBinContent(j, i, eigen_ucov_histo->GetBinContent(i, j));
      }
    }

    eigen_ucov_histo->GetXaxis()->SetTitleOffset(0.9);
    eigen_ucov_histo->GetYaxis()->SetTitleOffset(0.9);
    if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt") {
      eigen_ucov_histo->SetXTitle("p_{T} [GeV]");
      eigen_ucov_histo->SetYTitle("p_{T} [GeV]");
    }
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta") {
      eigen_ucov_histo->SetXTitle("#eta");
      eigen_ucov_histo->SetYTitle("#eta");
    }
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta") {
      eigen_ucov_histo->SetXTitle("#left|#eta#right|");
      eigen_ucov_histo->SetYTitle("#left|#eta#right|");
    }
    eigen_ucov_histo->GetXaxis()->SetRange(1, ucovMatVar.GetNrows());
    eigen_ucov_histo->GetYaxis()->SetRange(1, ucovMatVar.GetNrows());
    eigen_ucov_histo->SetStats(kFALSE);
    eigen_ucov_histo->Draw("COLZ");

    myText(0.2, 0.05, 0.05, kBlue, data.uInfoTag.Data());
    myText(0.2, 0.965, 0.035, kBlack, "Eigen-Variance Covariance Matrix (unsmoothed)");

    c1->SetLogx(setLog);
    c1->SetLogy(setLog);
    c1->SetName(Form("eigen_covariance_unsmoothed%s", data.uinfo_tag.Data()));

    if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    c1->Clear();

    TH1D *error_diff = new TH1D(Form("error_diff%s", data.info_tag.Data()), "",
        covMat.GetNrows(),
        1, covMat.GetNrows());
    error_diff->SetBins(data.result.GetNbinsX(),
        data.result.GetXaxis()->GetXbins()->GetArray());

    for (int i = 0; i < covMat.GetNrows(); ++i) {
      if (covMat(i, i) != 0.0) error_diff->SetBinContent(i + 1, 100.0 * (sqrt(covMatVar(i, i)) / sqrt(covMat(i, i)) - 1));
      else error_diff->SetBinContent(i + 1, 0.0);
    }

    error_diff->SetStats(kFALSE);
    error_diff->GetXaxis()->SetTitleOffset(0.9);
    error_diff->GetYaxis()->SetTitleOffset(1.2);
    error_diff->GetYaxis()->SetTitle("% Relative Difference in Total Error");
    if (TString(data.result.GetXaxis()->GetTitle()) == "pt")
      error_diff->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (TString(data.result.GetXaxis()->GetTitle()) == "eta")
      error_diff->GetXaxis()->SetTitle("#eta");
    else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta")
      error_diff->GetXaxis()->SetTitle("#left|#eta#right|");
    error_diff->Draw("hist");

    myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    myText(0.25, 0.25, 0.035, kBlue, pathName.Data());

    c1->SetLogx(setLog);
    c1->SetLogy(kFALSE);
    c1->SetName(error_diff->GetName());
    of->WriteTObject(c1);
    c1->Clear();

    TH1D *uerror_diff = new TH1D(Form("error_diff_unsmoothed_histo%s", data.uinfo_tag.Data()), "",
        ucovMat.GetNrows(),
        1, ucovMat.GetNrows());
    uerror_diff->SetBins(data.uresult.GetNbinsX(),
        data.uresult.GetXaxis()->GetXbins()->GetArray());

    for (int i = 0; i < ucovMat.GetNrows(); ++i) {
      if (ucovMat(i, i) != 0.0) uerror_diff->SetBinContent(i + 1, 100.0 * (sqrt(ucovMatVar(i, i)) / sqrt(ucovMat(i, i)) - 1));
      else uerror_diff->SetBinContent(i + 1, 0.0);
    }

    uerror_diff->SetStats(kFALSE);
    uerror_diff->GetXaxis()->SetTitleOffset(0.9);
    uerror_diff->GetYaxis()->SetTitleOffset(1.2);
    uerror_diff->GetYaxis()->SetTitle("% Relative Difference in Total Error");
    if (TString(data.uresult.GetXaxis()->GetTitle()) == "pt")
      uerror_diff->GetXaxis()->SetTitle("p_{T} [GeV]");
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "eta")
      uerror_diff->GetXaxis()->SetTitle("#eta");
    else if (TString(data.uresult.GetXaxis()->GetTitle()) == "abseta")
      uerror_diff->GetXaxis()->SetTitle("#left|#eta#right|");
    uerror_diff->Draw("hist");

    myText(0.2, 0.05, 0.05, kBlue, data.uInfoTag.Data());
    myText(0.25, 0.25, 0.035, kBlue, pathName.Data());

    c1->SetLogx(setLog);
    c1->SetLogy(kFALSE);
    c1->SetName(Form("error_diff_unsmoothed%s", data.uinfo_tag.Data()));
    if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    c1->Clear();

    // // ////////////////////////////////////////////////////////
    // // TODO: START TESTING
    // // ////////////////////////////////////////////////////////
    // // RMinimizer func(covMat, 1, "Minuit2", "Simplex");
    // RMinimizer func(covMat, 2, "Minuit2", "Simplex");
    // // RAmoeba aMin(covMat, 5);
    //
    // double x[1] = {10.0};
    // func(x);
    // // auto X = aMin(0.001);
    //
    // if (func.Successful()) {
    // // if (true) {
    //   auto covMatReduced = func.GetFunction().ComputeCovarianceMatrix(func.GetMinimizer().X());
    //   // auto covMatReduced = aMin.GetFunction().ComputeCovarianceMatrix(X);
    //
    //   TH2D *reduced_cov_histo = new TH2D(Form("reduced_cov_histo%s", data.info_tag.Data()), "",
    //   covMatReduced.GetNrows(), 1, covMatReduced.GetNrows(),
    //   covMatReduced.GetNrows(), 1, covMatReduced.GetNrows());
    //
    //
    //   if (data.result.GetNbinsX() == covMatReduced.GetNrows())
    //   reduced_cov_histo->SetBins(data.result.GetNbinsX(),
    //   data.result.GetXaxis()->GetXbins()->GetArray(),
    //   data.result.GetNbinsX(),
    //   data.result.GetXaxis()->GetXbins()->GetArray());
    //
    //   for (int i = 1; i <= reduced_cov_histo->GetNbinsX(); ++i) {
    //     for (int j = i; j <= reduced_cov_histo->GetNbinsY(); ++j) {
    //       reduced_cov_histo->SetBinContent(i, j, covMatReduced(i - 1, j - 1));
    //       reduced_cov_histo->SetBinContent(j, i, reduced_cov_histo->GetBinContent(i, j));
    //     }
    //   }
    //
    //   reduced_cov_histo->GetXaxis()->SetTitleOffset(0.9);
    //   reduced_cov_histo->GetYaxis()->SetTitleOffset(0.9);
    //
    //   if (data.result.GetNbinsX() == covMatReduced.GetNrows()) {
    //     if (TString(data.result.GetXaxis()->GetTitle()) == "pt") {
    //       reduced_cov_histo->SetXTitle("p_{T} [GeV]");
    //       reduced_cov_histo->SetYTitle("p_{T} [GeV]");
    //     }
    //     else if (TString(data.result.GetXaxis()->GetTitle()) == "eta") {
    //       reduced_cov_histo->SetXTitle("#eta");
    //       reduced_cov_histo->SetYTitle("#eta");
    //     }
    //     else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta") {
    //       reduced_cov_histo->SetXTitle("#left|#eta#right|");
    //       reduced_cov_histo->SetYTitle("#left|#eta#right|");
    //     }
    //   }
    //   reduced_cov_histo->GetXaxis()->SetRange(1, covMatReduced.GetNrows());
    //   reduced_cov_histo->GetYaxis()->SetRange(1, covMatReduced.GetNrows());
    //   reduced_cov_histo->SetStats(kFALSE);
    //   reduced_cov_histo->Draw("COLZ");
    //
    //   myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    //   myText(0.2, 0.965, 0.035, kBlack, "Reduced Covariance Matrix");
    //
    //   c1->SetLogx(kFALSE);
    //   c1->SetLogy(kFALSE);
    //
    //   if (data.result.GetNbinsX() == covMatReduced.GetNrows()) {
    //     c1->SetLogx(setLog);
    //     c1->SetLogy(setLog);
    //   }
    //   c1->SetName(Form("reduced_covariance%s", data.info_tag.Data()));
    //
    //   if (data.result.GetNbinsX() == covMatReduced.GetNrows()) {
    //     for (int ii = 0; ii < data.uresult.GetNbinsX() - 1; ++ii) {
    //       double x    = data.uresult.GetBinCenter(ii + 1) + data.uresult.GetBinWidth(ii + 1) / 2.0,
    //       xmin = data.uresult.GetXaxis()->GetXmin(),
    //       xmax = data.uresult.GetXaxis()->GetXmax();
    //       TLine *vl   = new TLine(x, xmin, x, xmax);
    //       vl->Draw();
    //       lines.Add(vl);
    //       TLine *hl = new TLine(xmin, x, xmax, x);
    //       hl->Draw();
    //       lines.Add(hl);
    //     }
    //   }
    //
    //   if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    //   c1->Clear();
    //
    //   covMatReduced = func.GetFunction().ComputeRelativeDifferenceCovarianceMatrix(func.GetMinimizer().X());
    //   // covMatReduced = aMin.GetFunction().ComputeRelativeDifferenceCovarianceMatrix(X);
    //
    //   TH2D *reduced_cov_diff_histo = new TH2D(Form("reduced_cov_diff_histo%s", data.info_tag.Data()), "",
    //   covMatReduced.GetNrows(), 1, covMatReduced.GetNrows(),
    //   covMatReduced.GetNrows(), 1, covMatReduced.GetNrows());
    //
    //
    //   if (data.result.GetNbinsX() == covMatReduced.GetNrows())
    //   reduced_cov_diff_histo->SetBins(data.result.GetNbinsX(),
    //   data.result.GetXaxis()->GetXbins()->GetArray(),
    //   data.result.GetNbinsX(),
    //   data.result.GetXaxis()->GetXbins()->GetArray());
    //
    //   for (int i = 1; i <= reduced_cov_diff_histo->GetNbinsX(); ++i) {
    //     for (int j = i; j <= reduced_cov_diff_histo->GetNbinsY(); ++j) {
    //       reduced_cov_diff_histo->SetBinContent(i, j, 100.0 * covMatReduced(i - 1, j - 1));
    //       reduced_cov_diff_histo->SetBinContent(j, i, reduced_cov_diff_histo->GetBinContent(i, j));
    //     }
    //   }
    //
    //   reduced_cov_diff_histo->GetXaxis()->SetTitleOffset(0.9);
    //   reduced_cov_diff_histo->GetYaxis()->SetTitleOffset(0.9);
    //
    //   if (data.result.GetNbinsX() == covMatReduced.GetNrows()) {
    //     if (TString(data.result.GetXaxis()->GetTitle()) == "pt") {
    //       reduced_cov_diff_histo->SetXTitle("p_{T} [GeV]");
    //       reduced_cov_diff_histo->SetYTitle("p_{T} [GeV]");
    //     }
    //     else if (TString(data.result.GetXaxis()->GetTitle()) == "eta") {
    //       reduced_cov_diff_histo->SetXTitle("#eta");
    //       reduced_cov_diff_histo->SetYTitle("#eta");
    //     }
    //     else if (TString(data.result.GetXaxis()->GetTitle()) == "abseta") {
    //       reduced_cov_diff_histo->SetXTitle("#left|#eta#right|");
    //       reduced_cov_diff_histo->SetYTitle("#left|#eta#right|");
    //     }
    //   }
    //   reduced_cov_diff_histo->GetXaxis()->SetRange(1, covMatReduced.GetNrows());
    //   reduced_cov_diff_histo->GetYaxis()->SetRange(1, covMatReduced.GetNrows());
    //   reduced_cov_diff_histo->SetStats(kFALSE);
    //   reduced_cov_diff_histo->Draw("COLZ");
    //
    //   myText(0.2, 0.05, 0.05, kBlue, data.InfoTag.Data());
    //   myText(0.2, 0.965, 0.035, kBlack, "Reduced Covariance Matrix % Relative Difference");
    //
    //   c1->SetLogx(kFALSE);
    //   c1->SetLogy(kFALSE);
    //
    //   if (data.result.GetNbinsX() == covMatReduced.GetNrows()) {
    //     c1->SetLogx(setLog);
    //     c1->SetLogy(setLog);
    //   }
    //   c1->SetName(Form("reduced_covariance_diff%s", data.info_tag.Data()));
    //
    //   if (data.result.GetNbinsX() == covMatReduced.GetNrows()) {
    //     for (int ii = 0; ii < data.uresult.GetNbinsX() - 1; ++ii) {
    //       double x    = data.uresult.GetBinCenter(ii + 1) + data.uresult.GetBinWidth(ii + 1) / 2.0,
    //       xmin = data.uresult.GetXaxis()->GetXmin(),
    //       xmax = data.uresult.GetXaxis()->GetXmax();
    //       TLine *vl   = new TLine(x, xmin, x, xmax);
    //       vl->Draw();
    //       lines.Add(vl);
    //       TLine *hl = new TLine(xmin, x, xmax, x);
    //       hl->Draw();
    //       lines.Add(hl);
    //     }
    //   }
    //
    //   if (!of->Get(c1->GetName())) of->WriteTObject(c1);
    //   c1->Clear();
    // }
    //
    // // ////////////////////////////////////////////////////////
    // // TODO: END TESTING
    // // ////////////////////////////////////////////////////////



    c1->SetRightMargin(old_right_margin);
    c1->SetLeftMargin(old_left_margin);
    c1->SetBottomMargin(old_bottom_margin);

    // ////////////////////////////////////////////////////////
    // Write Out Bulk Statistics to Data File (smoothed only)
    // ////////////////////////////////////////////////////////
    TString data_file_name(TString::Format("eigen_data_%lu_var_%s%s_%s.txt",
          data.eigen_set.size(),
          tag.Data(),
          data.info_tag.Data(),
          pathName.Copy().ReplaceAll("/", "-").Data()).Data());
    std::ofstream myfile(data_file_name.Data());

    Info("eigenDecompCalib", "statistical data will be stored in %s", data_file_name.Data());

    if (myfile.is_open()) {
      double super_maxdiff    = 0;
      double diagonal_maxdiff = 0;
      double avediff          = 0.0;
      long   denom_avediff    = 0;

      for (int u = 0; u < covarianceDiff.GetNrows(); ++u) {
        double maxdiff = 0.0;

        for (int v = 0; v < covarianceDiff.GetNcols(); ++v) {
          double diff = covarianceDiff(u, v);

          if (fabs(diff) > maxdiff) maxdiff = fabs(diff);

          if ((v == u) && (fabs(diff) > diagonal_maxdiff)) diagonal_maxdiff = fabs(diff);
          avediff += diff;
          ++denom_avediff;
        }

        if (maxdiff > 0.0) {
          if (super_maxdiff < maxdiff) super_maxdiff = maxdiff;
        }
      }

      if (denom_avediff != 0) myfile << "The average relative difference in the covariance matrix is: " << avediff / denom_avediff << std::endl;
      myfile << "The largest relative difference in the covariance matrix is: " << super_maxdiff << std::endl;
      myfile << "The largest relative difference in the diagonal elements of the covariance matrix is: " << diagonal_maxdiff << std::endl;

      super_maxdiff = 0;
      avediff       = 0.0;
      denom_avediff = 0;

      for (int u = 0; u < correlationDiff.GetNrows(); ++u) {
        for (int v = 0; v < correlationDiff.GetNcols(); ++v) {
          if (super_maxdiff < fabs(correlationDiff(u, v))) super_maxdiff = fabs(correlationDiff(u, v));

          if (fabs(correlationDiff(u, v)) != 0.0) {
            avediff += correlationDiff(u, v);
            ++denom_avediff;
          }
        }
      }

      if (denom_avediff != 0) myfile << "The average relative difference in the correlation matrix is: " << avediff / denom_avediff << std::endl;
      myfile << "The largest relative difference in the correlation matrix is: " << super_maxdiff << std::endl;
      myfile << "The largest relative difference in the diagonal elements of the correlation matrix is: " << diagonal_maxdiff << std::endl;

      super_maxdiff = 0;
      avediff       = 0.0;
      denom_avediff = 0;

      for (int i = 0; i < covMat.GetNrows(); ++i) {
        double diff = fabs(error_diff->GetBinContent(i));

        if (super_maxdiff < diff) super_maxdiff = diff;
        avediff += diff;
        ++denom_avediff;
      }

      if (denom_avediff != 0) myfile << "The average relative difference in the total error is: " << avediff / denom_avediff << std::endl;
      myfile << "The largest relative difference in the total error is: " << super_maxdiff << std::endl;


      myfile.close();
    }
  }
}

void MakeDataFrom2D (Int_t xbin,
    TH1 *result, TH1 *uresult,
    CalibrationDataHistogramContainer *c, CalibrationDataHistogramContainer *uc,
    CalibrationDataEigenVariations *eigen, CalibrationDataEigenVariations *ueigen,
    TMatrixDSym &original_cov, TMatrixDSym &uoriginal_cov,
    TMatrixDSym &eigen_cov, TMatrixDSym &ueigen_cov,
    std::function<const TAxis*(const TH1&)> GetAxis,
    std::function<Int_t(const TH1&)> GetNbinsAxis,
    std::function<Int_t(const TH1&, Int_t, Int_t, Int_t)> GetBin,
    std::function<TH1D*(TH2&, const char*, Int_t, Int_t, Option_t*)> Project,
    std::function<void(AnalysisData&)> RunAnalysis)
{
  auto xaxis      = GetAxis(*result),
       uxaxis     = GetAxis(*uresult);
  auto xcenter    = xaxis->GetBinCenter(xbin);
  auto uxbin      = uxaxis->FindFixBin(xcenter);
  auto uxcenter   = uxaxis->GetBinCenter(uxbin);
  TString info_tag(Form("_%s_%.3f", xaxis->GetTitle(), xcenter)),
          uinfo_tag(Form("_%s_%.3f", uxaxis->GetTitle(), uxcenter));
  TString InfoTag, uInfoTag;
  TH1D *result0   = Project(*static_cast<TH2*>(result), Form("smoothed_nominal%s", info_tag.Data()), xbin, xbin, ""),
       *uresult0  = Project(*static_cast<TH2*>(uresult), Form("unsmoothed_nominal%s", info_tag.Data()), uxbin, uxbin, "");
  vector<TH1*> stat_set, ustat_set, syst_set, usyst_set, eigen_set, ueigen_set;
  typename std::remove_reference<decltype(original_cov)>::type original_cov0;
  typename std::remove_reference<decltype(uoriginal_cov)>::type uoriginal_cov0;
  typename std::remove_reference<decltype(eigen_cov)>::type eigen_cov0;
  typename std::remove_reference<decltype(ueigen_cov)>::type ueigen_cov0;

  if (TString(xaxis->GetTitle()) == "abseta")
    InfoTag = Form("#left|#eta#right| = %.3f", xcenter);
  else if (TString(xaxis->GetTitle()) == "eta")
    InfoTag = Form("#eta = %.3f", xcenter);
  else if (TString(xaxis->GetTitle()) == "pt")
    InfoTag = Form("p_{T} = %.3f", xcenter);
  else
    InfoTag = Form("%s = %.3f", xaxis->GetTitle(), xcenter);

  if (TString(uxaxis->GetTitle()) == "abseta")
    uInfoTag = Form("#left|#eta#right| = %.3f", uxcenter);
  else if (TString(uxaxis->GetTitle()) == "eta")
    uInfoTag = Form("#eta = %.3f", uxcenter);
  else if (TString(uxaxis->GetTitle()) == "pt")
    uInfoTag = Form("p_{T} = %.3f", uxcenter);
  else
    uInfoTag = Form("%s = %.3f", uxaxis->GetTitle(), uxcenter);

  // construct stat, syst, and eigen sets
  // ustat, usyst
  {
    auto utemp_stat = static_cast<decltype(uresult0)>(uresult0->Clone());
    utemp_stat->SetDirectory(0);
    for (auto bin = 1; bin <= utemp_stat->GetNbinsX(); ++bin) utemp_stat->SetBinContent(bin, utemp_stat->GetBinError(bin));
    ustat_set.push_back(utemp_stat);

    for (auto un : uc->listUncertainties()) {
      if ((un == "result") ||
          (un == "extrapolation") ||
          (un == "systematics") ||
          !(*uc)(un.c_str())->InheritsFrom("TH1")) continue;
      auto syst = static_cast<TH2*>((*uc)(un.c_str()));
      auto h = Project(*syst, Form("%s%s", syst->GetName(), info_tag.Data()), uxbin, uxbin, "");
      h->SetDirectory(0);

      if (!uc->isBinCorrelated(un.c_str()))
        ustat_set.push_back(h);
      else
        usyst_set.push_back(h);
    }
  }
  // stat, syst
  {
    for (auto un : c->listUncertainties()) {
      if ((un == "result") ||
          (un == "extrapolation") ||
          (un == "systematics") ||
          !(*c)(un.c_str())->InheritsFrom("TH1")) continue;
      auto syst = static_cast<TH2*>((*c)(un.c_str()));
      auto h = Project(*syst, Form("%s%s", syst->GetName(), info_tag.Data()), xbin, xbin, "");
      h->SetDirectory(0);
      if ((un.find("stat_nuis") != string::npos) || (un.find("Stat Nuis") != string::npos) || (un.find("Nuis ") != string::npos))
        stat_set.push_back(h);
      else
        syst_set.push_back(h);
    }
  }
  // eigen ueigen
  {
    TH1 *up = nullptr, *down = nullptr,
        *uup = nullptr, *udown = nullptr;
    for (decltype(ueigen->getNumberOfEigenVariations()) i = 0; i < ueigen->getNumberOfEigenVariations(); ++i) {
      if (ueigen->getEigenvectorVariation(i, uup, udown)) {
        TH1 *up0 = Project(*static_cast<TH2*>(uup), Form("%s_%d%s", uup->GetTitle(), i, uinfo_tag.Data()), uxbin, uxbin, "");
        up0->Add(uresult0, -1.0);
        up0->SetDirectory(0);
        ueigen_set.push_back(up0);
      }
    }

    for (decltype(eigen->getNumberOfEigenVariations()) i = 0; i < eigen->getNumberOfEigenVariations(); ++i) {
      if (eigen->getEigenvectorVariation(i, up, down)) {
        TH1 *up0 = Project(*static_cast<TH2*>(up), Form("%s_%d%s", up->GetTitle(), i, info_tag.Data()), xbin, xbin, "");
        up0->Add(result0, -1.0);
        up0->SetDirectory(0);
        eigen_set.push_back(up0);
      }
    }
  }

  // construct reduced covariance matrix
  {
    vector<decltype(result->GetNbinsX())> valid_bins;
    for (auto ybin = 1; ybin <= GetNbinsAxis(*result); ++ybin)
      // valid_bins.push_back(result->GetBin(xbin, ybin));
      valid_bins.push_back(GetBin(*result, xbin, ybin, 0));
    original_cov0.ResizeTo(valid_bins.size(), valid_bins.size());
    eigen_cov0.ResizeTo(valid_bins.size(), valid_bins.size());
    for (decltype(valid_bins.size()) i = 0; i < valid_bins.size(); ++i) {
      for (decltype(valid_bins.size()) j = 0; j < valid_bins.size(); ++j) {
        original_cov0(i, j) = original_cov(valid_bins[i], valid_bins[j]);
        original_cov0(j, i) = original_cov(valid_bins[j], valid_bins[i]);
        eigen_cov0(i, j) = eigen_cov(valid_bins[i], valid_bins[j]);
        eigen_cov0(j, i) = eigen_cov(valid_bins[j], valid_bins[i]);
      }
    }
    valid_bins.clear();
    for (auto uybin = 1; uybin <= GetNbinsAxis(*uresult); ++uybin)
      // valid_bins.push_back(uresult->GetBin(uxbin, uybin));
      valid_bins.push_back(GetBin(*uresult, uxbin, uybin, 0));
    uoriginal_cov0.ResizeTo(valid_bins.size(), valid_bins.size());
    ueigen_cov0.ResizeTo(valid_bins.size(), valid_bins.size());
    for (decltype(valid_bins.size()) i = 0; i < valid_bins.size(); ++i) {
      for (decltype(valid_bins.size()) j = 0; j < valid_bins.size(); ++j) {
        uoriginal_cov0(i, j) = uoriginal_cov(valid_bins[i], valid_bins[j]);
        uoriginal_cov0(j, i) = uoriginal_cov(valid_bins[j], valid_bins[i]);
        ueigen_cov0(i, j) = ueigen_cov(valid_bins[i], valid_bins[j]);
        ueigen_cov0(j, i) = ueigen_cov(valid_bins[j], valid_bins[i]);
      }
    }
  }

  AnalysisData data{
    info_tag,
    uinfo_tag,
    InfoTag,
    uInfoTag,
    *result0,
    *uresult0,
    stat_set,
    ustat_set,
    syst_set,
    usyst_set,
    eigen_set,
    ueigen_set,
    original_cov0,
    uoriginal_cov0,
    eigen_cov0,
    ueigen_cov0
  };
  RunAnalysis(data);
}

void MakeDataFrom3D (Int_t xbin, Int_t ybin,
    TH1 *result, TH1 *uresult,
    CalibrationDataHistogramContainer *c, CalibrationDataHistogramContainer *uc,
    CalibrationDataEigenVariations *eigen, CalibrationDataEigenVariations *ueigen,
    TMatrixDSym &original_cov, TMatrixDSym &uoriginal_cov,
    TMatrixDSym &eigen_cov, TMatrixDSym &ueigen_cov,
    std::function<const TAxis*(const TH1&)> GetXAxis,
    std::function<const TAxis*(const TH1&)> GetYAxis,
    std::function<Int_t(const TH1&)> GetNbinsAxis,
    std::function<Int_t(const TH1&, Int_t, Int_t, Int_t)> GetBin,
    std::function<TH1D*(TH3&, const char*, Int_t, Int_t, Int_t, Int_t)> Project,
    std::function<void(AnalysisData&)> RunAnalysis)
{
  auto xaxis      = GetXAxis(*result),
       uxaxis     = GetXAxis(*uresult);
  auto yaxis      = GetYAxis(*result),
       uyaxis     = GetYAxis(*uresult);
  auto xcenter    = xaxis->GetBinCenter(xbin);
  auto uxbin      = uxaxis->FindFixBin(xcenter);
  auto uxcenter   = uxaxis->GetBinCenter(uxbin);
  auto ycenter    = yaxis->GetBinCenter(ybin);
  auto uybin      = uyaxis->FindFixBin(ycenter);
  auto uycenter   = uyaxis->GetBinCenter(uybin);
  TString info_tag(Form("_%s_%.3f_%s_%.3f", xaxis->GetTitle(), xcenter, yaxis->GetTitle(), ycenter)),
          uinfo_tag(Form("_%s_%.3f_%s_%.3f", uxaxis->GetTitle(), uxcenter, uyaxis->GetTitle(), uycenter));
  TString InfoTag, uInfoTag;
  TH1D *result0   = Project(*static_cast<TH3*>(result), Form("smoothed_nominal%s", info_tag.Data()), xbin, xbin, ybin, ybin),
       *uresult0  = Project(*static_cast<TH3*>(uresult), Form("unsmoothed_nominal%s", info_tag.Data()), uxbin, uxbin, uybin, uybin);
  vector<TH1*> stat_set, ustat_set, syst_set, usyst_set, eigen_set, ueigen_set;
  typename std::remove_reference<decltype(original_cov)>::type original_cov0;
  typename std::remove_reference<decltype(uoriginal_cov)>::type uoriginal_cov0;
  typename std::remove_reference<decltype(eigen_cov)>::type eigen_cov0;
  typename std::remove_reference<decltype(ueigen_cov)>::type ueigen_cov0;

  if (TString(xaxis->GetTitle()) == "abseta")
    InfoTag = Form("#left|#eta#right| = %.3f", xcenter);
  else if (TString(xaxis->GetTitle()) == "eta")
    InfoTag = Form("#eta = %.3f", xcenter);
  else if (TString(xaxis->GetTitle()) == "pt")
    InfoTag = Form("p_{T} = %.3f", xcenter);
  else
    InfoTag = Form("%s = %.3f", xaxis->GetTitle(), xcenter);

  if (TString(uxaxis->GetTitle()) == "abseta")
    uInfoTag = Form("#left|#eta#right| = %.3f", uxcenter);
  else if (TString(uxaxis->GetTitle()) == "eta")
    uInfoTag = Form("#eta = %.3f", uxcenter);
  else if (TString(uxaxis->GetTitle()) == "pt")
    uInfoTag = Form("p_{T} = %.3f", uxcenter);
  else
    uInfoTag = Form("%s = %.3f", uxaxis->GetTitle(), uxcenter);

  if (TString(yaxis->GetTitle()) == "abseta")
    InfoTag = Form("%s, #left|#eta#right| = %.3f", InfoTag.Data(), ycenter);
  else if (TString(yaxis->GetTitle()) == "eta")
    InfoTag = Form("%s, #eta = %.3f", InfoTag.Data(), ycenter);
  else if (TString(yaxis->GetTitle()) == "pt")
    InfoTag = Form("%s, p_{T} = %.3f", InfoTag.Data(), ycenter);
  else
    InfoTag = Form("%s, %s = %.3f", InfoTag.Data(), yaxis->GetTitle(), ycenter);

  if (TString(uyaxis->GetTitle()) == "abseta")
    uInfoTag = Form("%s, #left|#eta#right| = %.3f", InfoTag.Data(), uycenter);
  else if (TString(uyaxis->GetTitle()) == "eta")
    uInfoTag = Form("%s, #eta = %.3f", InfoTag.Data(), uycenter);
  else if (TString(uyaxis->GetTitle()) == "pt")
    uInfoTag = Form("%s, p_{T} = %.3f", InfoTag.Data(), uycenter);
  else
    uInfoTag = Form("%s, %s = %.3f", InfoTag.Data(), uyaxis->GetTitle(), uycenter);

  // construct stat, syst, and eigen sets
  // ustat, usyst
  {
    auto utemp_stat = static_cast<decltype(uresult0)>(uresult0->Clone());
    utemp_stat->SetDirectory(0);
    for (auto bin = 1; bin <= utemp_stat->GetNbinsX(); ++bin) utemp_stat->SetBinContent(bin, utemp_stat->GetBinError(bin));
    ustat_set.push_back(utemp_stat);

    for (auto un : uc->listUncertainties()) {
      if ((un == "result") ||
          (un == "extrapolation") ||
          (un == "systematics") ||
          !(*uc)(un.c_str())->InheritsFrom("TH1")) continue;
      auto syst = static_cast<TH3*>((*uc)(un.c_str()));
      auto h = Project(*syst, Form("%s%s", syst->GetName(), info_tag.Data()), uxbin, uxbin, uybin, uybin);
      h->SetDirectory(0);

      if (!uc->isBinCorrelated(un.c_str()))
        ustat_set.push_back(h);
      else
        usyst_set.push_back(h);
    }
  }
  // stat, syst
  {
    for (auto un : c->listUncertainties()) {
      if ((un == "result") ||
          (un == "extrapolation") ||
          (un == "systematics") ||
          !(*c)(un.c_str())->InheritsFrom("TH1")) continue;
      auto syst = static_cast<TH3*>((*c)(un.c_str()));
      auto h = Project(*syst, Form("%s%s", syst->GetName(), info_tag.Data()), xbin, xbin, ybin, ybin);
      h->SetDirectory(0);
      if ((un.find("stat_nuis") != string::npos) || (un.find("Stat Nuis") != string::npos) || (un.find("Nuis ") != string::npos))
        stat_set.push_back(h);
      else
        syst_set.push_back(h);
    }
  }
  // eigen ueigen
  {
    TH1 *up = nullptr, *down = nullptr,
        *uup = nullptr, *udown = nullptr;
    for (decltype(ueigen->getNumberOfEigenVariations()) i = 0; i < ueigen->getNumberOfEigenVariations(); ++i) {
      if (ueigen->getEigenvectorVariation(i, uup, udown)) {
        TH1 *up0 = Project(*static_cast<TH3*>(uup), Form("%s_%d%s", uup->GetTitle(), i, uinfo_tag.Data()), uxbin, uxbin, uybin, uybin);
        up0->Add(uresult0, -1.0);
        up0->SetDirectory(0);
        ueigen_set.push_back(up0);
      }
    }

    for (decltype(eigen->getNumberOfEigenVariations()) i = 0; i < eigen->getNumberOfEigenVariations(); ++i) {
      if (eigen->getEigenvectorVariation(i, up, down)) {
        TH1 *up0 = Project(*static_cast<TH3*>(up), Form("%s_%d%s", up->GetTitle(), i, info_tag.Data()), xbin, xbin, ybin, ybin);
        up0->Add(result0, -1.0);
        up0->SetDirectory(0);
        eigen_set.push_back(up0);
      }
    }
  }

  // construct reduced covariance matrix
  {
    vector<decltype(result->GetNbinsX())> valid_bins;
    for (auto zbin = 1; zbin <= GetNbinsAxis(*result); ++zbin)
      valid_bins.push_back(GetBin(*result, xbin, ybin, zbin));
    original_cov0.ResizeTo(valid_bins.size(), valid_bins.size());
    eigen_cov0.ResizeTo(valid_bins.size(), valid_bins.size());
    for (decltype(valid_bins.size()) i = 0; i < valid_bins.size(); ++i) {
      for (decltype(valid_bins.size()) j = 0; j < valid_bins.size(); ++j) {
        original_cov0(i, j) = original_cov(valid_bins[i], valid_bins[j]);
        original_cov0(j, i) = original_cov(valid_bins[j], valid_bins[i]);
        eigen_cov0(i, j) = eigen_cov(valid_bins[i], valid_bins[j]);
        eigen_cov0(j, i) = eigen_cov(valid_bins[j], valid_bins[i]);
      }
    }
    valid_bins.clear();
    for (auto uzbin = 1; uzbin <= GetNbinsAxis(*uresult); ++uzbin)
      valid_bins.push_back(GetBin(*uresult, uxbin, uybin, uzbin));
    uoriginal_cov0.ResizeTo(valid_bins.size(), valid_bins.size());
    ueigen_cov0.ResizeTo(valid_bins.size(), valid_bins.size());
    for (decltype(valid_bins.size()) i = 0; i < valid_bins.size(); ++i) {
      for (decltype(valid_bins.size()) j = 0; j < valid_bins.size(); ++j) {
        uoriginal_cov0(i, j) = uoriginal_cov(valid_bins[i], valid_bins[j]);
        uoriginal_cov0(j, i) = uoriginal_cov(valid_bins[j], valid_bins[i]);
        ueigen_cov0(i, j) = ueigen_cov(valid_bins[i], valid_bins[j]);
        ueigen_cov0(j, i) = ueigen_cov(valid_bins[j], valid_bins[i]);
      }
    }
  }

  AnalysisData data{
    info_tag,
    uinfo_tag,
    InfoTag,
    uInfoTag,
    *result0,
    *uresult0,
    stat_set,
    ustat_set,
    syst_set,
    usyst_set,
    eigen_set,
    ueigen_set,
    original_cov0,
    uoriginal_cov0,
    eigen_cov0,
    ueigen_cov0
  };
  RunAnalysis(data);
}

}

void Analysis::eigenDecompCalib (TString fileName,
                                 TString unsmoothedFileName,
                                 TString pathName,
                                 int     numberOfNPs,
                                 bool    preserveMergedError)
{
  gROOT->SetBatch(true);

#ifdef __CINT__
  //gSystem->Load("libCalibrationDataInterface.so");
#endif // ifdef __CINT__

  // if ((fileName == "help") || (unsmoothedFileName == "")) {
  //   cout << "usage: eigenDecompCalib(fileName, pathName)" << endl;
  //   cout << "where" << endl;
  //   cout << "            fileName:    name of the smoothedCDI ROOT file" << endl;
  //   cout << "  unsmoothedFileName:    name of the unsmoothed CDI ROOT file" << endl;
  //   cout << "            pathName:    name of the calibration object's path name within the ROOT file (see e.g. output from checkCalibrationFile())" << endl;
  //   cout << "         numberOfNPs:    number of desired smoothed nusiance parameters (+1 for merged NPs)" << endl;
  //   cout << " preserveMergedError:    whether or not to preserve error when merging NPs" << endl;
  //   return;
  // }

  TFile   *f  = TFile::Open(fileName.Data(), "READ");
  TFile   *uf = TFile::Open(unsmoothedFileName.Data(), "READ");
  TFile   *of = NULL;
  TCanvas *c1 = NULL;
  TLatex  *lt = new TLatex();
  TString  tag;
  vector<unsigned long> colors;
  CalibrationDataEigenVariations::IndexSuperSet variations_to_combine, uvariations_to_combine;
  CalibrationDataEigenVariations::IndexSuperSet variations_to_remove, uvariations_to_remove;
  CalibrationDataEigenVariations::IndexSet      variations;

  SetAtlasStyle();

  colors.push_back(kBlue);
  colors.push_back(kRed);
  colors.push_back(kGreen + 2);
  colors.push_back(kOrange + 7);
  colors.push_back(kCyan - 3);
  colors.push_back(kViolet);
  colors.push_back(kOrange + 4);

  lt->SetTextSize(0.04);
  lt->SetTextFont(42);

  Info("eigenDecompCalib", "looking for %s in %s", pathName.Data(), fileName.Data());

  CalibrationDataHistogramContainer *c, *uc;
  f->GetObject(pathName.Data(), c);
  uf->GetObject(pathName.Data(), uc);

  if ((c != NULL) && (uc != NULL)) {
    // TH2F *result = static_cast<TH2F*>((*c)("result")),
    //  *uresult = static_cast<TH2F*>((*uc)("result"));
    // TH2 *result = static_cast<TH2*>((*c)("result")),
    // *uresult    = static_cast<TH2*>((*uc)("result"));
    TH1 *result = static_cast<TH1*>((*c)("result")),
    *uresult    = static_cast<TH1*>((*uc)("result"));
    result->SetDirectory(0);
    uresult->SetDirectory(0);
    CalibrationDataEigenVariations *eigen = new CalibrationDataEigenVariations(c),
    *ueigen                               = new CalibrationDataEigenVariations(uc);
    TMatrixDSym original_cov              = eigen->getEigenCovarianceMatrix(),
                uoriginal_cov             = ueigen->getEigenCovarianceMatrix();
    vector<string> named_variations;

    if ((result == NULL) || (uresult == NULL)) {
      Fatal("eigenDecompCalib", "couldn't find central histogram (named \"result\") - aborting...");
      return;
    }

    Info("eigenDecompCalib", "found %s - starting analysis", pathName.Data());

    std::cout << "initialize" << std::endl;
    eigen->initialize();
    ueigen->initialize();

    // Combination
    if (numberOfNPs > 0) {
      for (int i = eigen->getNumberOfEigenVariations(); i > numberOfNPs; --i) variations.insert(i);
      variations_to_combine.insert(variations);
      variations.clear();
      --numberOfNPs;
    }

    //// Removal
    // for (int i = eigen->getNumberOfEigenVariations(); i > 6; --i)
    //  variations.insert(i);
    // variations_to_remove.insert(variations);
    // variations.clear();

    if (variations_to_remove.size() != 0) tag = "trim";
    else if ((variations_to_combine.size() != 0) && preserveMergedError) tag = "comb_preserve";
    else if ((variations_to_combine.size() != 0) && !preserveMergedError) tag = "comb_cov";
    else tag = "orginal";

    // eigen->removeVariations(variations_to_remove);
    // eigen->mergeVariations(variations_to_combine/*, preserveMergedError*/);
    eigen->mergeVariationsFrom(numberOfNPs);

    if (numberOfNPs <= static_cast<decltype(numberOfNPs)>(ueigen->getNumberOfEigenVariations())) ueigen->mergeVariationsFrom(numberOfNPs);

    // ueigen->removeVariations(uvariations_to_remove);
    // ueigen->mergeVariations(uvariations_to_combine);

    auto eigen_cov  = eigen->getEigenCovarianceMatrixFromVariations(),
         ueigen_cov = ueigen->getEigenCovarianceMatrixFromVariations();

    of = TFile::Open(TString::Format("eigen_plots_%d_var_%s_%s.root",
                                     eigen->getNumberOfEigenVariations(),
                                     tag.Data(),
                                     pathName.Copy().ReplaceAll("/", "-").Data()).Data(),
                     "RECREATE");
    Info("eigenDecompCalib", "plots will be stored in %s", of->GetName());

    c1 = new TCanvas("central", "central", 200, 10, 700, 500);

    {
      TH1 *up, *down;
      if ((eigen->getNumberOfEigenVariations() <= 0) ||
          !eigen->getEigenvectorVariation(0, up, down) ||
          (ueigen->getNumberOfEigenVariations() <= 0) ||
          !ueigen->getEigenvectorVariation(0, up, down)) {
        Fatal("eigenDecompCalib", "no eigen-variaitons - aborting...");
        return;
      }
    }

    Info("eigenDecompCalib", "list of named variations: ");
    for (auto var : eigen->listNamedVariations())
      Info("eigenDecompCalib", "\t%s", var.c_str());

    Info("eigenDecompCalib", "number of eigen-variations:");
    Info("eigenDecompCalib", "\tsmoothed: %d", eigen->getNumberOfEigenVariations());
    Info("eigenDecompCalib", "\tunsmoothed: %d", ueigen->getNumberOfEigenVariations());

    // loop over axes
    auto ndim = result->GetDimension();
    if (ndim == 1) {
    }
    else if (ndim == 2) {
      // project along x-axis
      if (uresult->GetNbinsX() > 1) {
        for (auto ybin = 1; ybin <= uresult->GetNbinsY(); ++ybin) {
          auto ycenter    = uresult->GetYaxis()->GetBinCenter(ybin);
          auto bin        = result->GetYaxis()->FindFixBin(ycenter);
          auto GetAxis = [](const TH1 &h)->const TAxis* {return h.GetYaxis();};
          auto GetNbinsAxis = [](const TH1 &h)->Int_t {return h.GetNbinsX();};
          auto GetBin = [](const TH1 &h, Int_t xbin, Int_t ybin, Int_t zbin)->Int_t {return h.GetBin(ybin, xbin, zbin);};
          auto Project = [](TH2 &h, const char *name, Int_t low, Int_t high, Option_t *o)->TH1D* {
            return h.ProjectionX(name, low, high, o); };
          auto RunAnalysis = [&of, &pathName, &tag, &c1, &colors](AnalysisData &data)->void {
            return DoEigenAnalysis(of, pathName, tag, c1, colors, data); };

          MakeDataFrom2D(bin, result, uresult, c, uc, eigen, ueigen,
              original_cov, uoriginal_cov, eigen_cov, ueigen_cov,
              GetAxis, GetNbinsAxis, GetBin, Project, RunAnalysis);
        }
      }
      // project along y-axis
      if (uresult->GetNbinsY() > 1) {
        for (auto xbin = 1; xbin <= uresult->GetNbinsX(); ++xbin) {
          auto xcenter    = uresult->GetXaxis()->GetBinCenter(xbin);
          auto bin        = result->GetXaxis()->FindFixBin(xcenter);
          auto GetAxis = [](const TH1 &h)->const TAxis* {return h.GetXaxis();};
          auto GetNbinsAxis = [](const TH1 &h)->Int_t {return h.GetNbinsY();};
          auto GetBin = [](const TH1 &h, Int_t xbin, Int_t ybin, Int_t zbin)->Int_t {return h.GetBin(xbin, ybin, zbin);};
          auto Project = [](TH2 &h, const char *name, Int_t low, Int_t high, Option_t *o)->TH1D* {
            return h.ProjectionY(name, low, high, o); };
          auto RunAnalysis = [&of, &pathName, &tag, &c1, &colors](AnalysisData &data)->void {
            return DoEigenAnalysis(of, pathName, tag, c1, colors, data); };

          MakeDataFrom2D(bin, result, uresult, c, uc, eigen, ueigen,
              original_cov, uoriginal_cov, eigen_cov, ueigen_cov,
              GetAxis, GetNbinsAxis, GetBin, Project, RunAnalysis);
        }
      }
    }
    else if (ndim == 3) {
      // project along x-axis
      if (uresult->GetNbinsX() > 1) {
        for (auto ybin = 1; ybin <= uresult->GetNbinsY(); ++ybin) {
          for (auto zbin = 1; zbin <= uresult->GetNbinsZ(); ++zbin) {
            auto ycenter      = uresult->GetYaxis()->GetBinCenter(ybin);
            auto zcenter      = uresult->GetZaxis()->GetBinCenter(zbin);
            auto ybin_center  = result->GetYaxis()->FindFixBin(ycenter);
            auto zbin_center  = result->GetZaxis()->FindFixBin(zcenter);
            auto GetYAxis = [](const TH1 &h)->const TAxis* {return h.GetYaxis();};
            auto GetZAxis = [](const TH1 &h)->const TAxis* {return h.GetZaxis();};
            auto GetNbinsAxis = [](const TH1 &h)->Int_t {return h.GetNbinsX();};
            auto GetBin = [](const TH1 &h, Int_t xbin, Int_t ybin, Int_t zbin)->Int_t {return h.GetBin(zbin, xbin, ybin);};
            auto Project = [](TH3 &h, const char *name, Int_t lowy, Int_t highy, Int_t lowz, Int_t highz)->TH1D* {
              h.GetYaxis()->SetRange(lowy, highy);
              h.GetZaxis()->SetRange(lowz, highz);
              auto proj = h.Project3D("x");
              proj->SetName(name);
              h.GetYaxis()->SetRange(0, 0);
              h.GetZaxis()->SetRange(0, 0);
              return static_cast<TH1D*>(proj); };
            auto RunAnalysis = [&of, &pathName, &tag, &c1, &colors](AnalysisData &data)->void {
              return DoEigenAnalysis(of, pathName, tag, c1, colors, data); };

            MakeDataFrom3D(ybin_center, zbin_center, result, uresult, c, uc, eigen, ueigen,
                original_cov, uoriginal_cov, eigen_cov, ueigen_cov,
                GetYAxis, GetZAxis, GetNbinsAxis, GetBin, Project, RunAnalysis);
          }
        }
      }
      // project along y-axis
      if (uresult->GetNbinsY() > 1) {
        for (auto xbin = 1; xbin <= uresult->GetNbinsX(); ++xbin) {
          for (auto zbin = 1; zbin <= uresult->GetNbinsZ(); ++zbin) {
            auto xcenter      = uresult->GetXaxis()->GetBinCenter(xbin);
            auto zcenter      = uresult->GetZaxis()->GetBinCenter(zbin);
            auto xbin_center  = result->GetXaxis()->FindFixBin(xcenter);
            auto zbin_center  = result->GetZaxis()->FindFixBin(zcenter);
            auto GetXAxis = [](const TH1 &h)->const TAxis* {return h.GetXaxis();};
            auto GetZAxis = [](const TH1 &h)->const TAxis* {return h.GetZaxis();};
            auto GetNbinsAxis = [](const TH1 &h)->Int_t {return h.GetNbinsY();};
            auto GetBin = [](const TH1 &h, Int_t xbin, Int_t ybin, Int_t zbin)->Int_t {return h.GetBin(xbin, zbin, ybin);};
            auto Project = [](TH3 &h, const char *name, Int_t lowx, Int_t highx, Int_t lowz, Int_t highz)->TH1D* {
              h.GetXaxis()->SetRange(lowx, highx);
              h.GetZaxis()->SetRange(lowz, highz);
              auto proj = h.Project3D("y");
              proj->SetName(name);
              h.GetXaxis()->SetRange(0, 0);
              h.GetZaxis()->SetRange(0, 0);
              return static_cast<TH1D*>(proj); };
            auto RunAnalysis = [&of, &pathName, &tag, &c1, &colors](AnalysisData &data)->void {
              return DoEigenAnalysis(of, pathName, tag, c1, colors, data); };

            MakeDataFrom3D(xbin_center, zbin_center, result, uresult, c, uc, eigen, ueigen,
                original_cov, uoriginal_cov, eigen_cov, ueigen_cov,
                GetXAxis, GetZAxis, GetNbinsAxis, GetBin, Project, RunAnalysis);
          }
        }
      }
      // project along z-axis
      if (uresult->GetNbinsZ() > 1) {
        for (auto xbin = 1; xbin <= uresult->GetNbinsX(); ++xbin) {
          for (auto ybin = 1; ybin <= uresult->GetNbinsY(); ++ybin) {
            auto xcenter      = uresult->GetXaxis()->GetBinCenter(xbin);
            auto ycenter      = uresult->GetYaxis()->GetBinCenter(ybin);
            auto xbin_center  = result->GetXaxis()->FindFixBin(xcenter);
            auto ybin_center  = result->GetYaxis()->FindFixBin(ycenter);
            auto GetXAxis = [](const TH1 &h)->const TAxis* {return h.GetXaxis();};
            auto GetYAxis = [](const TH1 &h)->const TAxis* {return h.GetYaxis();};
            auto GetNbinsAxis = [](const TH1 &h)->Int_t {return h.GetNbinsZ();};
            auto GetBin = [](const TH1 &h, Int_t xbin, Int_t ybin, Int_t zbin)->Int_t {return h.GetBin(xbin, ybin, zbin);};
            auto Project = [](TH3 &h, const char *name, Int_t lowx, Int_t highx, Int_t lowy, Int_t highy)->TH1D* {
              h.GetXaxis()->SetRange(lowx, highx);
              h.GetYaxis()->SetRange(lowy, highy);
              auto proj = h.Project3D("z");
              proj->SetName(name);
              h.GetXaxis()->SetRange(0, 0);
              h.GetYaxis()->SetRange(0, 0);
              return static_cast<TH1D*>(proj); };
            auto RunAnalysis = [&of, &pathName, &tag, &c1, &colors](AnalysisData &data)->void {
              return DoEigenAnalysis(of, pathName, tag, c1, colors, data); };

            MakeDataFrom3D(xbin_center, ybin_center, result, uresult, c, uc, eigen, ueigen,
                original_cov, uoriginal_cov, eigen_cov, ueigen_cov,
                GetXAxis, GetYAxis, GetNbinsAxis, GetBin, Project, RunAnalysis);
          }
        }
      }
    }

    delete result;
    delete eigen;
    delete uresult;
    delete ueigen;

    of->Close();
  } else
    Fatal("eigenDecompCalib", "calibration object (%s) cannot be found or is not a histogram container", pathName.Data());


  f->Close();
  uf->Close();
} // eigenDecompCalib

    // // ////////////////////////////////////////////////////////
    // // 2D Plots (experimental)
    // // ////////////////////////////////////////////////////////
    // //result->Draw("E");
    // //uresult->Draw("E SAME");
    // //uresult->Draw("E");
    // auto graph2d = new TGraph2DErrors(result->GetNbinsX()*result->GetNbinsY());
    // auto ugraph2d = new TGraph2DErrors(uresult->GetNbinsX()*uresult->GetNbinsY());
    //
    // auto point_count = 0;
    // for (auto xi = 1; xi <= result->GetNbinsX(); ++xi) {
    //   auto x = result->GetXaxis()->GetBinCenter(xi);
    //   for (auto yi = 1; yi <= result->GetNbinsY(); ++yi) {
    //     auto y = result->GetYaxis()->GetBinCenter(yi);
    //     auto z = result->GetBinContent(xi, yi);
    //     auto z_err = result->GetBinError(xi, yi);
    //     graph2d->SetPoint(point_count, x, y, z);
    //     graph2d->SetPointError(point_count, 0.0, 0.0, z_err);
    //     cout << x << ", " << y << ", " << z << " +/- " << z_err << endl;
    //     ++point_count;
    //   }
    // }
    // point_count = 0;
    // for (auto xi = 1; xi <= uresult->GetNbinsX(); ++xi) {
    //   auto x = uresult->GetXaxis()->GetBinCenter(xi);
    //   for (auto yi = 1; yi <= uresult->GetNbinsY(); ++yi) {
    //     auto y = uresult->GetYaxis()->GetBinCenter(yi);
    //     auto z = uresult->GetBinContent(xi, yi);
    //     auto z_err = uresult->GetBinError(xi, yi);
    //     ugraph2d->SetPoint(point_count, x, y, z);
    //     ugraph2d->SetPointError(point_count, 0.0, 0.0, z_err);
    //     cout << x << ", " << y << ", " << z << " +/- " << z_err << endl;
    //     ++point_count;
    //   }
    // }
    //
    // auto x_width = ugraph2d->GetXaxis()->GetXmax() - ugraph2d->GetXaxis()->GetXmin(),
    //      y_width = ugraph2d->GetYaxis()->GetXmax() - ugraph2d->GetYaxis()->GetXmin(),
    //      z_width = ugraph2d->GetZaxis()->GetXmax() - ugraph2d->GetZaxis()->GetXmin();
    // auto x_min = -0.2*x_width + ugraph2d->GetXaxis()->GetXmin();
    // auto x_max = 0.2*x_width + ugraph2d->GetXaxis()->GetXmax();
    // auto y_min = -0.2*y_width + ugraph2d->GetYaxis()->GetXmin();
    // auto y_max = 0.2*y_width + ugraph2d->GetYaxis()->GetXmax();
    // // auto z_min = ugraph2d->GetZaxis()->GetXmin() < 0 ?
    // // auto z_max = ugraph2d->GetZaxis()->GetXmax() < 0 ?
    //
    // ugraph2d->GetXaxis()->SetCanExtend(kTRUE);
    // ugraph2d->GetYaxis()->SetCanExtend(kTRUE);
    // ugraph2d->GetXaxis()->SetLimits(x_min, x_max);
    // ugraph2d->GetYaxis()->SetLimits(y_min, y_max);
    // graph2d->GetXaxis()->SetLimits(x_min, x_max);
    // graph2d->GetYaxis()->SetLimits(y_min, y_max);
    //
    // ugraph2d->SetMarkerStyle(25);
    // ugraph2d->SetMarkerSize(0.15);
    // ugraph2d->SetMarkerColor(kBlack);
    // graph2d->SetMarkerSize(0.25);
    // graph2d->SetMarkerColor(kBlue);
    //
    // ugraph2d->Draw("P ERR");
    // graph2d->Draw("P SAME");
    //
    // c1->SetLogy();
    // c1->SetName("result2D");
    // of->WriteTObject(c1);
    // c1->Clear();
    //
    // c1->SetLogy(kFALSE);
