#include <algorithm>
#include <string>
#include <sstream>
#include <vector>
#include <map>
#include <iostream>
#include <fstream>

#include <TFile.h>
#include <TH1.h>

#include "CalibrationDataInterface/CalibrationDataContainer.h"
#include "CalibrationDataInterface/CalibrationDataEigenVariations.h"

std::vector<std::string> &split(const std::string &s, char delim, std::vector<std::string> &elems)
{
  std::stringstream ss(s);
  std::string item;
  while (std::getline(ss, item, delim)) {
    elems.push_back(item);
  }
  return elems;
}

std::vector<std::string> split(const std::string &s, char delim)
{
  std::vector<std::string> elems;
  split(s, delim, elems);
  return elems;
}

void trim_number_string (std::string &num)
{
  using std::string;
  if (num.find('.') != string::npos) {
    auto decimal_position = num.find_last_not_of('0') + 1;
    if (std::all_of(num.begin() + decimal_position, num.end(), [](char &c) {return (c == '.' || c == '0');})) {
      num.replace(decimal_position, string::npos, "");
    }
    if (num.back() == '.') num.erase(num.size() - 1);
  }
}

void to_csv (std::ostream &ostream, const std::vector<std::vector<double>>&, const std::vector<std::string> &header, const std::vector<std::string> &indices = {});
void to_csv (std::ostream &ostream, const std::map<std::string, std::vector<double>>&, const std::vector<std::string> &header);

int main (int argc, const char *argv[])
{
  using std::string;
  using std::vector;
  using std::map;

  vector<string> args;

  for (int i = 0; i < argc; ++i) args.emplace_back(argv[i]);

  if (args.size() < 3) {
    Error("ConvertCDIToCSV", "This program needs the location of the input CDI file as it's first parameter and the path to the calibration object as it's second parameter. Aborting!");
    return 1;
  }
  string input_file_name = args[1];
  string CDI_object_path = args[2];
  int n_EV_components = -1;
  if (args.size() > 3) n_EV_components = std::stod(args[3]);
  if (n_EV_components < -1) n_EV_components = -1;

  // NOTE: integrated cases
  // string input_file_name = "CDIFiles/13TeV/2016-Winter-13TeV-MC15-CDI-March4_unofficial.root";
  // string input_file_name = "CDIFiles/13TeV/2016-Winter-13TeV-MC15-CDI-March4_unofficial_order_1_smoothing_0.4_ptbins_100.root";
  // NOTE: Light
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_60/Light/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_70/Light/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_77/Light/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_85/Light/default_SF";
  // NOTE: C
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_60/C/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_70/C/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_77/C/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_85/C/default_SF";
  // NOTE: B
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_60/B/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_70/B/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_77/B/default_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_85/B/default_SF";

  // NOTE: continuous cases
  // string input_file_name = "CDIFiles/13TeV/Continuous_240316.root";
  // string input_file_name = "CDIFiles/13TeV/Continuous_240316_order_1_smoothing_0.4_ptbins_25.root";
  // NOTE: Light
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/continuous/Light/IntegratedWP_5TagBins_SF";
  // NOTE: C
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/continuous/C/IntegratedWP_5TagBins_SF";
  // NOTE: B
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/continuous/B/IntegratedWP_5TagBins_SF";
  // string CDI_object_path = "MV2c20/AntiKt4EMTopoJets/continuous/B/IntegratedWP_7TagBins_SF";

  // int n_EV_components = 5;

  string input_file_name_only;
  {
    auto nix_spit_file_name = split(input_file_name, '/');
    nix_spit_file_name.back().replace(nix_spit_file_name.back().rfind(".root"), string::npos, "");
    input_file_name_only = nix_spit_file_name.back();
  }
  string CDI_object_path_name = CDI_object_path;
  std::replace_if(CDI_object_path_name.begin(), CDI_object_path_name.end(), [](const char &c)->bool {return c == '/';}, '_');

  map<string, vector<double>> data;
  vector<string> indices;
  vector<string> header;
  auto input_file = TFile::Open(input_file_name.c_str(), "READ");
  Analysis::CalibrationDataHistogramContainer *c;
  input_file->GetObject(CDI_object_path.c_str(), c);
  if (c) {
    auto result = static_cast<TH1*>((*c)("result"));
    result->SetDirectory(0);
    auto eigen = new Analysis::CalibrationDataEigenVariations(c);
    eigen->initialize();
    auto large_original_cov = eigen->getEigenCovarianceMatrix();
    auto jac = eigen->getJacobianReductionMatrix();
    auto original_cov = large_original_cov.Similarity(jac);
    if (n_EV_components > 0) eigen->mergeVariationsFrom(n_EV_components - 1);
    auto eigen_cov  = eigen->getEigenCovarianceMatrixFromVariations();
    auto nbins = (result->GetNbinsX() + 2)*(result->GetDimension() > 1 ? result->GetNbinsY() + 2 : 1)*(result->GetDimension() > 2 ? result->GetNbinsZ() + 2 : 1);

    // NOTE: central values
    vector<double> result_list;
    vector<double> result_error_list_up;
    vector<double> result_error_list_down;

    for (Int_t ibin = 0; ibin < nbins; ++ibin) {
      if (result->IsBinOverflow(ibin) or result->IsBinUnderflow(ibin)) continue;
      Int_t xbin, ybin, zbin;
      result->GetBinXYZ(ibin, xbin, ybin, zbin);

      auto xbin_center = std::to_string(result->GetXaxis()->GetBinCenter(xbin));
      trim_number_string(xbin_center);
      string xbin_title = result->GetXaxis()->GetTitle();
      auto bin_title = xbin_title + "=" + xbin_center;
      if (result->GetDimension() > 1) {
        auto ybin_center = std::to_string(result->GetYaxis()->GetBinCenter(ybin));
        trim_number_string(ybin_center);
        string ybin_title = result->GetYaxis()->GetTitle();
        bin_title += ":" + ybin_title + "=" + ybin_center;
      }
      if (result->GetDimension() > 2) {
        auto zbin_center = std::to_string(result->GetZaxis()->GetBinCenter(zbin));
        trim_number_string(zbin_center);
        string zbin_title = result->GetZaxis()->GetTitle();
        bin_title += ":" + zbin_title + "=" + zbin_center;
      }

      header.push_back(bin_title);
      auto up = result->GetBinContent(ibin);
      auto down = result->GetBinError(ibin);
      // std::cout << bin_title << ": " << up << " +/- " << down << std::endl;
      result_list.push_back(up);
      result_error_list_up.push_back(down);
      result_error_list_down.push_back(-down);
    }

    indices.push_back("result");
    indices.push_back("result_stat_error__up");
    indices.push_back("result_stat_error__down");
    data["result"] = result_list;
    data["result_stat_error__up"] = result_error_list_up;
    data["result_stat_error__down"] = result_error_list_down;

    // NOTE: EV components
    vector<vector<double>> eigen_vec_data_up;
    vector<vector<double>> eigen_vec_data_down;
    for (unsigned ieigen = 0; ieigen < eigen->getNumberOfEigenVariations(); ++ieigen) {
      TH1 *up = nullptr, *down = nullptr;
      if (eigen->getEigenvectorVariation(ieigen, up, down)) {
        vector<double> eigen_vec_up;
        vector<double> eigen_vec_down;
        for (Int_t ibin = 0; ibin < nbins; ++ibin) {
          if (result->IsBinOverflow(ibin) || result->IsBinUnderflow(ibin)) continue;
          eigen_vec_up.push_back(up->GetBinContent(ibin) - result->GetBinContent(ibin));
          eigen_vec_down.push_back(down->GetBinContent(ibin) - result->GetBinContent(ibin));
        }
        eigen_vec_data_up.push_back(eigen_vec_up);
        eigen_vec_data_down.push_back(eigen_vec_down);
      }
    }

    // NOTE: covariance matrix
    vector<vector<double>> cov_mat;
    for (Int_t irow = 0; irow < original_cov.GetNrows(); ++irow) {
      cov_mat.push_back({});
      for (Int_t icol = 0; icol < original_cov.GetNrows(); ++icol) {
        cov_mat.back().push_back(original_cov(irow, icol));
      }
    }

    std::ofstream csvfile;

    csvfile.open (input_file_name_only + "_" + CDI_object_path_name + "_n" + std::to_string(n_EV_components) + "_matrix.csv");
    to_csv(csvfile, cov_mat, header, header);
    csvfile.close();

    csvfile.open (input_file_name_only + "_" + CDI_object_path_name + "_n" + std::to_string(n_EV_components) + "_eigen_vec_up.csv");
    to_csv(csvfile, eigen_vec_data_up, header);
    csvfile.close();

    csvfile.open (input_file_name_only + "_" + CDI_object_path_name + "_n" + std::to_string(n_EV_components) + "_eigen_vec_down.csv");
    to_csv(csvfile, eigen_vec_data_down, header);
    csvfile.close();

    csvfile.open (input_file_name_only + "_" + CDI_object_path_name + "_n" + std::to_string(n_EV_components) + "_systs.csv");
    to_csv(csvfile, data, header);
    csvfile.close();

    delete eigen;
    delete result;
  }

  return 0;
}

void to_csv (std::ostream &ostream, const std::vector<std::vector<double>> &data, const std::vector<std::string> &header, const std::vector<std::string> &indices)
{
  // NOTE: header
  for (auto &h : header)
    ostream << ',' << h;
  ostream << std::endl;

  // NOTE: actual data
  auto do_indices = indices.size() == data.size();
  unsigned data_counter = 0;
  for (auto &d : data) {
    if (do_indices)
      ostream << indices[data_counter];
    else
      ostream << data_counter;
    for (auto &dd : d)
      ostream << ',' << dd;
    ostream << std::endl;
    ++data_counter;
  }
}

void to_csv (std::ostream &ostream, const std::map<std::string, std::vector<double>> &data, const std::vector<std::string> &header)
{
  // NOTE: header
  for (auto &h : header)
    ostream << ',' << h;
  ostream << std::endl;

  // NOTE: actual data
  for (auto &p : data) {
    ostream << p.first;
    for (auto &dd : p.second)
      ostream << ',' << dd;
    ostream << std::endl;
  }
}
