import logging
import os
import pandas as pd
import ROOT

def main():
    logging.basicConfig(#filename='NP_reduction.log',
                        level=logging.DEBUG)
    rootCore_import_result = ROOT.gROOT.Macro('$ROOTCOREDIR/scripts/load_packages.C')
    if rootCore_import_result != 0 and rootCore_import_result != 1:
        logging.critical("Couldn't import RootCore package libraries. Aborting...")
        exit(1)
    else:
        from ROOT import Analysis, TH1D#, TMatrixDSym
        from ROOT.Analysis import CalibrationDataHistogramContainer, CalibrationDataEigenVariations

    # NOTE: integrated cases
    input_file_name = "CDIFiles/13TeV/2016-Winter-13TeV-MC15-CDI-March4_unofficial.root"
    # input_file_name = "CDIFiles/13TeV/2016-Winter-13TeV-MC15-CDI-March4_unofficial_order_1_smoothing_0.4_ptbins_100.root"
    # NOTE: Light
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_60/Light/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_70/Light/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_77/Light/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_85/Light/default_SF"
    # NOTE: C
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_60/C/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_70/C/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_77/C/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_85/C/default_SF"
    # NOTE: B
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_60/B/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_70/B/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_77/B/default_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/FixedCutBEff_85/B/default_SF"

    # NOTE: continuous cases
    # input_file_name = "CDIFiles/13TeV/Continuous_240316.root"
    # input_file_name = "CDIFiles/13TeV/Continuous_240316_order_1_smoothing_0.4_ptbins_25.root"
    # NOTE: Light
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/continuous/Light/IntegratedWP_5TagBins_SF"
    # NOTE: C
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/continuous/C/IntegratedWP_5TagBins_SF"
    # NOTE: B
    CDI_object_path = "MV2c20/AntiKt4EMTopoJets/continuous/B/IntegratedWP_5TagBins_SF"
    # CDI_object_path = "MV2c20/AntiKt4EMTopoJets/continuous/B/IntegratedWP_7TagBins_SF"

    n_EV_components = 3
    input_file_name_only = input_file_name.split(os.path.sep)[-1].replace(".root", "")
    CDI_object_path_name = CDI_object_path.replace(os.path.sep, "_")
    data = {}
    indices = []
    header = []
    input_file = ROOT.TFile.Open(input_file_name, "READ")
    c = input_file.Get(CDI_object_path)
    print type(c)
    eigen = CalibrationDataEigenVariations(c)
    eigen.initialize()
    large_original_cov = eigen.getEigenCovarianceMatrix()
    jac = eigen.getJacobianReductionMatrix()
    original_cov = large_original_cov.Similarity(jac)
    eigen.mergeVariationsFrom(n_EV_components - 1)

    result = c("result")
    result_list = []
    result_error_list_up = []
    result_error_list_down = []

    nbins = (result.GetNbinsX() + 2)*(result.GetNbinsY() + 2 if result.GetDimension() > 1 else 1)*(result.GetNbinsZ() + 2 if result.GetDimension() > 2  else 1)

    for ibin in range(0, nbins):
        if result.IsBinOverflow(ibin) or result.IsBinUnderflow(ibin): continue
        xbin, ybin, zbin = ROOT.Long(0), ROOT.Long(0), ROOT.Long(0)
        result.GetBinXYZ(ibin, xbin, ybin, zbin)

        xbin_center = str(result.GetXaxis().GetBinCenter(xbin))
        xbin_title = result.GetXaxis().GetTitle()
        bin_title = xbin_title + "=" + xbin_center

        if result.GetDimension() > 1:
            ybin_center = str(result.GetYaxis().GetBinCenter(ybin))
            ybin_title = result.GetYaxis().GetTitle()
            bin_title = bin_title + ":" + ybin_title + "=" + ybin_center

        if result.GetDimension() > 2:
            zbin_center = str(result.GetZaxis().GetBinCenter(zbin))
            zbin_title = result.GetZaxis().GetTitle()
            bin_title = bin_title + ":" + zbin_title + "=" + zbin_center

        header.append(bin_title)
        up = result.GetBinContent(ibin)
        down = result.GetBinError(ibin)
        result_list.append(up)
        result_error_list_up.append(down)
        result_error_list_down.append(-down)

    indices.append("result")
    indices.append("result_stat_error__up")
    indices.append("result_stat_error__down")
    data[indices[-3]] = result_list
    data[indices[-2]] = result_error_list_up
    data[indices[-1]] = result_error_list_down

    eigen_vec_data_up = []
    eigen_vec_data_down = []
    for ieigen in range(eigen.getNumberOfEigenVariations()):
        up = ROOT.TH1D()
        down = ROOT.TH1D()
        eigen.getEigenvectorVariation(ieigen, up, down)
        eigen_vec_up = []
        eigen_vec_down = []
        for ibin in range(0, nbins):
            if result.IsBinOverflow(ibin) or result.IsBinUnderflow(ibin): continue
            # xbin, ybin, zbin = ROOT.Long(0), ROOT.Long(0), ROOT.Long(0)
            # result.GetBinXYZ(ibin, xbin, ybin, zbin)
            eigen_vec_up.append(up.GetBinContent(ibin) - result.GetBinContent(ibin))
            eigen_vec_down.append(down.GetBinContent(ibin) - result.GetBinContent(ibin))
        eigen_vec_data_up.append(eigen_vec_up)
        eigen_vec_data_down.append(eigen_vec_down)


    for syst in c.listUncertainties():
        hist = c(syst)
        if (syst == "result") or (syst == "extrapolation") or \
           (syst == "comment") or (syst == "combined") or \
           (syst == "MCreference") or (syst == "MChadronisation") or \
           (syst == "ReducedSets") or (syst == "statistics") or \
           (syst == "systematics") or not hist.InheritsFrom("TH1"):
            continue
        indices.append(hist.GetName().replace(" ", "_") + "__up")
        indices.append(hist.GetName().replace(" ", "_") + "__down")
        up_list = []
        down_list = []
        for ibin in range(0, nbins):
            if result.IsBinOverflow(ibin) or result.IsBinUnderflow(ibin): continue
            xbin, ybin, zbin = ROOT.Long(0), ROOT.Long(0), ROOT.Long(0)
            result.GetBinXYZ(ibin, xbin, ybin, zbin)
            up = hist.GetBinContent(ibin)# + result.GetBinContent(ixbin, iybin, izbin)
            down = hist.GetBinError(ibin)# + result.GetBinContent(ixbin, iybin, izbin)
            up_list.append(up)
            down_list.append(down)
        data[indices[-2]] = up_list
        data[indices[-1]] = down_list

    df = pd.DataFrame(data, index=header).transpose()
    df.to_csv(input_file_name_only + "_" + CDI_object_path_name + "_n" + str(n_EV_components) + "_systs.csv")

    cov_mat = []
    for irow in range(original_cov.GetNrows()):
        cov_mat.append([])
        for icol in range(original_cov.GetNrows()):
            cov_mat[-1].append(original_cov(irow, icol))
    df = pd.DataFrame(data=cov_mat, index=header, columns=header)
    df.to_csv(input_file_name_only + "_" + CDI_object_path_name + "_n" + str(n_EV_components) + "_matrix.csv")

    df = pd.DataFrame(data=eigen_vec_data_up, columns=header)
    df.to_csv(input_file_name_only + "_" + CDI_object_path_name + "_n" + str(n_EV_components) + "_eigen_vec_up.csv")

    df = pd.DataFrame(data=eigen_vec_data_down, columns=header)
    df.to_csv(input_file_name_only + "_" + CDI_object_path_name + "_n" + str(n_EV_components) + "_eigen_vec_down.csv")

    logging.debug("Ignore any seg-faults past this point. This is due to something in the CalibrationDataEigenVariations.")


if __name__ == "__main__":
    main()
