pt_bins=25
# pt_bins=100
bandwidth="0.4"
# bandwidth="30"
order=1
smoothing_suffix="_order_"$order"_smoothing_"$bandwidth"_ptbins_"$pt_bins
nB_eigens=-1
nC_eigens=-1
nLight_eigens=-1

# AnalyzeEigenVariation CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI-July12_v1.root CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI-July12_BEFORE_v1.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/C/DStar_extrap_pre_ttbar_PDF_7b_rebin_SF -1
# AnalyzeEigenVariation CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI-July12_v1.root CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI-July12_BEFORE_v1.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/C/default_SF -1
#
# AnalyzeEigenVariation CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI-July12_v1.root CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI-July12_BEFORE_v1.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/Light/negative_tags_extrap_SF -1
# AnalyzeEigenVariation CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI-July12_v1.root CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI-July12_BEFORE_v1.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/Light/default_SF -1

# 13TeV
# CDI_basename="CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI"
CDI_basename="CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI_order_1_smoothing_0.4_ptbins_100"

# # continuous calibrations
AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/Continuous/B/default_SF $nB_eigens
AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/Continuous/B/default_SF $nB_eigens
AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/Continuous/C/default_SF $nC_eigens
AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/Continuous/Light/default_SF $nLight_eigens

# # B
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/B/default_SF $nB_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/B/default_SF $nB_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/B/default_SF $nB_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/B/default_SF $nB_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/B/default_SF $nB_eigens

# # C
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/C/default_SF $nC_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/C/default_SF $nC_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/C/default_SF $nC_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/C/default_SF $nC_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/C/default_SF $nC_eigens
#
# # Light
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/Light/default_SF $nLight_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/Light/default_SF $nLight_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/Light/default_SF $nLight_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/Light/default_SF $nLight_eigens
# AnalyzeEigenVariation ${CDI_basename}${smoothing_suffix}.root ${CDI_basename}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/Light/default_SF $nLight_eigens
