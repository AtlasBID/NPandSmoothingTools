pt_bins=25
bandwidth="0.4"
order=1
smoothing_suffix="_smooth"
# nB_eigens=-1
nB_eigens=15
# nC_eigens=-1
nC_eigens=5
# nLight_eigens=-1
nLight_eigens=30

echo "First arg: $1"
arg="$(echo $1)"
# arg="$(echo $1 | cut -d "." -f 1)"
# echo "First arg with processing: $arg"

# continuous calibrations
b_bandwidth="0.5"
c_bandwidth="0.6"
light_bandwidth="0.7"
smoothing_suffix="_order_"$order"_smoothing_"$b_bandwidth"_"$c_bandwidth"_"$light_bandwidth"_ptbins_"$pt_bins
# CDI_basename="CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI_order_1_smoothing_0.4_ptbins_100"

AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/Continuous/B/IntegratedWP_TopoJets_SF $nB_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/Continuous/C/IntegratedWP_TopoJets_SF $nC_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/Continuous/Light/IntegratedWP_TopoJets_SF $nLight_eigens

# # B
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/B/default_SF $nB_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/B/default_SF $nB_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/B/default_SF $nB_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/B/default_SF $nB_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/B/default_SF $nB_eigens
#
# # C
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/C/default_SF $nC_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/C/default_SF $nC_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/C/default_SF $nC_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/C/default_SF $nC_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/C/default_SF $nC_eigens
#
# # Light
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/Light/default_SF $nLight_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/Light/default_SF $nLight_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/Light/default_SF $nLight_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/Light/default_SF $nLight_eigens
# AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/Light/default_SF $nLight_eigens
