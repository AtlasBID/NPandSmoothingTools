pt_bins=25
bandwidth="0.4"
order=1
smoothing_suffix="_smooth"
nB_eigens=15
nC_eigens=5
nLight_eigens=30

echo "First arg: $1"
arg="$(echo $1 | cut -d "." -f 1)"
echo "First arg with processing: $arg"

# continuous calibrations
b_bandwidth="0.5"
c_bandwidth="0.6"
light_bandwidth="0.7"
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/Continuous/B/default_SF $nB_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/Continuous/C/default_SF $nC_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/Continuous/Light/default_SF $nLight_eigens
#AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/Continuous/B/default_SF $nB_eigens
#AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/Continuous/C/default_SF $nC_eigens
#AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/Continuous/Light/default_SF $nLight_eigens

# B
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/B/default_SF $nB_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/B/default_SF $nB_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/B/default_SF $nB_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/B/default_SF $nB_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/B/default_SF $nB_eigens

# C
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/C/default_SF $nC_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/C/default_SF $nC_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/C/default_SF $nC_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/C/default_SF $nC_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/C/default_SF $nC_eigens

# Light
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_60/Light/default_SF $nLight_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_70/Light/default_SF $nLight_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_77/Light/default_SF $nLight_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt4EMTopoJets/FixedCutBEff_85/Light/default_SF $nLight_eigens
AnalyzeEigenVariation ${arg}${smoothing_suffix}.root ${arg}.root MV2c10/AntiKt2PV0TrackJets/FixedCutBEff_70/Light/default_SF $nLight_eigens
