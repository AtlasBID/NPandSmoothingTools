pt_bins=100
#bandwidth="2.25"
#bandwidth="1.75"
# bandwidth="0.45"
bandwidth="0.4"
# bandwidth="30"
order=1
#smooth_b_jets="false"
smooth_b_jets="true"
b_NP_set="5,3,2"
#smooth_c_jets="false"
smooth_c_jets="true"
c_NP_set="4,4,3"
#smooth_light_jets="false"
smooth_light_jets="true"
light_NP_set="12,5,4"

# 13TeV
#SmoothCDI CDIFiles/13TeV/2015-PreRecomm-13TeV-MC12-CDI-July10-v1.root $pt_bins $bandwidth $order $smooth_b_jets $smooth_c_jets $smooth_light_jets
#SmoothCDI CDIFiles/13TeV/2015-PreRecomm-13TeV-MC12-CDI_August3-v1.root $pt_bins $bandwidth $order $smooth_b_jets $smooth_c_jets $smooth_light_jets
#SmoothCDI CDIFiles/13TeV/2015-PreRecomm-13TeV-MC12-CDI.root $pt_bins $bandwidth $order $smooth_b_jets $smooth_c_jets $smooth_light_jets
#SmoothCDI CDIFiles/13TeV/2015-PreRecomm-13TeV-MC12-CDI-2.root $pt_bins $bandwidth $order $smooth_b_jets $smooth_c_jets $smooth_light_jets
#SmoothCDI CDIFiles/13TeV/2016-Winter-13TeV-MC15-CDI-March4_unofficial.root ".*_SF" $pt_bins $bandwidth $order $smooth_b_jets $b_NP_set $smooth_c_jets $c_NP_set $smooth_light_jets $light_NP_set
#SmoothCDI CDIFiles/13TeV/Continuous_240316.root $pt_bins $bandwidth $order $smooth_b_jets $smooth_c_jets $smooth_light_jets
SmoothCDI CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI.root ".*_SF" $pt_bins $bandwidth $order $smooth_b_jets $b_NP_set $smooth_c_jets $c_NP_set $smooth_light_jets $light_NP_set

# causes (harmless) seg-fault on file close with recent versions of CDI package
## 8TeV
#SmoothCDI CDIFiles/8TeV/Cumulative/2014-Winter-8TeV-MC12-CDI.root $pt_bins $bandwidth $order $smooth_b_jets $smooth_c_jets $smooth_light_jets
