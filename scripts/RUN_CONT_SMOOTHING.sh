pt_bins=25
#bandwidth="2.25"
#bandwidth="1.75"
bandwidth="0.4"
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

#SmoothContCDI CDIFiles/13TeV/Continuous_240316.root "IntegratedWP_.*_SF" $pt_bins $bandwidth $order $smooth_b_jets $b_NP_set $smooth_c_jets $c_NP_set $smooth_light_jets $light_NP_set
SmoothContCDI CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI_order_1_smoothing_0.4_ptbins_100.root "IntegratedWP_.*_SF" $pt_bins $bandwidth $order $smooth_b_jets $b_NP_set $smooth_c_jets $c_NP_set $smooth_light_jets $light_NP_set
