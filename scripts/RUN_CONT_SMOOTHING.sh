pt_bins=25
#bandwidth="2.25"
#bandwidth="1.75"
b_bandwidth="0.5"
c_bandwidth="0.6"
light_bandwidth="0.7"
order=1
# smooth_b_jets="false"
smooth_b_jets="true"
b_NP_set="30,27,15"
# smooth_c_jets="false"
smooth_c_jets="true"
c_NP_set="20,15,10"
# smooth_light_jets="false"
smooth_light_jets="true"
light_NP_set="60,50,30"

# 13TeV

#SmoothContCDI CDIFiles/13TeV/Continuous_240316.root "IntegratedWP_.*_SF" $pt_bins $bandwidth $order $smooth_b_jets $b_NP_set $smooth_c_jets $c_NP_set $smooth_light_jets $light_NP_set
SmoothContCDI CDIFiles/13TeV/2016-20_7-13TeV-MC15-CDI_order_1_smoothing_0.4_ptbins_100.root "IntegratedWP_.*_SF" $pt_bins $b_bandwidth $c_bandwidth $light_bandwidth $order $smooth_b_jets $b_NP_set $smooth_c_jets $c_NP_set $smooth_light_jets $light_NP_set
