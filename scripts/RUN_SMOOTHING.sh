pt_bins=100
bandwidth="0.4"
order=1
smooth_b_jets="true"
b_NP_set="5,3,2"
smooth_c_jets="true"
c_NP_set="4,4,3"
smooth_light_jets="true"
light_NP_set="12,5,4"

echo "First arg: $1"
echo "Second arg (order): $2"
echo "Third arg (smoothing): $3"
echo "Fourth arg (ptbins): $4"

# SmoothCDI $1 ".*_SF" $4 $bandwidth $2 $smooth_b_jets $b_NP_set $smooth_c_jets $c_NP_set $smooth_light_jets $light_NP_set
SmoothCDI $1 ".*_SF" $4 $3 $2 $smooth_b_jets $b_NP_set $smooth_c_jets $c_NP_set $smooth_light_jets $light_NP_set
