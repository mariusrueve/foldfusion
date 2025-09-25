reinitialize everything
set ray_opaque_background, on
bg_color white
set ray_shadows, 0
set ambient, 0.5
set ray_trace_mode, 1
set ray_trace_gain, 0.1
set cartoon_transparency, 0.2
set stick_radius, 0.12
set sphere_scale, 0.3

load analysis/data/foldfusion_output/Results/P80882/AlphaFold/AF-P80882-F1-model_v4_processed.pdb, P80882_AF
load analysis/data/foldfusion_output/Results/P80882/JamdaScorer/1HLQ/SF4_A_76.sdf, P80882_SF4
remove hydro

hide everything
show cartoon, P80882_AF
color 0xc0d6ff, P80882_AF
show surface, P80882_AF
set transparency, 0.4, P80882_AF
show sticks, P80882_SF4
color 0xd62728, P80882_SF4

select P80882_clash_zone, P80882_AF within 4 of P80882_SF4
color 0xff6666, P80882_clash_zone
set surface_color, 0xff6666, P80882_clash_zone

# Visualize clash distances
remove P80882_overlap
distance P80882_overlap, P80882_SF4, P80882_clash_zone
set dash_width, 3, P80882_overlap
set dash_color, red, P80882_overlap
label P80882_overlap, "%s%s: %.2f Å" % (resn, resi, length)
set label_color, black, P80882_overlap
set label_bg_color, white, P80882_overlap
set label_size, 18, P80882_overlap

orient P80882_SF4
zoom P80882_SF4, 14

# Rotate manually and run `ray 2000,1500` + `png ...` as desired.
