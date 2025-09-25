reinitialize everything
set ray_opaque_background, on
bg_color white
set ray_shadows, 0
set ambient, 0.5
set ray_trace_mode, 1
set ray_trace_gain, 0.1
set cartoon_transparency, 0.2
set stick_radius, 0.12
set sphere_scale, 0.4

load analysis/data/foldfusion_output/Results/P08306/AlphaFold/AF-P08306-F1-model_v4_processed.pdb, P08306_AF
load analysis/data/foldfusion_output/Results/P08306/JamdaScorer/1AR1/CU_B_270.sdf, P08306_CU
remove hydro

hide everything
show cartoon, P08306_AF
color 0xf2f2f2, P08306_AF
show sticks, P08306_CU
color 0xffa600, P08306_CU
show spheres, P08306_CU
set sphere_scale, 0.4, P08306_CU

select P08306_shell, (P08306_AF within 3.2 of P08306_CU) and (name ND1+NE2+OD1+OD2+OE1+OE2+SG+SD)
show sticks, P08306_shell
color 0x1f77b4, P08306_shell

remove P08306_contacts
distance P08306_contacts, P08306_CU, P08306_shell
set dash_width, 3, P08306_contacts
set dash_color, grey50, P08306_contacts
label P08306_contacts, "%s%s: %.2f Å" % (resn, resi, length)
set label_color, black, P08306_contacts
set label_bg_color, white, P08306_contacts
set label_size, 18, P08306_contacts

orient P08306_CU
zoom P08306_CU, 12

# Rotate manually and run `ray 2000,1500` + `png ...` as desired.
