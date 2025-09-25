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

load analysis/data/foldfusion_output/Results/P00383/AlphaFold/AF-P00383-F1-model_v4_processed.pdb, P00383_AF
load analysis/data/foldfusion_output/Results/P00383/JamdaScorer/2P4T/NAP_A_157.sdf, P00383_NAP_refined
load analysis/data/foldfusion_output/Results/P00383/LigandExtractor/2P4T/NAP_A_157.sdf, P00383_NAP_donor
remove hydro

hide everything
show cartoon, P00383_AF
color 0x8fb3ff, P00383_AF
show sticks, P00383_NAP_refined
color orange, P00383_NAP_refined
show sticks, P00383_NAP_donor
color 0x2ca02c, P00383_NAP_donor
set stick_transparency, 0.35, P00383_NAP_donor

align P00383_NAP_donor, P00383_NAP_refined

select P00383_contact_shell, P00383_AF within 4 of P00383_NAP_refined
show lines, P00383_contact_shell
color grey70, P00383_contact_shell
set line_width, 2, P00383_contact_shell

orient P00383_NAP_refined
zoom P00383_NAP_refined, 12

# Rotate manually and run `ray 2000,1500` + `png ...` as desired.
