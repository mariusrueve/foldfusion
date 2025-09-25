reinitialize
# ---------- INPUTS ----------
# AlphaFold model
load analysis/data/foldfusion_output/Results/A4D0H5/AlphaFold/AF-A4D0H5-F1-model_v4_processed.pdb, target
# SIENA-aligned donor protein
load analysis/data/foldfusion_output/Results/A4D0H5/Siena/ensemble/2WET_5.pdb, donor
remove donor and (chain B+C)
# Ligands (pre/post JAMDAscorer)
load analysis/data/foldfusion_output/Results/A4D0H5/LigandExtractor/2WET/FAD_A_1513.sdf, lig_pre
load analysis/data/foldfusion_output/Results/A4D0H5/JamdaScorer/2WET/FAD_A_1513.sdf, lig_post

# ---------- GLOBAL STYLE ----------
bg_color white
set ray_opaque_background, off
set antialias, 2
set ambient, 0.2
set spec_reflect, 0
set spec_power, 80
set stick_radius, 0.25
set cartoon_sampling, 14
set ray_shadow, off
hide everything, all

# Target protein
show cartoon, target
color gray80, target
set cartoon_transparency, 0.25, target

# Donor protein (from SIENA)
show cartoon, donor and polymer
color marine, donor and polymer
set cartoon_transparency, 0.25, donor

# Ligands
show sticks, lig_pre or lig_post
util.cbag, lig_pre
util.cbag, lig_post
color firebrick, lig_pre
color forest, lig_post

# ---------- POCKET SELECTION FROM RESIDUE LIST ----------
select pocket_sel, chain A and resi 8+9+10+11+12+13+34+35+36+39+40+41+42+44+45+46+47+193+194+195+225+226+227+228+229+231+253+281+283+294+325+327+343+344+345+349+352+355+356+357+358+361
# Visualize pocket as its own surface/mesh
create pocket, target and pocket_sel
as surface, pocket
color orange, pocket
# pocket transparency (surface + mesh)
set transparency, 0.7, pocket    # 0=opaque, 1=invisible (surface)
set mesh_alpha, 0.4, pocket      # 0=opaque, 1=invisible (mesh)

as mesh, pocket
set mesh_width, 0.5
set mesh_cutoff, 9999

# ---------- CLEANUP ----------
# Keep view consistent. You rotate once, then trigger panels below and screenshot.
orient target

# ---------- QUICK PANEL TOGGLES (type one command then screenshot) ----------
# 2) AlphaFold structure only
alias panel_af_only, disable all; enable target; show cartoon, target; zoom target, 0
# 3) Pocket on AF structure
alias panel_pocket, disable all; enable target pocket; show cartoon, target; show surface, pocket; show mesh, pocket; zoom pocket, 12
# 4) SIENA donor overlay
alias panel_siena, disable all; enable target donor pocket; show cartoon, target or donor; show mesh, pocket; zoom pocket, 12
# 5) Ligand before JAMDAscorer
alias panel_pre, disable all; enable target lig_pre pocket; show cartoon, target; show sticks, lig_pre; show mesh, pocket; zoom lig_pre, 12
# 6) Ligand after JAMDAscorer
alias panel_post, disable all; enable target lig_post pocket; show cartoon, target; show sticks, lig_post; show mesh, pocket; zoom lig_post, 12

# ---------- OPTIONAL: OUTLINE AND CLASH HINTS ----------
# Outline protein for print clarity
set cartoon_side_chain_helper, on
set surface_quality, 1
set two_sided_lighting, on
set depth_cue, 0

# Uncomment to compute simple contacts around lig_post (visual cue, not a metric)
# find_pairs (target and name C+N+O+S), lig_post, cutoff=3.6

# ---------- NOTES ----------
# Use: run script, rotate to your preferred angle once.
# Then execute: panel_af_only  -> screenshot
#                panel_pocket  -> screenshot
#                panel_siena   -> screenshot
#                panel_pre     -> screenshot
#                panel_post    -> screenshot
# For ray-traced export: add 'ray 2000,1500' before 'png file.png, dpi=300'.
