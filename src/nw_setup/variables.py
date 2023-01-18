###########
# basic_cond = [ver_cognac, blank_udf, base_udf, core]
ver_cognac = ''
blank_udf = ''
base_udf = ''
core = 1

# nw_cond = [nw_model, strand, n_strand, n_segments, n_cell, n_sc, l_bond, c_n]
nw_model = ''
strand = ''
n_strand = 1
n_segments = 1
strand_type = ''
n_cell = 1
n_sc = 1
l_bond = 1.0
c_n = 1.0

# rnd_cond = [restart, cond_top, histgram_bins]
restart = ''
cond_top = []
histgram_bins = 0

# sim_cond = [entanglement, multi_org, density_org, shrinkage, expand, step_press, press_time, step_rfc, step_rfc_time, equilib_repeat, equilib_time, greenkubo, greenkubo_repeat, greenkubo_time, calc]
entanglement = ''
multi_org = 0
density_org = 0
shrinkage = 0
expand = 1.0
step_press = []
press_time = []
step_rfc = []
step_rfc_time = []
equilib_repeat = 1
equilib_time = []
greenkubo = ''
greenkubo_repeat = 1
greenkubo_time = []
calc = ''

# mod_cond = [n_chains, n_beads_unit, e2e, org_unitcell]
n_chains = 1
n_beads_unit = 1
e2e = 1.0
org_unitcell = 1.0

# sim_cond2 = [multi_mod, density_mod]
multi_mod = 1
density_mod = 1.0

# target_cond = [system, unit_cell, total_net_atom, nu]
system = 1.0
unit_cell = 1.0
total_net_atom = 1
nu = 1.0

target_dir = ''
target_name = ''


# Cognac用の名称設定
nw_name = "Network"
atom_name = ["JP_A", "End_A", "Strand_A", "Side_A", "Solvent"]
bond_name = ["bond_JP-Chn", "bond_Strand", "bond_Side"]
angle_name = ["angle_AAA"]
site_name = ["site_JP", "site_End", "site_Strand", "site_Solvent"]
pair_name = ["site_JP-site_JP", "site_Strand-site_JP", "site_Strand-site_Strand", 
                "site_JP-site_End", "site_Strand-site_End", "site_End-site_End",
                "site_Solvent-site_Solvent", "site_Solvent-site_JP", "site_Solvent-site_End",
                "site_Solvent-site_Strand"]
site_pair_name = [ 
                ["site_JP", "site_JP"], 
                ["site_Strand", "site_JP"], 
                ["site_Strand", "site_Strand"],
                ["site_JP", "site_End"], 
                ["site_Strand", "site_End"], 
                ["site_End", "site_End"],
                ["site_Solvent", "site_Solvent"],
                ["site_Solvent", "site_JP"],
                ["site_Solvent", "site_End"],
                ["site_Solvent", "site_Strand"],
                ]
# Set [R0, K]
harmonic = [0.97, 1000]
# [Potential_Type, theta0, K]		
angle = ['Theta2', 74, 10.0]
# [Cutoff, Scale_1_4_Pair, sigma, epsilon, range]
lj_cond = [2**(1/6), 1.0, 1.0, 1.0, 1.0]	