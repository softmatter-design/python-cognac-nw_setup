#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
import numpy as np
import os
from UDFManager import UDFManager
import nw_setup.variables as var
##########################################
# Initial_UDF の作成
##########################################
################################################################################
# UDFファイルを設定し、バッチ処理を作成
def setup_baseudf(calcd_data_dic):
	# 初期udfの内容を作成する
	make_base()
	# すべてのアトムの位置座標及びボンド情報を設定
	setup_atoms(calcd_data_dic)
	return
################################################################################
def make_base():
	#--- create an empty UDF file ---
	target_udf = os.path.join(var.target_dir, var.base_udf)
	with open(target_udf, 'w') as f:
		f.write(r'\include{"%s"}' % var.blank_udf)

	u = UDFManager(target_udf)
	# goto global data
	u.jump(-1)

	#--- Simulation_Conditions ---
	# Solver
	p = 'Simulation_Conditions.Solver.'
	u.put('Dynamics', p + 'Solver_Type')
	if var.entanglement == "NO_Entangled":
		u.put('NPT_Andersen_Kremer_Grest', 	p + 'Dynamics.Dynamics_Algorithm')
		u.put(var.total_net_atom, 				p + 'Dynamics.NPT_Andersen_Kremer_Grest.Cell_Mass')
		u.put(0.5, 							p + 'Dynamics.NPT_Andersen_Kremer_Grest.Friction')
	elif var.entanglement == "Entangled":
		u.put('NVT_Kremer_Grest', 			p + 'Dynamics.Dynamics_Algorithm')
		u.put(0.5, 							p + 'Dynamics.NVT_Kremer_Grest.Friction')
	# Boundary_Conditions
	p = 'Simulation_Conditions.Boundary_Conditions'
	u.put(['PERIODIC', 'PERIODIC', 'PERIODIC', 1], p)
	#
	p = "Simulation_Conditions.Dynamics_Conditions.Moment."
	u.put(10000, p + "Interval_of_Calc_Moment")
	u.put(1, p + "Calc_Moment")
	u.put(1, p + "Stop_Translation")
	u.put(1, p + "Stop_Rotation")

	# Calc_Potential_Flags
	p = 'Simulation_Conditions.Calc_Potential_Flags.'
	u.put(1, p + 'Bond')
	u.put(0, p + 'Angle')
	u.put(1, p + 'Non_Bonding_Interchain')
	u.put(1, p + 'Non_Bonding_1_3')
	u.put(1, p + 'Non_Bonding_1_4')
	u.put(1, p + 'Non_Bonding_Intrachain')

	# Output_Flags.Statistics
	p = 'Simulation_Conditions.Output_Flags.Statistics.'
	u.put(1, p + 'Energy')
	u.put(1, p + 'Temperature')
	u.put(1, p + 'Pressure')
	u.put(1, p + 'Stress')
	u.put(1, p + 'Volume')
	u.put(1, p + 'Density')
	u.put(1, p + 'Cell')
	u.put(0, p + 'Wall_Pressure')
	u.put(0, p + 'Energy_Flow')

	# Output_Flags.Structure
	p = 'Simulation_Conditions.Output_Flags.Structure.'
	u.put(1, p + 'Position')
	u.put(0, p + 'Velocity')
	u.put(0, p + 'Force')

	#--- Initial_Structure ---
	# Initial_Unit_Cell
	p = 'Initial_Structure.Initial_Unit_Cell.'
	if var.entanglement == "NO_Entangled":
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0, p + 'Density')
		u.put([var.system*var.expand, var.system*var.expand, var.system*var.expand, 90.0, 90.0, 90.0], p + 'Cell_Size')
	else:
		u.put(var.density_mod, p + 'Density')
		u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
	
	#--- Molecular_Attributes ---
	# Atomes
	for i, atomname in enumerate(var.atom_name):
		p = 'Molecular_Attributes.Atom_Type[].'
		u.put(atomname, 	p + 'Name', [i])
		u.put(1.0, 			p + 'Mass', [i])
	# Bond
	for i, bondname in enumerate(var.bond_name):
		p = 'Molecular_Attributes.Bond_Potential[].'
		u.put(bondname, 			p + 'Name', [i])
		u.put('Harmonic', 			p + 'Potential_Type', [i])
		u.put(var.harmonic[0], 	p + 'R0', [i])
		u.put(var.harmonic[1], 	p + 'Harmonic.K', [i])
	# Angle
	for i, anglename in enumerate(var.angle_name):
		p = 'Molecular_Attributes.Angle_Potential[].'
		u.put(anglename, 			p + 'Name', [i])
		u.put(var.angle[0], 	p + 'Potential_Type', [i])
		u.put(var.angle[1], 	p + 'theta0', [i])
		u.put(var.angle[2], 	p + 'Theta2.K', [i])

	# Site
	for i, sitename in enumerate(var.site_name):
		p = 'Molecular_Attributes.Interaction_Site_Type[].'
		u.put(sitename, 		p + 'Name', [i])
		u.put(1, 				p + 'Num_of_Atoms', [i])
		u.put(var.lj_cond[4], 	p + 'Range', [i])

	#--- Pair_Interaction[] ---
	for i, pairname in enumerate(var.pair_name):
		p = 'Interactions.Pair_Interaction[].'
		u.put(pairname,   					p + 'Name', [i])
		u.put('Lennard_Jones', 				p + 'Potential_Type', [i])
		u.put(var.site_pair_name[i][0],	p + 'Site1_Name', [i])
		u.put(var.site_pair_name[i][1],	p + 'Site2_Name', [i])
		u.put(var.lj_cond[0],				p + 'Cutoff', [i])
		u.put(var.lj_cond[1],				p + 'Scale_1_4_Pair', [i])
		u.put(var.lj_cond[2],				p + 'Lennard_Jones.sigma', [i])
		u.put(var.lj_cond[3],				p + 'Lennard_Jones.epsilon', [i])

	#--- Write UDF ---
	u.write(target_udf)
	return

################################################################################
# すべてのアトムの位置座標及びボンド情報を設定
def setup_atoms(calcd_data_dic):
	# 多重配置の場合の位置シフト量
	shift = 2
	#
	target_udf = os.path.join(var.target_dir, var.base_udf)
	u = UDFManager(target_udf)
	u.jump(-1)

	#--- Set_of_Molecules の入力
	p = 'Set_of_Molecules.molecule[].'
	pa = p + 'atom[].'
	pi = p + 'interaction_Site[].'
	pb = p + 'bond[].'
	pang = p +  'angle[].'
	#
	count = 0
	for mul in range(len(calcd_data_dic)):
		atom_all = calcd_data_dic[mul]["atom_all"]
		bond_all = calcd_data_dic[mul]["bond_all"]
		pos_all = calcd_data_dic[mul]["pos_all"]
		angle_all = calcd_data_dic[mul]["angle_all"]
		#
		u.put(var.nw_name + '_' + str(mul), p + 'Mol_Name', [count])
		# beads
		n_atom = 0
		for atom in atom_all:
			# atom
			id_shift = len(atom_all) # + var.n_solvent
			atom_id = atom[0] + count*id_shift
			u.put(atom_id, 						pa + 'Atom_ID', [count, n_atom])
			u.put(var.atom_name[atom[1]], 		pa + 'Atom_Name', [count, n_atom])
			u.put(var.atom_name[atom[1]], 		pa + 'Atom_Type_Name', [count, n_atom])
			u.put(0, 							pa + 'Chirality', [count, n_atom])
			u.put(1, 							pa + 'Main_Chain', [count, n_atom])
			# interaction site
			u.put(var.site_name[atom[2]], 		pi + 'Type_Name', [count, n_atom])
			u.put(n_atom, 						pi + 'atom[]', [count, n_atom, 0])
			n_atom += 1
		# bonds
		n_bond = 0
		ang_list=[]
		for bond in range(len(bond_all)):
			u.put(var.bond_name[bond_all[bond][0]], 	pb + 'Potential_Name', [count, n_bond])
			u.put(bond_all[bond][1][0], 				pb + 'atom1', [count, n_bond])
			ang_list.append(bond_all[bond][1][0])
			u.put(bond_all[bond][1][1], 				pb + 'atom2', [count, n_bond])
			n_bond += 1

		# angles
		n_ang = 0
		for a_list in angle_all:
			for j in range(len(a_list)-2):
				u.put(var.angle_name[0], 	pang + 'Potential_Name', [count, n_ang])
				u.put(a_list[j], 			pang + 'atom1', [count, n_ang])
				u.put(a_list[j + 1], 		pang + 'atom2', [count, n_ang])
				u.put(a_list[j + 2], 		pang + 'atom3', [count, n_ang])
				n_ang += 1

		# Draw_Attributes
		color = ["Red", "Green", "Blue", "Magenta", "Cyan", "Yellow", "White", "Black", "Gray"]
		mm = mul % 9
		u.put([var.nw_name + '_' + str(mul), color[mm], 1.0, 1.0], 'Draw_Attributes.Molecule[]', [count])
		count += 1

	# アトムの座標位置を、シフトしながら、設定
	sp = 'Structure.Position.mol[].atom[]'
	count = 0
	totalatom = 0
	# tmp_set = []
	# ネットワークのセグメントをセット
	for mul in range(len(calcd_data_dic)):
		shift_vec = var.expand*count*shift*np.array(np.random.rand(3))
		pos_all = calcd_data_dic[mul]["pos_all"]
		for i in range(len(pos_all)):
			if var.entanglement == "NO_Entangled":
				mod_pos = var.expand*var.unit_cell*np.array(list(pos_all[i])) + var.expand*shift_vec
			else:
				mod_pos = var.unit_cell*np.array(list(pos_all[i])) + shift_vec
			# print(mod_pos)
			u.put(list(mod_pos), sp, [count, i])
			# tmp_set.append(list(mod_pos))
			totalatom +=1
		count+=1

	#--- Write UDF ---
	u.write(target_udf)

	return
