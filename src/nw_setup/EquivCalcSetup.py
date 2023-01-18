#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
# import numpy as np
import platform
import os
from UDFManager import UDFManager
import nw_setup.variables as var
##########################################
# UDF の作成
####################################
# UDFファイルを設定し、バッチ処理を作成
def setup_udf():
	if platform.system() == "Windows":
		batch = ""
	elif platform.system() == "Linux":
		batch = "#!/bin/bash\n"

	###############
	# entanglementに応じて計算条件を選択
	if var.entanglement == "Entangled":
		batch = entangle_calc(batch)
	elif var.entanglement == "NO_Entangled":
		batch = npt_calc(batch)

	#####################
	# バッチファイルを作成
	f_batch = os.path.join(var.target_dir, '_Calc_all.bat')
	with open(f_batch,'w') as f:
		# f.write(batch_all)
		f.write(batch)
		if platform.system() == "Linux":
			os.chmod(f_batch, 0o777)
	return

#####################################################################################
# ターミナルのタイトルを設定
def make_title(batch, title):
	if platform.system() == "Windows":
		batch += "title " + title + "\n"
	elif platform.system() == "Linux":
		batch += r'echo -ne "\033]0; ' + title + ' \007"' + '\n'
	return batch

# ファイル名の処理
def make_step(fn_ext, batch, f_eval):
	present_udf = fn_ext[0] + fn_ext[1]
	out_udf = present_udf.replace("uin", "out")
	batch += var.ver_cognac + ' -I ' + present_udf + ' -O ' + out_udf + ' -n ' + str(var.core) + ' \n'
	if f_eval:
		batch += 'evaluate_nw' + ' ' + out_udf + '\n'
	read_udf = out_udf
	return present_udf, read_udf, batch

######################
# 各種バッチ条件を設定
######################
######################################################################
# KG鎖の計算
def entangle_calc(batch):
	if var.strand_type == 'KG':
		# Force Capped LJ によりステップワイズに初期化
		# r = 1.1
		batch = make_title(batch, var.target_name + "Calculating-Init")
		fn_ext = ['Init_', "uin.udf"]
		time = [0.01, 100000, 1000]
		f_eval = 1
		present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
		step_nonbond_setup(var.base_udf, '', present_udf, time, 1.1)
		pre = read_udf
		template = present_udf
		#
		for r in var.step_rfc:
			# 平衡化
			batch = make_title(batch, var.target_name + "Calculating-Pre_" + str(round(r, 3)).replace('.', '_'))
			fn_ext = ['Pre_rfc_' + str(round(r, 3)).replace('.', '_') + "_", "uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
			step_nonbond_setup(template, pre, present_udf, var.step_rfc_time, r)
			pre = read_udf
			template = present_udf
		# KG 鎖に設定
		time = [0.01, 1000000, 10000]
		batch = make_title(batch, var.target_name + "Calculating-KG")
		fn_ext = ['KG_', "uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
		kg_setup(template, pre, present_udf, time)
		pre = read_udf
		template = present_udf
		# 平衡化計算
		for i in range(var.equilib_repeat):
			# 平衡化
			batch = make_title(batch, var.target_name + "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
			eq_setup(template, pre, present_udf, var.equilib_time)
			pre = read_udf
			template = present_udf
		# グリーン久保
		if var.greenkubo_repeat != 0:
			for i in range(var.greenkubo_repeat):
				# 平衡化
				batch = make_title(batch, var.target_name + "Calculating-GK_" + str(i))
				fn_ext = ['GK_' + str(i) + "_", "uin.udf"]
				f_eval = 1
				present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
				greenkubo_setup(template, pre, present_udf, var.greenkubo_time)
				pre = read_udf
				template = present_udf
		return batch
	elif var.strand_type == 'FENE':
		# Force Capped LJ によりステップワイズに初期化
		# r = 1.1
		batch = make_title(batch, var.target_name + "Calculating-Init")
		fn_ext = ['Init_', "uin.udf"]
		time = [0.01, 100000, 1000]
		f_eval = 1
		present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
		step_nonbond_setup(var.base_udf, '', present_udf, time, r=1.1)
		pre = read_udf
		template = present_udf
		#
		for r in var.step_rfc:
			# 平衡化
			batch = make_title(batch, var.target_name + "Calculating-Pre_" + str(round(r, 3)).replace('.', '_'))
			fn_ext = ['Pre_rfc_' + str(round(r, 3)).replace('.', '_') + "_", "uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
			step_nonbond_setup(template, pre, present_udf, var.step_rfc_time, r)
			pre = read_udf
			template = present_udf
		# KG 鎖に設定
		time = [0.01, 1000000, 10000]
		batch = make_title(batch, var.target_name + "Calculating-FENE")
		fn_ext = ['FENE_', "uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
		FENE_setup(template, pre, present_udf, time)
		pre = read_udf
		template = present_udf
		# 平衡化計算
		for i in range(var.equilib_repeat):
			# 平衡化
			batch = make_title(batch, var.target_name + "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
			eq_setup(template, pre, present_udf, var.equilib_time)
			pre = read_udf
			template = present_udf
		# グリーン久保
		if var.greenkubo_repeat != 0:
			for i in range(var.greenkubo_repeat):
				# 平衡化
				batch = make_title(batch, var.target_name + "Calculating-GK_" + str(i))
				fn_ext = ['GK_' + str(i) + "_", "uin.udf"]
				f_eval = 1
				present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
				greenkubo_setup(template, pre, present_udf, var.greenkubo_time)
				pre = read_udf
				template = present_udf
		return batch

	elif var.strand_type == 'Harmonic':
		# 初期状態で弱いセグメント間相互作用を入れて緩和
		batch = make_title(batch, var.target_name + "Calculating-Init")
		fn_ext = ['Init_', "uin.udf"]
		time = [0.01, 100000, 1000]
		f_eval = 1
		present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
		harmonic_setup(var.base_udf, '', present_udf, time, 1.1)
		pre = read_udf
		template = present_udf
		# セグメント間相互作用をなくす
		time = [0.01, 1000000, 10000]
		batch = make_title(batch, var.target_name + "Calculating-Harmonic")
		fn_ext = ['Harmonic_', "uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
		harmonic(template, pre, present_udf, time)
		pre = read_udf
		template = present_udf
		# 平衡化計算
		for i in range(var.equilib_repeat):
			# 平衡化
			batch = make_title(batch, var.target_name + "Calculating-Eq_" + str(i))
			fn_ext = ['Eq_' + str(i) + "_", "uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
			eq_setup(template, pre, present_udf, var.equilib_time)
			pre = read_udf
			template = present_udf
		# グリーン久保
		if var.greenkubo_repeat != 0:
			for i in range(var.greenkubo_repeat):
				# 平衡化
				batch = make_title(batch, var.target_name + "Calculating-GK_" + str(i))
				fn_ext = ['GK_' + str(i) + "_", "uin.udf"]
				f_eval = 1
				present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
				greenkubo_setup(template, pre, present_udf, var.greenkubo_time)
				pre = read_udf
				template = present_udf
		return batch

###########################################
# NPT 条件で、設定密度まで圧縮
def npt_calc(batch):
	# NPTの設定
	pres = var.step_press[0]
	batch = make_title(batch, var.target_name + "Calculating-Ini_pres_" + str(pres).replace('.', '_'))
	fn_ext = ['Init_pres_' + str(pres).replace('.', '_') + '_', "uin.udf"]
	time = [0.001, 1000, 100]
	f_eval = 0
	present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
	npt_setup(var.base_udf, '', present_udf, time, pres)
	pre = read_udf
	template = present_udf
	# ステップワイズに圧力増加
	for pres in var.step_press[1:]:
		batch = make_title(batch, var.target_name + "Calculating-Compress_" + str(pres).replace('.', '_'))
		fn_ext = ['Compress_pres_' + str(pres).replace('.', '_') + '_', "uin.udf"]
		time = var.press_time 
		f_eval = 1
		present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
		npt_setup(template, pre, present_udf, time, pres)
		pre = read_udf
		template = present_udf
	# KG 鎖への遷移
	time = [0.01, 1000, 100]
	batch = make_title(batch, var.target_name + "Calculating-pre_KG")
	fn_ext = ['PreKG_', "uin.udf"]
	f_eval = 1
	present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
	pre_kg_setup(template, pre, present_udf, time)
	pre = read_udf
	template = present_udf
	# KG 鎖に設定
	time = [0.01, 100000, 1000]
	batch = make_title(batch, var.target_name + "Calculating-KG")
	fn_ext = ['SetupKG_', "uin.udf"]
	f_eval = 1
	present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
	kg_setup(template, pre, present_udf, time)
	pre = read_udf
	template = present_udf
	# 平衡化計算
	for i in range(var.equilib_repeat):
		# 平衡化
		batch = make_title(batch, var.target_name + "Calculating-Eq_" + str(i))
		fn_ext = ['Eq_' + str(i) + "_", "uin.udf"]
		f_eval = 1
		present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
		eq_setup(template, pre, present_udf, var.equilib_time)
		pre = read_udf
		template = present_udf
	# グリーン久保
	if var.greenkubo_repeat != 0:
		for i in range(var.greenkubo_repeat):
			# 平衡化
			batch = make_title(batch, var.target_name + "Calculating-GK_" + str(i))
			fn_ext = ['GK_' + str(i) + "_", "uin.udf"]
			f_eval = 1
			present_udf, read_udf, batch = make_step(fn_ext, batch, f_eval)
			greenkubo_setup(template, pre, present_udf, var.greenkubo_time)
			pre = read_udf
			template = present_udf

	return batch

###############
# UDF 作成条件
###############
##############################################################################
# ユーザーポテンシャルにより、Force Capped LJ で、ステップワイズにノンボンドを増加
def step_nonbond_setup(template, read_udf, present_udf, time, r_fc):
	u = UDFManager(os.path.join(var.target_dir, template))
	# goto global data
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(100000000., 	p + 'Max_Force')
	u.put(time[0], 		p + 'Time.delta_T')
	u.put(time[1], 		p + 'Time.Total_Steps')
	u.put(time[2], 		p + 'Time.Output_Interval_Steps')
	u.put(1.0, 			p + 'Temperature.Temperature')
	u.put(0., 			p + 'Pressure_Stress.Pressure')
	# Calc_Potential_Flags
	p = 'Simulation_Conditions.Calc_Potential_Flags.'
	u.put(1, p + 'Bond')
	u.put(1, p + 'Angle')
	u.put(1, p + 'Non_Bonding_Interchain')
	u.put(1, p + 'Non_Bonding_1_3')
	u.put(1, p + 'Non_Bonding_1_4')
	u.put(1, p + 'Non_Bonding_Intrachain')
	#--- Initial_Structure ---
	# Initial_Unit_Cell
	p = 'Initial_Structure.Initial_Unit_Cell.'
	u.put(var.density_mod, p + 'Density')
	u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
	# Generate_Method
	if read_udf == '':
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', p+'Method')
		u.put(['', -1, 0, 0], p+'Restart')
		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		u.put(1, p + 'Relaxation')
		u.put('DYNAMICS', p + 'Method')
		u.put(100, p + 'Max_Relax_Force')
		u.put(100000, p + 'Max_Relax_Steps')
	else:
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', p+'Method')
		u.put([read_udf, -1, 1, 0], p+'Restart')
		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		u.put(0, p + 'Relaxation')
	#--- Simulation_Conditions ---
	# bond
	for i, b_name in enumerate(var.bond_name):
		p = 'Molecular_Attributes.Bond_Potential[].'
		u.put(b_name, 		p + 'Name', [i])
		if read_udf == '':
			u.put('Harmonic', 	p + 'Potential_Type', [i])
			u.put(0.97,			p + 'R0', [i])
			u.put(1000, 		p + 'Harmonic.K', [i])
		else:
			u.put('Bond_Polynomial', 	p + 'Potential_Type', [i])
			u.put(0.9609,		p + 'R0', [i])
			u.put(3,			p + 'Bond_Polynomial.N', [i])
			u.put(20.2026,	p + 'Bond_Polynomial.p[]', [i, 0])
			u.put(490.628,	p + 'Bond_Polynomial.p[]', [i, 1])
			u.put(2256.76,	p + 'Bond_Polynomial.p[]', [i, 2])
			u.put(9685.31,	p + 'Bond_Polynomial.p[]', [i, 3])
			u.put('YES',			p + 'Bond_Polynomial.Use_Equilibrated_Value', [i])
	# Angle
	for i, anglename in enumerate(var.angle_name):
		p = 'Molecular_Attributes.Angle_Potential[].'
		u.put(anglename, 		p + 'Name', [i])
		u.put('Force_Cap_LJ', 	p + 'Potential_Type', [i])
		u.put(73.0, 				p + 'theta0', [i])
		u.put(1.0, 			p + 'Force_Cap_LJ.sigma', [i])
		u.put(1.0, 			p + 'Force_Cap_LJ.epsilon', [i])
		u.put(1.122462, 		p + 'Force_Cap_LJ.cutoff', [i])
		u.put(0.8, 			p + 'Force_Cap_LJ.r_FC', [i])
	#--- Pair_Interaction[] ---
	for i, pairname in enumerate(var.pair_name):
		p = 'Interactions.Pair_Interaction[].'
		u.put(pairname,   					p + 'Name', [i])
		u.put('Force_Cap_LJ', 		p + 'Potential_Type', [i])
		u.put(var.site_pair_name[i][0],	p + 'Site1_Name', [i])
		u.put(var.site_pair_name[i][1],	p + 'Site2_Name', [i])
		u.put(1.12246204830937,				p + 'Cutoff', [i])
		u.put(1.0,							p + 'Scale_1_4_Pair', [i])
		u.put(1.0,							p + 'Force_Cap_LJ.sigma', [i])
		u.put(1.0,							p + 'Force_Cap_LJ.epsilon', [i])
		u.put(r_fc,							p + 'Force_Cap_LJ.r_FC', [i])
	#--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return

###############################################
# ボンドをHarmonic、ノンボンドをLJとしてKG鎖を設定
def pre_kg_setup(template, read_udf, present_udf, time):
	u = UDFManager(os.path.join(var.target_dir, template))
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(time[0],  p+'Time.delta_T')
	u.put(time[1],  p+'Time.Total_Steps')
	u.put(time[2],  p+'Time.Output_Interval_Steps')
	u.put(1.0, 		p + 'Temperature.Temperature')
	u.put(0., 		p + 'Pressure_Stress.Pressure')
	# Solver
	p = 'Simulation_Conditions.Solver.'
	u.put('Dynamics', 			p + 'Solver_Type')
	u.put('NVT_Kremer_Grest', 	p + 'Dynamics.Dynamics_Algorithm')
	u.put(0.5, 					p + 'Dynamics.NVT_Kremer_Grest.Friction')
	# Calc_Potential_Flags
	p = 'Simulation_Conditions.Calc_Potential_Flags.'
	u.put(1, p + 'Bond')
	u.put(0, p + 'Angle')
	u.put(1, p + 'Non_Bonding_Interchain')
	u.put(1, p + 'Non_Bonding_1_3')
	u.put(1, p + 'Non_Bonding_1_4')
	u.put(1, p + 'Non_Bonding_Intrachain')
	#--- Initial_Structure ---
	# Initial_Unit_Cell
	p = 'Initial_Structure.Initial_Unit_Cell.'
	u.put(var.density_mod , p + 'Density')
	u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', p+'Method')
	u.put([read_udf, -1, 0, 0], p+'Restart')
	p = 'Initial_Structure.Relaxation.'
	u.put(1, p + 'Relaxation')
	u.put(100.0, p + 'Max_Relax_Force')
	#--- Simulation_Conditions ---
	# Bond
	for i, bondname in enumerate(var.bond_name):
		p = 'Molecular_Attributes.Bond_Potential[].'
		u.put(bondname, 	p + 'Name', [i])
		u.put('Harmonic', 	p + 'Potential_Type', [i])
		u.put(0.97, 		p + 'R0', [i])
		u.put(1000, 		p + 'Harmonic.K', [i])
	# Site
	for i, sitename in enumerate(var.site_name):
		p = 'Molecular_Attributes.Interaction_Site_Type[].'
		u.put(sitename, 	p + 'Name', [i])
		u.put(1, 			p + 'Num_of_Atoms', [i])
		u.put(1.0, 			p + 'Range', [i])
	#--- Pair_Interaction[] ---
	for i, pairname in enumerate(var.pair_name):
		p = 'Interactions.Pair_Interaction[].'
		u.put(pairname,   					p + 'Name', [i])
		u.put('Lennard_Jones', 				p + 'Potential_Type', [i])
		u.put(var.site_pair_name[i][0],	p + 'Site1_Name', [i])
		u.put(var.site_pair_name[i][1],	p + 'Site2_Name', [i])
		u.put(2**(1/6),						p + 'Cutoff', [i])
		u.put(1.0,							p + 'Scale_1_4_Pair', [i])
		u.put(1.0,							p + 'Lennard_Jones.sigma', [i])
		u.put(1.0,							p + 'Lennard_Jones.epsilon', [i])
	#--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return


###############################################
# ボンドをFENE、ノンボンドをLJとしてKG鎖を設定
def kg_setup(template, read_udf, present_udf, time):
	u = UDFManager(os.path.join(var.target_dir, template))
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(time[0],  p+'Time.delta_T')
	u.put(time[1],  p+'Time.Total_Steps')
	u.put(time[2],  p+'Time.Output_Interval_Steps')
	u.put(1.0, 		p + 'Temperature.Temperature')
	u.put(0., 		p + 'Pressure_Stress.Pressure')
	# Solver
	p = 'Simulation_Conditions.Solver.'
	u.put('Dynamics', 			p + 'Solver_Type')
	u.put('NVT_Kremer_Grest', 	p + 'Dynamics.Dynamics_Algorithm')
	u.put(0.5, 					p + 'Dynamics.NVT_Kremer_Grest.Friction')
	# Calc_Potential_Flags
	p = 'Simulation_Conditions.Calc_Potential_Flags.'
	u.put(1, p + 'Bond')
	u.put(0, p + 'Angle')
	u.put(1, p + 'Non_Bonding_Interchain')
	u.put(1, p + 'Non_Bonding_1_3')
	u.put(1, p + 'Non_Bonding_1_4')
	u.put(1, p + 'Non_Bonding_Intrachain')
	#--- Initial_Structure ---
	# Initial_Unit_Cell
	p = 'Initial_Structure.Initial_Unit_Cell.'
	u.put(var.density_mod , p + 'Density')
	u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', p+'Method')
	u.put([read_udf, -1, 0, 0], p+'Restart')
	p = 'Initial_Structure.Relaxation.'
	u.put(1, p + 'Relaxation')
	#--- Simulation_Conditions ---
	# Bond
	for i, b_name in enumerate(var.bond_name):
		p = 'Molecular_Attributes.Bond_Potential[].'
		u.put(b_name, 		p + 'Name', [i])
		u.put('FENE_LJ', 	p + 'Potential_Type', [i])
		u.put(1.0,			p + 'R0', [i])
		u.put(1.5,			p + 'FENE_LJ.R_max', [i])
		u.put(30,			p + 'FENE_LJ.K', [i])
		u.put(1.0,			p + 'FENE_LJ.sigma', [i])
		u.put(1.0,			p + 'FENE_LJ.epsilon', [i])
	# Site
	for i, sitename in enumerate(var.site_name):
		p = 'Molecular_Attributes.Interaction_Site_Type[].'
		u.put(sitename, 	p + 'Name', [i])
		u.put(1, 			p + 'Num_of_Atoms', [i])
		u.put(1.0, 			p + 'Range', [i])
	#--- Pair_Interaction[] ---
	for i, pairname in enumerate(var.pair_name):
		p = 'Interactions.Pair_Interaction[].'
		u.put(pairname,   					p + 'Name', [i])
		u.put('Lennard_Jones', 				p + 'Potential_Type', [i])
		u.put(var.site_pair_name[i][0],	p + 'Site1_Name', [i])
		u.put(var.site_pair_name[i][1],	p + 'Site2_Name', [i])
		u.put(2**(1/6),						p + 'Cutoff', [i])
		u.put(1.0,							p + 'Scale_1_4_Pair', [i])
		u.put(1.0,							p + 'Lennard_Jones.sigma', [i])
		u.put(1.0,							p + 'Lennard_Jones.epsilon', [i])
	#--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return

###############################################
# ボンドをFENE、ノンボンドをLJとしてKG鎖を設定
def FENE_setup(template, read_udf, present_udf, time):
	u = UDFManager(os.path.join(var.target_dir, template))
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(time[0],  p+'Time.delta_T')
	u.put(time[1],  p+'Time.Total_Steps')
	u.put(time[2],  p+'Time.Output_Interval_Steps')
	u.put(1.0, 		p + 'Temperature.Temperature')
	u.put(0., 		p + 'Pressure_Stress.Pressure')
	# Calc_Potential_Flags
	p = 'Simulation_Conditions.Calc_Potential_Flags.'
	u.put(1, p + 'Bond')
	u.put(1, p + 'Angle')
	u.put(0, p + 'Non_Bonding_Interchain')
	u.put(0, p + 'Non_Bonding_1_3')
	u.put(0, p + 'Non_Bonding_1_4')
	u.put(0, p + 'Non_Bonding_Intrachain')
	#--- Initial_Structure ---
	# Initial_Unit_Cell
	p = 'Initial_Structure.Initial_Unit_Cell.'
	u.put(var.density_mod , p + 'Density')
	u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', p+'Method')
	u.put([read_udf, -1, 0, 0], p+'Restart')
	p = 'Initial_Structure.Relaxation.'
	u.put(1, p + 'Relaxation')
	#--- Simulation_Conditions ---
	# Bond
	for i, b_name in enumerate(var.bond_name):
		p = 'Molecular_Attributes.Bond_Potential[].'
		u.put(b_name, 		p + 'Name', [i])
		u.put('FENE_LJ', 	p + 'Potential_Type', [i])
		u.put(1.0,			p + 'R0', [i])
		u.put(1.5,			p + 'FENE_LJ.R_max', [i])
		u.put(30,			p + 'FENE_LJ.K', [i])
		u.put(1.0,			p + 'FENE_LJ.sigma', [i])
		u.put(1.0,			p + 'FENE_LJ.epsilon', [i])
	# Angle
	for i, anglename in enumerate(var.angle_name):
		p = 'Molecular_Attributes.Angle_Potential[].'
		u.put(anglename, 		p + 'Name', [i])
		u.put('Force_Cap_LJ', 	p + 'Potential_Type', [i])
		u.put(73.0, 				p + 'theta0', [i])
		u.put(1.0, 			p + 'Force_Cap_LJ.sigma', [i])
		u.put(1.0, 			p + 'Force_Cap_LJ.epsilon', [i])
		u.put(1.122462, 		p + 'Force_Cap_LJ.cutoff', [i])
		u.put(0.8, 			p + 'Force_Cap_LJ.r_FC', [i])
	#--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return

##############################################################################
# Force Capped LJ で弱い相互作用を入れて初期緩和
def harmonic_setup(template, read_udf, present_udf, time, r_fc):
	u = UDFManager(os.path.join(var.target_dir, template))
	# goto global data
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(100000000., 	p + 'Max_Force')
	u.put(time[0], 		p + 'Time.delta_T')
	u.put(time[1], 		p + 'Time.Total_Steps')
	u.put(time[2], 		p + 'Time.Output_Interval_Steps')
	u.put(1.0, 			p + 'Temperature.Temperature')
	u.put(0., 			p + 'Pressure_Stress.Pressure')
	# Calc_Potential_Flags
	p = 'Simulation_Conditions.Calc_Potential_Flags.'
	u.put(1, p + 'Bond')
	u.put(1, p + 'Angle')
	u.put(1, p + 'Non_Bonding_Interchain')
	u.put(1, p + 'Non_Bonding_1_3')
	u.put(1, p + 'Non_Bonding_1_4')
	u.put(1, p + 'Non_Bonding_Intrachain')
	#--- Initial_Structure ---
	# Initial_Unit_Cell
	p = 'Initial_Structure.Initial_Unit_Cell.'
	u.put(var.density_mod, p + 'Density')
	u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', p+'Method')
	u.put(['', -1, 0, 0], p+'Restart')
	# Relaxation
	p = 'Initial_Structure.Relaxation.'
	u.put(1, p + 'Relaxation')
	u.put('DYNAMICS', p + 'Method')
	u.put(100, p + 'Max_Relax_Force')
	u.put(100000, p + 'Max_Relax_Steps')
	#--- Simulation_Conditions ---
	# bond
	for i, b_name in enumerate(var.bond_name):
		p = 'Molecular_Attributes.Bond_Potential[].'
		u.put(b_name, 		p + 'Name', [i])
		u.put('Harmonic', 	p + 'Potential_Type', [i])
		u.put(0.97,			p + 'R0', [i])
		u.put(1000, 		p + 'Harmonic.K', [i])
	# Angle
	for i, anglename in enumerate(var.angle_name):
		p = 'Molecular_Attributes.Angle_Potential[].'
		u.put(anglename, 		p + 'Name', [i])
		u.put('Force_Cap_LJ', 	p + 'Potential_Type', [i])
		u.put(73.0, 				p + 'theta0', [i])
		u.put(1.0, 			p + 'Force_Cap_LJ.sigma', [i])
		u.put(1.0, 			p + 'Force_Cap_LJ.epsilon', [i])
		u.put(1.122462, 		p + 'Force_Cap_LJ.cutoff', [i])
		u.put(0.8, 			p + 'Force_Cap_LJ.r_FC', [i])
	#--- Pair_Interaction[] ---
	for i, pairname in enumerate(var.pair_name):
		p = 'Interactions.Pair_Interaction[].'
		u.put(pairname,   					p + 'Name', [i])
		u.put('Force_Cap_LJ', 		p + 'Potential_Type', [i])
		u.put(var.site_pair_name[i][0],	p + 'Site1_Name', [i])
		u.put(var.site_pair_name[i][1],	p + 'Site2_Name', [i])
		u.put(1.12246204830937,				p + 'Cutoff', [i])
		u.put(1.0,							p + 'Scale_1_4_Pair', [i])
		u.put(1.0,							p + 'Force_Cap_LJ.sigma', [i])
		u.put(1.0,							p + 'Force_Cap_LJ.epsilon', [i])
		u.put(r_fc,							p + 'Force_Cap_LJ.r_FC', [i])
	#--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return

##############################################################################
# Force Capped LJ で弱い相互作用を入れて初期緩和
def harmonic(template, read_udf, present_udf, time):
	u = UDFManager(os.path.join(var.target_dir, template))
	# goto global data
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(100000000., 	p + 'Max_Force')
	u.put(time[0], 		p + 'Time.delta_T')
	u.put(time[1], 		p + 'Time.Total_Steps')
	u.put(time[2], 		p + 'Time.Output_Interval_Steps')
	u.put(1.0, 			p + 'Temperature.Temperature')
	u.put(0., 			p + 'Pressure_Stress.Pressure')
	# Calc_Potential_Flags
	p = 'Simulation_Conditions.Calc_Potential_Flags.'
	u.put(1, p + 'Bond')
	u.put(1, p + 'Angle')
	u.put(0, p + 'Non_Bonding_Interchain')
	u.put(0, p + 'Non_Bonding_1_3')
	u.put(0, p + 'Non_Bonding_1_4')
	u.put(0, p + 'Non_Bonding_Intrachain')
	#--- Initial_Structure ---
	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', p+'Method')
	u.put([read_udf, -1, 1, 0], p+'Restart')
	# Relaxation
	p = 'Initial_Structure.Relaxation.'
	u.put(0, p + 'Relaxation')
	#--- Simulation_Conditions ---
	# bond
	for i, b_name in enumerate(var.bond_name):
		p = 'Molecular_Attributes.Bond_Potential[].'
		u.put(b_name, 		p + 'Name', [i])
		u.put('Harmonic', 	p + 'Potential_Type', [i])
		u.put(0.97,			p + 'R0', [i])
		u.put(1000, 		p + 'Harmonic.K', [i])
	# Angle
	for i, anglename in enumerate(var.angle_name):
		p = 'Molecular_Attributes.Angle_Potential[].'
		u.put(anglename, 		p + 'Name', [i])
		u.put('Force_Cap_LJ', 	p + 'Potential_Type', [i])
		u.put(73.0, 				p + 'theta0', [i])
		u.put(1.0, 			p + 'Force_Cap_LJ.sigma', [i])
		u.put(1.0, 			p + 'Force_Cap_LJ.epsilon', [i])
		u.put(1.122462, 		p + 'Force_Cap_LJ.cutoff', [i])
		u.put(0.8, 			p + 'Force_Cap_LJ.r_FC', [i])
	#--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return

###################################
# NPT 条件をハーモニックボンドで設定
def npt_setup(template, read_udf, present_udf, time, pressure):
	u = UDFManager(os.path.join(var.target_dir, template))
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(100000000.,	p + 'Max_Force')
	u.put(time[0],		p + 'Time.delta_T')
	u.put(time[1], 		p + 'Time.Total_Steps')
	u.put(time[2], 		p + 'Time.Output_Interval_Steps')
	u.put(1.0, 			p + 'Temperature.Temperature')
	# Pressure
	u.put(pressure, 	p + 'Pressure_Stress.Pressure')
	# Calc_Potential_Flags
	p = 'Simulation_Conditions.Calc_Potential_Flags.'
	u.put(1, p + 'Bond')
	u.put(0, p + 'Angle')
	u.put(1, p + 'Non_Bonding_Interchain')
	u.put(1, p + 'Non_Bonding_1_3')
	u.put(1, p + 'Non_Bonding_1_4')
	u.put(1, p + 'Non_Bonding_Intrachain')
	#--- Initial_Structure ---
	if read_udf == '':
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 		p + 'Method')
		u.put(['', -1, 0, 0], 	p + 'Restart')
	else:
		# Initial_Unit_Cell
		p = 'Initial_Structure.Initial_Unit_Cell.'
		u.put(0, p + 'Density')
		u.put([0, 0, 0, 90.0, 90.0, 90.0], p + 'Cell_Size')
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 			p + 'Method')
		u.put([read_udf, -1, 1, 0], p + 'Restart')
	# Relaxation
	p = 'Initial_Structure.Relaxation.'
	u.put(1, p + 'Relaxation')
	u.put('DYNAMICS', p + 'Method')
	u.put(100, p + 'Max_Relax_Force')
	u.put(10000, p + 'Max_Relax_Steps')
	#--- Simulation_Conditions ---
	# Bond
	p = 'Molecular_Attributes.Bond_Potential[].'		
	for i, bondname in enumerate(var.bond_name):
		u.put(bondname, 	p + 'Name', [i])
		if read_udf == '':
			u.put('Harmonic', 	p + 'Potential_Type', [i])
			u.put(0.97,			p + 'R0', [i])
			u.put(1000, 		p + 'Harmonic.K', [i])
		else:
			u.put('FENE_LJ', 	p + 'Potential_Type', [i])
			u.put(1.0,			p + 'R0', [i])
			u.put(1.5,			p + 'FENE_LJ.R_max', [i])
			u.put(30,			p + 'FENE_LJ.K', [i])
			u.put(1.0,			p + 'FENE_LJ.sigma', [i])
			u.put(1.0,			p + 'FENE_LJ.epsilon', [i])
	# #--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return

#####################################################
# 直前の条件を維持して平衡化
def eq_setup(template, read_udf, present_udf, time):
	u = UDFManager(os.path.join(var.target_dir, template))
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(time[0],  p+'Time.delta_T')
	u.put(time[1],  p+'Time.Total_Steps')
	u.put(time[2],  p+'Time.Output_Interval_Steps')
	# Moment
	u.put(0, p + "Moment.Interval_of_Calc_Moment")
	u.put(0, p + "Moment.Calc_Moment")
	u.put(0, p + "Moment.Stop_Translation")
	u.put(0, p + "Moment.Stop_Rotation")

	#--- Initial_Structure ---
	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', p+'Method')
	u.put([read_udf, -1, 1, 0], p+'Restart')
	p = 'Initial_Structure.Relaxation.'
	u.put(0, p + 'Relaxation')

	#--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return

###########################################################
# グリーン久保計算
def greenkubo_setup(template, read_udf, present_udf, time):
	u = UDFManager(os.path.join(var.target_dir, template))
	u.jump(-1)
	#--- Simulation_Conditions ---
	# Dynamics_Conditions
	p = 'Simulation_Conditions.Dynamics_Conditions.'
	u.put(time[0],  p+'Time.delta_T')
	u.put(time[1],  p+'Time.Total_Steps')
	u.put(time[2],  p+'Time.Output_Interval_Steps')
	# Calc Correlation
	u.put(1, 'Simulation_Conditions.Output_Flags.Correlation_Function.Stress')
	#--- Initial_Structure ---
	# Generate_Method
	p = 'Initial_Structure.Generate_Method.'
	u.put('Restart', p+'Method')
	u.put([read_udf, -1, 1, 0], p+'Restart')
	p = 'Initial_Structure.Relaxation.'
	u.put(0, p + 'Relaxation')
	#--- Write UDF ---
	u.write(os.path.join(var.target_dir, present_udf))
	return
