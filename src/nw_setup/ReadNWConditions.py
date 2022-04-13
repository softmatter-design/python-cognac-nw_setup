#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################
import os
import sys
import codecs
import numpy as np
from UDFManager import UDFManager
#######################################################
#
def setupcondition():
	# check 'calc_condition.udf' and make it.
	findudf()
	# Read udf and setup initial conditions.
	basic_cond, nw_cond, sim_cond, sim_cond2, rnd_cond, target_cond, target_dir = read_and_setcondition()

	return basic_cond, nw_cond, sim_cond, sim_cond2, rnd_cond, target_cond, target_dir
	
###########################################
# check 'calc_condition.udf' and make it.
def findudf():
	if not os.path.isfile('./calc_condition.udf'):
		print()
		print('In this directory, no "calc_condition.udf" is found !')
		print('New one will be generated.')
		print('Please, modify and save it !\n')
		makenewudf()
		input('Press ENTER to continue...')
	return

###########################################
# make new udf when not found.
def makenewudf():
	contents = '''
	\\begin{def}
		CalcCond:{
			Cognac_ver:select{"cognac112"} "使用する Cognac のバージョン",
			Cores: int "計算に使用するコア数を指定"
			} "計算の条件を設定"
		TargetCond:{
			Model:{TargetModel:select{"Regular", "Random"} "ネットワークのモデルを選択",
				Regular:{chains:select{"3_Chain_S", "3_Chain_D", "4_Chain", "6_Chain", "8_Chain"} "分岐の数と種類を選択"
					} "規則構造での条件を入力",
				Random:{chains:select{"3_Chain", "4_Chain", "5_Chain", "6_Chain", "7_Chain"} "分岐の数と種類を選択",
					Calc_Topolpgy:select{"Calc", "Read"} "ランダムネットワークの「計算を行うか、読み込むか」を選択",
						Calc:{pre_sampling:int "プレサンプリング数", 
						sampling:int "サンプリング数", 
						try:int "サンプリング時の再トライ数", 
						repeat:int "探索計算の繰り返し数"
						n_parallel:int "並行計算のCPU数"} "ランダムサーチ計算する場合の条件を設定",
						Read:{dir_name:string} "過去の計算結果のディレクトリを記入",
					histgram_bins:int "ヒストグラムの分割数"
					} "ランダム構造での条件を入力"
				} "シミュレーションの条件を設定"
			NetWork:{N_Segments: int "ストランド中のセグメント数", 
					N_Subchain: int "各セグメントの側鎖の数", 
					N_UnitCells: int "一辺あたりのユニットセルの数"
				} "ネットワークの条件を設定"
			Multiplisity:{Set_or_Calc:select{"Set", "Calc"} "多重度を設定するかどうかのフラッグ",
					Set:{Multiplicity: int} "多重度を設定",
					Calc:{TargetDensity:float} "多重度を自動設定した場合の密度を設定 \\n設定した密度になるように多重度を設定"
				} "多重度設定に関する設定"
			Shrinkage:{Shrink:select{"Yes", "No"} "ストランドを自然長から圧縮するかどうかのフラッグ \\n非圧縮時には、多重度に応じて密度が変化",
				Yes:{Control:select{"Density", "Shrink"} "圧縮する場合に、密度コントロールにするか、圧縮率を決めるかを設定", 
				Density:{target_density: float} "目標とする密度を設定", 
				Shrinkage:{value: float} "ストランドの圧縮比率を設定"
				}
				} "ストランドを自然長から圧縮するかどうかを設定"
			Entanglement:{
				Type:select{"Entangled", "NO_Entangled"} "ネットワーク・トポロジーを選択",
					Entangled:{Step_rfc[]: float "Slow Push Off での rfc 条件",
						Time:{delta_T: double, Total_Steps: int, Output_Interval_Steps: int} "時間条件を入力"} "密度、末端間距離を設定値に合わせるように多重度を自動設定。\\n絡み合いが入るように初期化",
					NO_Entangled:{
						ExpansionRatio: float "NPT 計算での初期膨張率", 
						StepPress[]: float "NPT 計算での圧力変化",
						Time:{delta_T: double, Total_Steps: int, Output_Interval_Steps: int} "時間条件を入力"
						} "密度、末端間距離を設定値に合わせるように多重度を自動設定。\\n絡み合いが入らないようにNPTで縮める。"
				} "ネットワーク・トポロジーを選択",
			} "計算ターゲットの条件を設定"
		SimulationCond:{
			Equilib_Condition:{
					Repeat: int "平衡化計算の繰り返し数",
					Time:{delta_T: double, Total_Steps: int, Output_Interval_Steps: int} "平衡化計算の時間条件を入力"
				} "平衡化計算の時間条件を入力",
			GreenKubo:{
				Calc:select{"Yes", "No"},
				Yes:{
					Repeat:int "計算の繰り返し数",
					Time:{delta_T: double, Total_Steps: int, Output_Interval_Steps: int} "時間条件を入力"
					} "GreenKubo により、応力緩和関数を計算するかどうかを決める。"
				}
			l_bond: float "シミュレーションでのボンドの自然長"
			} "シミュレーションの条件を設定"
	\end{def}

	\\begin{data}
		CalcCond:{"cognac112",1}
TargetCond:{
	{"Regular", {"4_Chain"}{"4_Chain","Read",{1000,100,100,10,1}{"4_chains_3_cells_100_trials_100_sampling"}100}}
	{20, 0, 3}
	{"Set", {1}{0.85}}
	{"No", {"Density", {0.85}{1.0}}}
	{"NO_Entangled",
		{[1.073,1.0,0.9,0.8], {1.0e-02,300000,2000}},
		{2.0, [0.2,0.5,1.0,2.0,3.0,4.5], {1.0e-02,300000,2000}}
		}
	}
SimulationCond:{
	{4,{1.0e-02,1000000,10000}}
	{"Yes",{5,{1.0e-02,1000000,10000}}}
	0.97
	}

\end{data}
	'''
	###
	with codecs.open('./calc_condition.udf', 'w', 'utf_8') as f:
		f.write(contents)
	return

#######################################
# # Read udf and setup initial conditions
def read_and_setcondition():
	dic={'y':True,'yes':True,'q':False,'quit':False}
	while True:
		# read udf
		basic_cond, rnd_cond, nw_cond, sim_cond = readconditionudf()
		# select
		mod_cond, target_cond, sim_cond2 = calc_conditions(basic_cond, nw_cond, rnd_cond, sim_cond)
		print('Change UDF: type [r]eload')
		print('Quit input process: type [q]uit')
		inp = input('Condition is OK ==> [y]es >> ').lower()
		if inp in dic:
			inp = dic[inp]
			break
		print('##### \nRead Condition UDF again \n#####\n\n')
	if inp:
		print("\n\nSetting UP progress !!")
		# 計算用のディレクトリーを作成
		target_dir = make_dir(basic_cond, nw_cond, sim_cond, sim_cond2, mod_cond, target_cond)

		return basic_cond, nw_cond, sim_cond, sim_cond2, rnd_cond, target_cond, target_dir
	else:
		sys.exit("##### \nQuit !!")

####################################
# Read condition udf
def readconditionudf():
	u = UDFManager('calc_condition.udf')
	u.jump(-1)
	##################
	# 使用するCognacのバージョン
	ver_cognac = u.get('CalcCond.Cognac_ver')
	# 計算に使用するコア数
	core = u.get('CalcCond.Cores')
	# ベースとするUDFの名前
	base_udf = "base_uin.udf"
	blank_udf = ver_cognac + '.udf'
	###########
	basic_cond = [ver_cognac, blank_udf, base_udf, core]
	#######################################################
	## 計算ターゲット

	###################
	## Networkモデルの設定
	nw_model = u.get('TargetCond.Model.TargetModel')
	###################
	## Networkモデルの設定
	restart = ''
	cond_top = []
	histgram_bins = 0
	#
	if nw_model == "Regular":
		strand = u.get('TargetCond.Model.Regular.chains')
	elif nw_model == "Random":
		strand = u.get('TargetCond.Model.Random.chains')
	################
	if strand == "3_Chain" or strand == "3_Chain_S" or strand == "3_Chain_D":
		n_strand = 3
	elif strand == "4_Chain":
		n_strand = 4
	elif strand == "5_Chain":
		n_strand = 5
	elif strand == "6_Chain":
		n_strand = 6
	elif strand == "7_Chain":
		n_strand = 7
	elif strand == "8_Chain":
		n_strand = 8
	###################
	## ポリマー鎖の設定
	n_segments = u.get('TargetCond.NetWork.N_Segments')
	n_sc = u.get('TargetCond.NetWork.N_Subchain')
	n_cell = u.get('TargetCond.NetWork.N_UnitCells')
	###################
	calc = ''
	if nw_model == "Random":
		calc = u.get('TargetCond.Model.Random.Calc_Topolpgy')
		histgram_bins = u.get('TargetCond.Model.Random.histgram_bins')
		if calc == 'Read':
			restart = u.get('TargetCond.Model.Random.Read.dir_name')
			if not os.path.exists(os.path.join(restart, 'init.pickle')):
				exit("##########\ntarget directory does not exists.")
			elif n_strand != int(restart.split('_')[0]):
				sys.exit("##########\nnumber of strands: selected n_strand is different from original Calculation.")
			elif n_cell != int(restart.split('_')[2]):
				sys.exit("##########\nnumber of cells: selected n_cell is different from original Calculation.")
		elif calc == 'Calc':
			cond_top = u.get('TargetCond.Model.Random.Calc')
	###################
	## 多重度の設定
	if u.get('TargetCond.Multiplisity.Set_or_Calc') == 'Set':
		multi_org = u.get('TargetCond.Multiplisity.Set.Multiplicity')
		density_org = 0
	elif u.get('TargetCond.Multiplisity.Set_or_Calc') == 'Calc':
		multi_org = 0
		density_org = u.get('TargetCond.Multiplisity.Calc.TargetDensity')
	## 収縮に関する設定
	if u.get('TargetCond.Shrinkage.Shrink') == 'Yes':
		if u.get('TargetCond.Shrinkage.Yes.Control') == 'Density':
			density_org = u.get('TargetCond.Shrinkage.Yes.Density.target_density')
			shrinkage = 1.
		elif u.get('TargetCond.Shrinkage.Yes.Control') == 'Shrink':
			shrinkage = u.get('TargetCond.Shrinkage.Yes.Shrinkage.value')
			density_org = 0.
	elif u.get('TargetCond.Shrinkage.Shrink') == 'No':
		shrinkage = 0.
	#####
	entanglement = u.get('TargetCond.Entanglement.Type')
	if entanglement == 'Entangled':
		step_rfc = u.get('TargetCond.Entanglement.Entangled.Step_rfc[]')
		step_rfc_time = u.get('TargetCond.Entanglement.Entangled.Time')
		expand = 1.0
		step_press = []
		press_time = []
	elif entanglement == 'NO_Entangled':
		step_rfc = []
		step_rfc_time = []
		expand = u.get('TargetCond.Entanglement.NO_Entangled.ExpansionRatio')
		step_press = np.round(np.array(u.get('TargetCond.Entanglement.NO_Entangled.StepPress[]')), 5)
		press_time = u.get('TargetCond.Entanglement.NO_Entangled.Time')
	##########
	## シミュレーションの条件
	equilib_repeat = u.get('SimulationCond.Equilib_Condition.Repeat')
	equilib_time = u.get('SimulationCond.Equilib_Condition.Time')
	#####
	greenkubo = u.get('SimulationCond.GreenKubo.Calc')
	greenkubo_repeat = 0
	greenkubo_time = []
	if greenkubo == 'Yes':
		greenkubo_repeat = u.get('SimulationCond.GreenKubo.Yes.Repeat')
		greenkubo_time = u.get('SimulationCond.GreenKubo.Yes.Time')
	#####
	l_bond = u.get('SimulationCond.l_bond')
	#####
	if n_segments <= 10:
		c_n  = 1.5
	elif n_segments <= 20:
		c_n  = 1.65
	elif n_segments <= 40:
		c_n  = 1.7
	else:
		c_n  = 1.75
	#########################################################################################
	
	rnd_cond = [restart, cond_top, histgram_bins]

	nw_cond = [nw_model, strand, n_strand, n_segments, n_cell, n_sc, l_bond, c_n]

	sim_cond = [entanglement, multi_org, density_org, shrinkage, expand, step_press, press_time, step_rfc, step_rfc_time, equilib_repeat, equilib_time, greenkubo, greenkubo_repeat, greenkubo_time, calc]

	return basic_cond, rnd_cond, nw_cond, sim_cond

############################################
#-----ネットワークポリマーの諸量を計算
def calc_conditions(basic_cond, nw_cond, rnd_cond, sim_cond):
	## 計算システムの諸量を計算して、出力
	mod_cond = set_length(nw_cond)
	target_cond, sim_cond2 = init_calc(basic_cond, nw_cond, rnd_cond, sim_cond, mod_cond)

	return mod_cond, target_cond, sim_cond2

#####################
#
def set_length(nw_cond):
	nw_model = nw_cond[0]
	strand = nw_cond[1]
	n_strand = nw_cond[2]
	n_segments = nw_cond[3]
	n_cell = nw_cond[4]
	n_sc = nw_cond[5]
	l_bond = nw_cond[6]
	c_n = nw_cond[7]

	e2e = l_bond*((n_segments + 1)*c_n )**0.5					# 理想鎖状態での末端間距離

	if nw_model == "Regular":
		if strand == "3_Chain_S":
			n_chains = 12						        					# サブチェインの本数
			n_beads_unit = 8 + n_segments*(1 + n_sc)*n_chains		# ユニットセル当たりの粒子数
			org_unitcell = (2*2**0.5)*e2e				        			# 理想鎖状態でのユニットセル長
		elif strand == "3_Chain_D":
			n_chains = 24						       
			n_beads_unit = 16 + n_segments*(1 + n_sc)*n_chains	
			org_unitcell = (2*2**0.5)*e2e		
		elif strand == "4_Chain":
			n_chains = 16						      
			n_beads_unit = 8 + n_segments*(1 + n_sc)*n_chains	
			org_unitcell = (4*3**0.5)*e2e/3			 
		elif strand == "6_Chain":
			n_chains = 3						  
			n_beads_unit = 1 + n_segments*(1 + n_sc)*n_chains		
			org_unitcell = e2e						      
		elif strand == "8_Chain":
			n_chains = 8						   
			n_beads_unit = 2 + n_segments*(1 + n_sc)*n_chains   
			org_unitcell = (2*3**0.5)*e2e/3	

	elif nw_model == "Random":
		n_chains = n_strand
		n_beads_unit = 2 + n_segments*(1 + n_sc)*n_chains
		org_unitcell = (2*3**0.5)*e2e/3	

	mod_cond = [n_chains, n_beads_unit, e2e, org_unitcell]
	return mod_cond

###############################################################
def init_calc(basic_cond, nw_cond, rnd_cond, sim_cond, mod_cond):
	ver_cognac = basic_cond[0]
	blank_udf = basic_cond[1]
	base_udf = basic_cond[2]
	core = basic_cond[3]

	nw_model = nw_cond[0]
	strand = nw_cond[1]
	n_strand = nw_cond[2]
	n_segments = nw_cond[3]
	n_cell = nw_cond[4]
	n_sc = nw_cond[5]
	l_bond = nw_cond[6]
	c_n = nw_cond[7]

	restart = rnd_cond[0]
	cond_top = rnd_cond[1]

	entanglement  = sim_cond[0]
	multi_org = sim_cond[1]
	density_org = sim_cond[2] 
	shrinkage = sim_cond[3]
	expand = sim_cond[4]
	step_press = sim_cond[5]
	press_time = sim_cond[6]
	step_rfc = sim_cond[7]
	step_rfc_time = sim_cond[8]
	equilib_repeat = sim_cond[9]
	equilib_time = sim_cond[10]
	greenkubo = sim_cond[11]
	greenkubo_repeat = sim_cond[12]
	greenkubo_time = sim_cond[13]
	calc = sim_cond[14]

	mod_cond[0]
	n_chains = mod_cond[0]
	n_beads_unit = mod_cond[1]
	e2e = mod_cond[2]
	org_unitcell = mod_cond[3]

	calc_flag = 0
	err_dens = 0
	if multi_org == 0:
		calc_flag = 1
		if shrinkage == 0.:
			calcd_multi = round(density_org*org_unitcell**3/n_beads_unit)	# 密度を設定値とした場合に必要な多重度
			calcd_density = n_beads_unit*calcd_multi/org_unitcell**3		# 上記の多重度での密度
			err_dens = round((calcd_density/density_org - 1)*100, 2) 		# 設定密度との誤差(%)
			single_net_atom = int(n_beads_unit*n_cell**3.)	    		# 一つのネットワーク中の粒子数
			total_net_atom = int(calcd_multi*single_net_atom)    			# 全システム中のネットワーク粒子数
			mod_unitcell = (calcd_multi*n_beads_unit/density_org)**(1/3)					# その際のシステムサイズ
			shrinkage = mod_unitcell/org_unitcell								# 収縮比
			system = mod_unitcell*n_cell
			# vol = system**3									    			# システム体積
			# print(n_cell)
			multi_mod = calcd_multi
			mod_e2e = shrinkage*e2e											# 収縮後の末端間距離
			unit_cell = mod_unitcell
			nu = n_chains*multi_org/unit_cell**3
			density_mod = density_org
		elif shrinkage != 0.:
			sys.exit(u"\n############################################## \n多重度を自動計算にした場合、収縮条件は選択できません\n条件設定を見直してください。\n##############################################\n")
	elif multi_org != 0:
		multi_mod = multi_org
		if shrinkage == 0.:
			err_dens = 0.
			system = org_unitcell*n_cell						# e2e から決めたシステムサイズ
			density_mod = n_beads_unit*multi_org/org_unitcell**3	# 多重度での密度
			single_net_atom = int(n_beads_unit*n_cell**3.)	    # 一つのネットワーク中の粒子数
			total_net_atom = int(multi_org*single_net_atom)    	# 全システム中のネットワーク粒子数
			# vol = system**3									    	# システム体積
			mod_e2e = e2e											# 収縮後の末端間距離
			unit_cell = org_unitcell
			nu = n_chains*multi_org/unit_cell**3
		elif shrinkage != 0.:
			if density_org == 0.:
				err_dens = 0.
				mod_unitcell = org_unitcell*shrinkage
				system = mod_unitcell*n_cell						# e2e から決めたシステムサイズ
				density_mod = n_beads_unit*multi_org/mod_unitcell**3	# 多重度での密度
				single_net_atom = int(n_beads_unit*n_cell**3.)	    # 一つのネットワーク中の粒子数
				total_net_atom = int(multi_org*single_net_atom)    	# 全システム中のネットワーク粒子数
				# vol = system**3									    	# システム体積
				mod_e2e = shrinkage*e2e		
				unit_cell = mod_unitcell									# 収縮後の末端間距離
				nu = n_chains*multi_org/unit_cell**3
			elif density_org > 0.:
				err_dens = 0.
				mod_unitcell = (n_beads_unit*multi_org/density_org)**(1/3)
				shrinkage = mod_unitcell/org_unitcell
				single_net_atom = int(n_beads_unit*n_cell**3.)	    # 一つのネットワーク中の粒子数
				total_net_atom = int(multi_org*single_net_atom)    	# 全システム中のネットワーク粒子数
				system = mod_unitcell*n_cell						# e2e から決めたシステムサイズ
				# vol = system**3									    	# システム体積
				mod_e2e = shrinkage*e2e											# 収縮後の末端間距離
				unit_cell = mod_unitcell
				nu = n_chains*multi_org/unit_cell**3
				density_mod = density_org
	else:
		sys.exit("Something Wrong!!")
	#
	text = "#########################################" + "\n"
	text += "計算に使用するコア数\t\t" + str(core ) + "\n"
	text += "#########################################" + "\n"
	text += "ネットワークトポロジー\t\t" + str(nw_model) + "\n"
	text += "ネットワークモデル\t\t" + str(strand) + "\n"
	if nw_model == "Random":
		if calc == 'Read':
			text += "\t** 過去の計算を読み込み **\n"
			text += "Directory:" + str(restart) + "\n"
		elif calc == 'Calc':
			text += "\t** ランダム構造を計算 **\n"
			text += "ランダム構造の計算条件\t" + str(cond_top) + "\n"
	text += "#########################################" + "\n"
	text += "ストランド中のセグメント数:\t" + str(n_segments) + "\n"
	text += "特性比:\t\t\t\t" + str(round(c_n , 2)) + "\n"
	text += "初期の末端間距離:\t\t" + str(round(e2e, 4)) + "\n"
	text += "当初の単位ユニット:\t\t" + str(round(org_unitcell, 4)) + "\n"
	text += "一辺当たりの単位ユニット数:\t" + str(n_cell) + "\n"
	# text += "当初のシステムサイズ:\t\t" + str(round(org_system, 4)) + "\n"
	text += "#########################################" + "\n"
	if calc_flag == 1:
		text += "設定密度:\t\t\t" + str(density_mod) + "\n"
		text += "算出された多重度:\t\t" + str(multi_mod) + "\n"
		text += "上記の多重度での密度:\t\t" + str(round(calcd_density, 4)) + "\n"
		text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
		text += "収縮比:\t\t\t\t" + str(round(shrinkage, 4)) + "\n"
	else:
		text += "多重度:\t\t\t\t" + str(multi_org) + "\n"
		text += "密度:\t\t\t\t" + str(round(density_mod, 4)) + "\n"
		text += "収縮比:\t\t\t\t" + str(round(shrinkage, 4)) + "\n"
	text += "NW の全セグメント数:\t\t" + str(total_net_atom) + "\n"
	text += "システムサイズ:\t\t\t" + str(round(system, 4)) + "\n"
	text += "#########################################" + "\n"
	text += "絡み合いの有無:\t\t\t" + str(entanglement) + "\n"
	if entanglement == 'Entangled':
		text += "Slow Push Off 条件: " + ', '.join(map(str, step_rfc)) + "\n"
		text += "Slow Push Off 時間条件: " + str(step_rfc_time) + "\n"
	else:
		text += "NPT 計算時の初期膨張率:\t\t" + str(expand) + "\n"
		text += "ステップ圧力:\t" + ', '.join(map(str, step_press)) + "\n"
		text += "圧力時間条件:\t\t" + str(press_time) + "\n"
	text += "#########################################" + "\n"
	text += "平衡化計算繰り返し:\t\t" + str(equilib_repeat) + "\n"
	text += "平衡化時間条件:\t\t" + str(equilib_time ) + "\n"
	if greenkubo == 'Yes':
		text += "応力緩和計算繰り返し:\t\t" + str(greenkubo_repeat) + "\n"
		text += "応力緩和時間条件:\t" + str(greenkubo_time) + "\n"
	text += "#########################################" + "\n"
	text += "ストランドの数密度:\t\t" + str(round(nu, 5)) + "\n"
	text += "#########################################" + "\n"
	print(text)

	if abs(err_dens) > 1:
		print(u"############################################## \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？\n##############################################\n")

	target_cond = [system, unit_cell, total_net_atom, nu]

	sim_cond2 = [multi_mod, density_mod]

	return target_cond, sim_cond2

################################################################################
# 計算用のディレクトリーを作成
def make_dir(basic_cond, nw_cond, sim_cond, sim_cond2, mod_cond, target_cond):
	nw_model = nw_cond[0]
	strand = nw_cond[1]
	n_strand = nw_cond[2]
	n_segments = nw_cond[3]
	n_cell = nw_cond[4]
	n_sc = nw_cond[5]
	l_bond = nw_cond[6]

	entanglement  = sim_cond[0]

	multi_mod = sim_cond2[0]

	nu = target_cond[3]
	target_name = nw_model + "_" + entanglement + "_" + strand + '_N_' + str(n_segments) + "_Cells_" + str(n_cell) + "_Multi_" + str(multi_mod)
	os.makedirs(target_name, exist_ok = True)
	# with open(os.path.join(target_name, "calc.dat"), "w") as f:
	# 	f.write("# segments\tbond_length\tCN\tfunc\tnu\tNW_type\n" + str(n_segments) + '\t' + str(l_bond) + '\t' + str(c_n ) + "\t" + str(round(nu, 5)) + '\t' + nw_model)
	make_cond_udf(basic_cond, target_name, nw_cond, sim_cond, sim_cond2, mod_cond, target_cond)
	return target_name

###########################################
# make new udf when not found.
def make_cond_udf(basic_cond, target_name, nw_cond, sim_cond, sim_cond2, mod_cond, target_cond):
	ver_cognac = basic_cond[0]
	blank_udf = basic_cond[1]
	base_udf = basic_cond[2]
	core = basic_cond[3]

	nw_model = nw_cond[0]
	strand = nw_cond[1]
	n_strand = nw_cond[2]
	n_segments = nw_cond[3]
	n_cell = nw_cond[4]
	n_sc = nw_cond[5]
	l_bond = nw_cond[6]
	c_n = nw_cond[7]

	entanglement  = sim_cond[0]
	shrinkage = sim_cond[3]
	equilib_repeat = sim_cond[9]
	equilib_time = sim_cond[10]

	multi_mod = sim_cond2[0]
	density_mod = sim_cond2[1] 

	e2e = mod_cond[2]
	org_unitcell = mod_cond[3]

	system = target_cond[0]
	unit_cell = target_cond[1]
	total_net_atom = target_cond[2]
	nu = target_cond[3]
	
	contents = '''
	\\begin{def}
		CalcCond:{
			Cognac_ver:string "使用する Cognac のバージョン",
			Cores: int "計算に使用するコア数を指定"
			} "計算の条件を設定"
		TargetCond:{
			Model:{
				TargetModel: string "ネットワークのモデル",
				RandomData: string "ランダム計算のディレクトリ",
				SelectedValue[]: float "ランダム構造の場合の選択"
				} "シミュレーションの条件を設定"
			NetWork:{
				Strand: string "分岐の数と種類",
				N_Strands: int "ストランドの数"
				N_Segments: int "ストランド中のセグメント数", 
				N_Subchain: int "各セグメントの側鎖の数", 
				N_UnitCells: int "一辺あたりのユニットセルの数"
				} "ネットワークの条件を設定"
			Strand:{Characteristic_Ratio: float "特性比",
					R: float "自然長",
					Initial_Unit_Cell: float "自然長でのユニットセル長さ"
				}
			Multiplisity:{
				Multiplicity: int
				} "多重度設定に関する設定"
			Shrinkage:{
				Shrinkage: float, 
				Density:float
				} "ストランドを自然長から圧縮するかどうかを設定"
			Entanglement:{
				Type: string
				} "ネットワーク・トポロジーを選択",
			System:{
				Total_Segments: int "全セグメント数",
				SystemSize: float,
				Nu: float
			}
			} "計算ターゲットの条件を設定"
		SimulationCond:{
			Equilib_Condition:{
					repeat: int "平衡化計算の繰り返し数",
					Time:{delta_T: double, Total_Steps: int, Output_Interval_Steps: int
					} "平衡化計算の時間条件を入力"
				} "平衡化計算の時間条件を入力",
			l_bond: float "シミュレーションでのボンドの自然長"
			} "シミュレーションの条件を設定"
	\end{def}

	\\begin{data}
		CalcCond:{"",1}
		TargetCond:{
			{"", "", []}
			{"", 4, 0, 0, 0}
			{0., 1., 1.}
			{1}
			{1.0, 0.85}
			{""}
			{1, 1., 1.}
			}
		SimulationCond:{
			{4,{1.0e-02,100000,1000}}
			0.97
			}

		\end{data}
	'''
	###
	with codecs.open(os.path.join(target_name, 'target_condition.udf'), 'w', 'utf_8') as f:
		f.write(contents)

	###
	u = UDFManager(os.path.join(target_name, 'target_condition.udf'))
	u.jump(-1)
	##################
	u.put(ver_cognac, 'CalcCond.Cognac_ver')
	u.put(core, 'CalcCond.Cores')
	##################
	u.put(nw_model, 'TargetCond.Model.TargetModel')
	u.put('none', 'TargetCond.Model.RandomData')
	u.put([], 'TargetCond.Model.SelectedValue[]')

	u.put(strand, 'TargetCond.NetWork.Strand')
	u.put(n_strand, 'TargetCond.NetWork.N_Strands')
	u.put(n_segments, 'TargetCond.NetWork.N_Segments')
	u.put(n_sc, 'TargetCond.NetWork.N_Subchain')
	u.put(n_cell, 'TargetCond.NetWork.N_UnitCells')

	u.put(c_n , 'TargetCond.Strand.Characteristic_Ratio')
	u.put(e2e, 'TargetCond.Strand.R')
	u.put(org_unitcell, 'TargetCond.Strand.Initial_Unit_Cell')

	u.put(multi_mod, 'TargetCond.Multiplisity.Multiplicity')

	u.put(shrinkage, 'TargetCond.Shrinkage.Shrinkage')
	u.put(density_mod, 'TargetCond.Shrinkage.Density')

	u.put(entanglement, 'TargetCond.Entanglement.Type')

	u.put(total_net_atom, 'TargetCond.System.Total_Segments')
	u.put(system, 'TargetCond.System.SystemSize')
	u.put(nu, 'TargetCond.System.Nu')
	###################
	u.put(equilib_repeat, 'SimulationCond.Equilib_Condition.repeat')
	u.put(equilib_time[0], 'SimulationCond.Equilib_Condition.Time.delta_T')
	u.put(equilib_time[1], 'SimulationCond.Equilib_Condition.Time.Total_Steps')
	u.put(equilib_time[2], 'SimulationCond.Equilib_Condition.Time.Output_Interval_Steps')

	u.put(l_bond, 'SimulationCond.l_bond')
	##################
	u.write()

	return
