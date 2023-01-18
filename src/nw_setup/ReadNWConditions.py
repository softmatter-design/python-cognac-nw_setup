#!/usr/bin/env python
# -*- coding: utf-8 -*-
#######################################################
import os
import sys
import codecs
import numpy as np
from UDFManager import UDFManager
import nw_setup.variables as var
#######################################################
#
def setupcondition():
	# check 'NW_setup.udf' and make it.
	findudf()
	# Read udf and setup initial conditions.
	read_and_setcondition()

	return
	
###########################################
# check 'NW_setup.udf' and make it.
def findudf():
	if not os.path.isfile('./NW_setup.udf'):
		print()
		print('In this directory, no "NW_setup.udf" is found !')
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
			Strand:{Char:select{"KG", "FENE", "Harmonic"} "ストランドの種類"
				} "ストランドの条件を設定"
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
	{"KG"}
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
	with codecs.open('./NW_setup.udf', 'w', 'utf_8') as f:
		f.write(contents)
	return

#######################################
# # Read udf and setup initial conditions
def read_and_setcondition():
	dic={'y':True,'yes':True,'q':False,'quit':False}
	while True:
		# read udf
		readconditionudf()
		# select
		calc_conditions()
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
		make_dir()

		return
	else:
		sys.exit("##### \nQuit !!")

####################################
# Read condition udf
def readconditionudf():
	u = UDFManager('NW_setup.udf')
	u.jump(-1)
	##################
	# 使用するCognacのバージョン
	var.ver_cognac = u.get('CalcCond.Cognac_ver')
	# 計算に使用するコア数
	var.core = u.get('CalcCond.Cores')
	# ベースとするUDFの名前
	var.base_udf = "base_uin.udf"
	var.blank_udf = var.ver_cognac + '.udf'

	#######################################################
	## 計算ターゲット

	###################
	## Networkモデルの設定
	var.nw_model = u.get('TargetCond.Model.TargetModel')
	###################
	## Networkモデルの設定
	if var.nw_model == "Regular":
		var.strand = u.get('TargetCond.Model.Regular.chains')
	elif var.nw_model == "Random":
		var.strand = u.get('TargetCond.Model.Random.chains')
	################
	if var.strand == "3_Chain" or var.strand == "3_Chain_S" or var.strand == "3_Chain_D":
		var.n_strand = 3
	elif var.strand == "4_Chain":
		var.n_strand = 4
	elif var.strand == "5_Chain":
		var.n_strand = 5
	elif var.strand == "6_Chain":
		var.n_strand = 6
	elif var.strand == "7_Chain":
		var.n_strand = 7
	elif var.strand == "8_Chain":
		var.n_strand = 8
	###################
	## ポリマー鎖の設定
	var.n_segments = u.get('TargetCond.NetWork.N_Segments')
	var.n_sc = u.get('TargetCond.NetWork.N_Subchain')
	var.n_cell = u.get('TargetCond.NetWork.N_UnitCells')
	###################
	if var.nw_model == "Random":
		var.calc = u.get('TargetCond.Model.Random.Calc_Topolpgy')
		var.histgram_bins = u.get('TargetCond.Model.Random.histgram_bins')
		if var.calc == 'Read':
			var.restart = u.get('TargetCond.Model.Random.Read.dir_name')
			if not os.path.exists(os.path.join(var.restart, 'init.pickle')):
				exit("##########\ntarget directory does not exists.")
			elif var.n_strand != int(var.restart.split('_')[0]):
				sys.exit("##########\nnumber of strands: selected n_strand is different from original Calculation.")
			elif var.n_cell != int(var.restart.split('_')[2]):
				sys.exit("##########\nnumber of cells: selected n_cell is different from original Calculation.")
		elif var.calc == 'Calc':
			var.cond_top = u.get('TargetCond.Model.Random.Calc')
	#####
	var.strand_type = u.get('TargetCond.Strand.Char')
	###################
	## 多重度の設定
	if u.get('TargetCond.Multiplisity.Set_or_Calc') == 'Set':
		var.multi_org = u.get('TargetCond.Multiplisity.Set.Multiplicity')
		var.density_org = 0
	elif u.get('TargetCond.Multiplisity.Set_or_Calc') == 'Calc':
		var.multi_org = 0
		var.density_org = u.get('TargetCond.Multiplisity.Calc.TargetDensity')
	## 収縮に関する設定
	if u.get('TargetCond.Shrinkage.Shrink') == 'Yes':
		if u.get('TargetCond.Shrinkage.Yes.Control') == 'Density':
			var.density_org = u.get('TargetCond.Shrinkage.Yes.Density.target_density')
			var.shrinkage = 1.
		elif u.get('TargetCond.Shrinkage.Yes.Control') == 'Shrink':
			var.shrinkage = u.get('TargetCond.Shrinkage.Yes.Shrinkage.value')
			var.density_org = 0.
	elif u.get('TargetCond.Shrinkage.Shrink') == 'No':
		var.shrinkage = 0.
	#####
	var.entanglement = u.get('TargetCond.Entanglement.Type')
	if var.entanglement == 'Entangled':
		var.step_rfc = u.get('TargetCond.Entanglement.Entangled.Step_rfc[]')
		var.step_rfc_time = u.get('TargetCond.Entanglement.Entangled.Time')
		var.expand = 1.0
	elif var.entanglement == 'NO_Entangled':
		var.expand = u.get('TargetCond.Entanglement.NO_Entangled.ExpansionRatio')
		var.step_press = np.round(np.array(u.get('TargetCond.Entanglement.NO_Entangled.StepPress[]')), 5)
		var.press_time = u.get('TargetCond.Entanglement.NO_Entangled.Time')
	##########
	## シミュレーションの条件
	var.equilib_repeat = u.get('SimulationCond.Equilib_Condition.Repeat')
	var.equilib_time = u.get('SimulationCond.Equilib_Condition.Time')
	#####
	var.greenkubo = u.get('SimulationCond.GreenKubo.Calc')
	var.greenkubo_repeat = 0
	var.greenkubo_time = []
	if var.greenkubo == 'Yes':
		var.greenkubo_repeat = u.get('SimulationCond.GreenKubo.Yes.Repeat')
		var.greenkubo_time = u.get('SimulationCond.GreenKubo.Yes.Time')
	#####
	var.l_bond = u.get('SimulationCond.l_bond')
	#####
	if var.n_segments <= 10:
		var.c_n  = 1.5
	elif var.n_segments <= 20:
		var.c_n  = 1.65
	elif var.n_segments <= 40:
		var.c_n  = 1.7
	else:
		var.c_n  = 1.75

	return

############################################
#-----ネットワークポリマーの諸量を計算
def calc_conditions():
	## 計算システムの諸量を計算して、出力
	set_length()
	init_calc()
	return

#####################
#
def set_length():
	var.e2e = var.l_bond*((var.n_segments + 1)*var.c_n )**0.5					# 理想鎖状態での末端間距離

	if var.nw_model == "Regular":
		if var.strand == "3_Chain_S":
			var.n_chains = 12						        					# サブチェインの本数
			var.n_beads_unit = 8 + var.n_segments*(1 + var.n_sc)*var.n_chains		# ユニットセル当たりの粒子数
			var.org_unitcell = (2*2**0.5)*var.e2e				        			# 理想鎖状態でのユニットセル長
		elif var.strand == "3_Chain_D":
			var.n_chains = 24						       
			var.n_beads_unit = 16 + var.n_segments*(1 + var.n_sc)*var.n_chains	
			var.org_unitcell = (2*2**0.5)*var.e2e		
		elif var.strand == "4_Chain":
			var.n_chains = 16						      
			var.n_beads_unit = 8 + var.n_segments*(1 + var.n_sc)*var.n_chains	
			var.org_unitcell = (4*3**0.5)*var.e2e/3			 
		elif var.strand == "6_Chain":
			var.n_chains = 3						  
			var.n_beads_unit = 1 + var.n_segments*(1 + var.n_sc)*var.n_chains		
			var.org_unitcell = var.e2e						      
		elif var.strand == "8_Chain":
			var.n_chains = 8						   
			var.n_beads_unit = 2 + var.n_segments*(1 + var.n_sc)*var.n_chains   
			var.org_unitcell = (2*3**0.5)*var.e2e/3	

	elif var.nw_model == "Random":
		var.n_chains = var.n_strand
		var.n_beads_unit = 2 + var.n_segments*(1 + var.n_sc)*var.n_chains
		var.org_unitcell = (2*3**0.5)*var.e2e/3	

	return

###############################################################
def init_calc():
	calc_flag = 0
	err_dens = 0
	if var.multi_org == 0:
		calc_flag = 1
		if var.shrinkage == 0.:
			calcd_multi = round(var.density_org*var.org_unitcell**3/var.n_beads_unit)	# 密度を設定値とした場合に必要な多重度
			calcd_density = var.n_beads_unit*calcd_multi/var.org_unitcell**3		# 上記の多重度での密度
			err_dens = round((calcd_density/var.density_org - 1)*100, 2) 		# 設定密度との誤差(%)
			single_net_atom = int(var.n_beads_unit*var.n_cell**3.)	    		# 一つのネットワーク中の粒子数
			var.total_net_atom = int(calcd_multi*single_net_atom)    			# 全システム中のネットワーク粒子数
			mod_unitcell = (calcd_multi*var.n_beads_unit/var.density_org)**(1/3)					# その際のシステムサイズ
			var.shrinkage = mod_unitcell/var.org_unitcell								# 収縮比
			var.system = mod_unitcell*var.n_cell
			# vol = var.system**3									    			# システム体積
			# print(var.n_cell)
			var.multi_mod = calcd_multi
			mod_e2e = var.shrinkage*var.e2e											# 収縮後の末端間距離
			var.unit_cell = mod_unitcell
			var.nu = var.n_chains*var.multi_mod/var.unit_cell**3
			# print(var.nu, var.n_chains, var.multi_mod, var.unit_cell)
			var.density_mod = var.density_org
		elif var.shrinkage != 0.:
			sys.exit(u"\n############################################## \n多重度を自動計算にした場合、収縮条件は選択できません\n条件設定を見直してください。\n##############################################\n")
	elif var.multi_org != 0:
		var.multi_mod = var.multi_org
		if var.shrinkage == 0.:
			err_dens = 0.
			var.system = var.org_unitcell*var.n_cell						# e2e から決めたシステムサイズ
			var.density_mod = var.n_beads_unit*var.multi_org/var.org_unitcell**3	# 多重度での密度
			single_net_atom = int(var.n_beads_unit*var.n_cell**3.)	    # 一つのネットワーク中の粒子数
			var.total_net_atom = int(var.multi_org*single_net_atom)    	# 全システム中のネットワーク粒子数
			# vol = var.system**3									    	# システム体積
			mod_e2e = var.e2e											# 収縮後の末端間距離
			var.unit_cell = var.org_unitcell
			var.nu = var.n_chains*var.multi_org/var.unit_cell**3
		elif var.shrinkage != 0.:
			if var.density_org == 0.:
				err_dens = 0.
				mod_unitcell = var.org_unitcell*var.shrinkage
				var.system = mod_unitcell*var.n_cell						# e2e から決めたシステムサイズ
				var.density_mod = var.n_beads_unit*var.multi_org/mod_unitcell**3	# 多重度での密度
				single_net_atom = int(var.n_beads_unit*var.n_cell**3.)	    # 一つのネットワーク中の粒子数
				var.total_net_atom = int(var.multi_org*single_net_atom)    	# 全システム中のネットワーク粒子数
				# vol = var.system**3									    	# システム体積
				mod_e2e = var.shrinkage*var.e2e		
				var.unit_cell = mod_unitcell									# 収縮後の末端間距離
				var.nu = var.n_chains*var.multi_org/var.unit_cell**3
			elif var.density_org > 0.:
				err_dens = 0.
				mod_unitcell = (var.n_beads_unit*var.multi_org/var.density_org)**(1/3)
				var.shrinkage = mod_unitcell/var.org_unitcell
				single_net_atom = int(var.n_beads_unit*var.n_cell**3.)	    # 一つのネットワーク中の粒子数
				var.total_net_atom = int(var.multi_org*single_net_atom)    	# 全システム中のネットワーク粒子数
				var.system = mod_unitcell*var.n_cell						# e2e から決めたシステムサイズ
				# vol = var.system**3									    	# システム体積
				mod_e2e = var.shrinkage*var.e2e											# 収縮後の末端間距離
				var.unit_cell = mod_unitcell
				var.nu = var.n_chains*var.multi_org/var.unit_cell**3
				var.density_mod = var.density_org
	else:
		sys.exit("Something Wrong!!")
	#
	text = "#########################################" + "\n"
	text += "計算に使用するコア数\t\t" + str(var.core ) + "\n"
	text += "#########################################" + "\n"
	text += "ネットワークトポロジー\t\t" + str(var.nw_model) + "\n"
	text += "ネットワークモデル\t\t" + str(var.strand) + "\n"
	if var.nw_model == "Random":
		if var.calc == 'Read':
			text += "\t** 過去の計算を読み込み **\n"
			text += "Directory:" + str(var.restart) + "\n"
		elif var.calc == 'Calc':
			text += "\t** ランダム構造を計算 **\n"
			text += "ランダム構造の計算条件\t" + str(var.cond_top) + "\n"
	text += "#########################################" + "\n"
	text += "ストランド中のセグメント数:\t" + str(var.n_segments) + "\n"
	text += "ストランドの種類:\t\t" + str(var.strand_type) + "\n"
	if var.strand_type == "KG":
		text += "\tbond: FENE-LJ\n"
		text += "\tangle: LJ\n"
		text += "\tsegment Int.: LJ\n"
	elif var.strand_type == "FENE":
		text += "\tbond: FENE-LJ\n"
		text += "\tangle: force Capped LJ\n"
		text += "\tsegment Int.: No\n"
	elif var.strand_type == "Harmonic":
		text += "\tbond: Harmonic\n"
		text += "\tangle: force Capped LJ\n"
		text += "\tsegment Int.: No\n"
	text += "特性比:\t\t\t\t" + str(round(var.c_n , 2)) + "\n"
	text += "初期の末端間距離:\t\t" + str(round(var.e2e, 4)) + "\n"
	text += "当初の単位ユニット:\t\t" + str(round(var.org_unitcell, 4)) + "\n"
	text += "一辺当たりの単位ユニット数:\t" + str(var.n_cell) + "\n"
	text += "#########################################" + "\n"
	if calc_flag == 1:
		text += "設定密度:\t\t\t" + str(var.density_mod) + "\n"
		text += "算出された多重度:\t\t" + str(var.multi_mod) + "\n"
		text += "上記の多重度での密度:\t\t" + str(round(calcd_density, 4)) + "\n"
		text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
		text += "収縮比:\t\t\t\t" + str(round(var.shrinkage, 4)) + "\n"
	else:
		text += "多重度:\t\t\t\t" + str(var.multi_org) + "\n"
		text += "密度:\t\t\t\t" + str(round(var.density_mod, 4)) + "\n"
		text += "収縮比:\t\t\t\t" + str(round(var.shrinkage, 4)) + "\n"
	text += "NW の全セグメント数:\t\t" + str(var.total_net_atom) + "\n"
	text += "システムサイズ:\t\t\t" + str(round(var.system, 4)) + "\n"
	text += "#########################################" + "\n"
	text += "絡み合いの有無:\t\t\t" + str(var.entanglement) + "\n"
	if var.entanglement == 'Entangled':
		text += "Slow Push Off 条件: " + ', '.join(map(str, var.step_rfc)) + "\n"
		text += "Slow Push Off 時間条件: " + str(var.step_rfc_time) + "\n"
	else:
		text += "NPT 計算時の初期膨張率:\t\t" + str(var.expand) + "\n"
		text += "ステップ圧力:\t" + ', '.join(map(str, var.step_press)) + "\n"
		text += "圧力時間条件:\t\t" + str(var.press_time) + "\n"
	text += "#########################################" + "\n"
	text += "平衡化計算繰り返し:\t\t" + str(var.equilib_repeat) + "\n"
	text += "平衡化時間条件:\t\t" + str(var.equilib_time ) + "\n"
	if var.greenkubo == 'Yes':
		text += "応力緩和計算繰り返し:\t\t" + str(var.greenkubo_repeat) + "\n"
		text += "応力緩和時間条件:\t" + str(var.greenkubo_time) + "\n"
	text += "#########################################" + "\n"
	text += "ストランドの数密度:\t\t" + str(round(var.nu, 5)) + "\n"
	text += "#########################################" + "\n"
	print(text)

	if abs(err_dens) > 1:
		print(u"############################################## \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？\n##############################################\n")
	return

################################################################################
# 計算用のディレクトリーを作成
def make_dir():
	var.target_dir = f"{var.nw_model:}_{var.entanglement:}_{var.strand:}_{var.strand_type:}_N_{var.n_segments:}_Cells_{var.n_cell:}_Multi_{var.multi_mod:}"
	var.target_name = f"{var.strand:}_{var.strand_type:}_N_{var.n_segments:}_"
	os.makedirs(var.target_dir, exist_ok = True)
	make_cond_udf()
	return

###########################################
# make new udf when not found.
def make_cond_udf():
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
			StrandType:{Type: string "ストランドの種類",
					Bond: string "ボンドの種類",
					Interaction: string "セグメント間相互作用"
					Angle: string "アングルポテンシャル"
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
			{"", "", "", ""}
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
	with codecs.open(os.path.join(var.target_dir, 'target_condition.udf'), 'w', 'utf_8') as f:
		f.write(contents)

	###
	u = UDFManager(os.path.join(var.target_dir, 'target_condition.udf'))
	u.jump(-1)
	##################
	u.put(var.ver_cognac, 'CalcCond.Cognac_ver')
	u.put(var.core, 'CalcCond.Cores')
	##################
	u.put(var.nw_model, 'TargetCond.Model.TargetModel')
	u.put('none', 'TargetCond.Model.RandomData')
	u.put([], 'TargetCond.Model.SelectedValue[]')

	u.put(var.strand, 'TargetCond.NetWork.Strand')
	u.put(var.n_strand, 'TargetCond.NetWork.N_Strands')
	u.put(var.n_segments, 'TargetCond.NetWork.N_Segments')
	u.put(var.n_sc, 'TargetCond.NetWork.N_Subchain')
	u.put(var.n_cell, 'TargetCond.NetWork.N_UnitCells')

	u.put(var.c_n , 'TargetCond.Strand.Characteristic_Ratio')
	u.put(var.e2e, 'TargetCond.Strand.R')
	u.put(var.org_unitcell, 'TargetCond.Strand.Initial_Unit_Cell')

	u.put(var.strand_type , 'TargetCond.StrandType.Type')
	if var.strand_type == 'KG':
		u.put('FENE-LJ' , 'TargetCond.StrandType.Bond')
		u.put('LJ' , 'TargetCond.StrandType.Interaction')
		u.put('LJ' , 'TargetCond.StrandType.Angle')
	elif var.strand_type == 'FENE':
		u.put('FENE-LJ' , 'TargetCond.StrandType.Bond')
		u.put('No' , 'TargetCond.StrandType.Interaction')
		u.put('forceCappedLJ' , 'TargetCond.StrandType.Angle')
	elif var.strand_type == 'Harmonic':
		u.put('Harmonic' , 'TargetCond.StrandType.Bond')
		u.put('No' , 'TargetCond.StrandType.Interaction')
		u.put('forceCappedLJ' , 'TargetCond.StrandType.Angle')

	u.put(var.multi_mod, 'TargetCond.Multiplisity.Multiplicity')

	u.put(var.shrinkage, 'TargetCond.Shrinkage.Shrinkage')
	u.put(var.density_mod, 'TargetCond.Shrinkage.Density')

	u.put(var.entanglement, 'TargetCond.Entanglement.Type')

	u.put(var.total_net_atom, 'TargetCond.System.Total_Segments')
	u.put(var.system, 'TargetCond.System.SystemSize')
	u.put(var.nu, 'TargetCond.System.Nu')
	###################
	u.put(var.equilib_repeat, 'SimulationCond.Equilib_Condition.repeat')
	u.put(var.equilib_time[0], 'SimulationCond.Equilib_Condition.Time.delta_T')
	u.put(var.equilib_time[1], 'SimulationCond.Equilib_Condition.Time.Total_Steps')
	u.put(var.equilib_time[2], 'SimulationCond.Equilib_Condition.Time.Output_Interval_Steps')

	u.put(var.l_bond, 'SimulationCond.l_bond')
	##################
	u.write()

	return
