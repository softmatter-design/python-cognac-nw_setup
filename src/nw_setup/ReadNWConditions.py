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
	condsetup = ReadCondSetup()
	basic_cond, nw_cond, sim_cond, rnd_cond, target_cond, target_dir = condsetup.read_and_setcondition()

	return basic_cond, nw_cond, sim_cond, rnd_cond, target_cond, target_dir
	
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
					N_histgram:int "ヒストグラムの分割数"
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

###########################################
##### Read and Setup Calculation Conditions
###########################################
class ReadCondSetup:
###########################################
	def __init__(self):
		pass
	#######################################
	# # Read udf and setup initial conditions
	def read_and_setcondition(self):
		dic={'y':True,'yes':True,'q':False,'quit':False}
		while True:
			# read udf
			basic_cond, rnd_cond = self.readconditionudf()
			# select
			target_cond = self.calc_conditions()
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
			target_dir = self.make_dir()
			nw_cond = [self.nw_model, self.strand, self.n_strand, self.n_segments, self.n_cell, self.n_sc]
			sim_cond = [self.entanglement, self.multi, self.density, self.shrinkage, self.expand, self.step_press, self.press_time, self.step_rfc, self.step_rfc_time, self.equilib_repeat, self.equilib_time, self.greenkubo_repeat, self.greenkubo_time]
			return basic_cond, nw_cond, sim_cond, rnd_cond, target_cond, target_dir
		else:
			sys.exit("##### \nQuit !!")

	####################################
	# Read condition udf
	def readconditionudf(self):
		u = UDFManager('calc_condition.udf')
		u.jump(-1)
		##################
		# 使用するCognacのバージョン
		self.ver_cognac = u.get('CalcCond.Cognac_ver')
		# 計算に使用するコア数
		self.core = u.get('CalcCond.Cores')
		# ベースとするUDFの名前
		base_udf = "base_uin.udf"
		blank_udf = self.ver_cognac + '.udf'
		###########
		basic_cond = [self.ver_cognac, blank_udf, base_udf, self.core]
		#######################################################
		## 計算ターゲット

		###################
		## Networkモデルの設定
		self.nw_model = u.get('TargetCond.Model.TargetModel')
		###################
		## Networkモデルの設定
		self.restart = ''
		self.cond_top = []
		n_hist = 0
		#
		if self.nw_model == "Regular":
			self.strand = u.get('TargetCond.Model.Regular.chains')
		elif self.nw_model == "Random":
			self.strand = u.get('TargetCond.Model.Random.chains')
		################
		if self.strand == "3_Chain" or self.strand == "3_Chain_S" or self.strand == "3_Chain_D":
			self.n_strand = 3
		elif self.strand == "4_Chain":
			self.n_strand = 4
		elif self.strand == "5_Chain":
			self.n_strand = 5
		elif self.strand == "6_Chain":
			self.n_strand = 6
		elif self.strand == "7_Chain":
			self.n_strand = 7
		elif self.strand == "8_Chain":
			self.n_strand = 8
		###################
		## ポリマー鎖の設定
		self.n_segments = u.get('TargetCond.NetWork.N_Segments')
		self.n_sc = u.get('TargetCond.NetWork.N_Subchain')
		self.n_cell = u.get('TargetCond.NetWork.N_UnitCells')
		###################
		if self.nw_model == "Random":
			self.calc = u.get('TargetCond.Model.Random.Calc_Topolpgy')
			n_hist = u.get('TargetCond.Model.Random.N_histgram')
			if self.calc == 'Read':
				self.restart = u.get('TargetCond.Model.Random.Read.dir_name')
				if not os.path.exists(os.path.join(self.restart, 'init.pickle')):
					exit("##########\ntarget directory does not exists.")
				elif self.n_strand != int(self.restart.split('_')[0]):
					sys.exit("##########\nnumber of strands: selected n_strand is different from original Calculation.")
				elif self.n_cell != int(self.restart.split('_')[2]):
					sys.exit("##########\nnumber of cells: selected n_cell is different from original Calculation.")
			elif self.calc == 'Calc':
				self.cond_top = u.get('TargetCond.Model.Random.Calc')
		###################
		## 多重度の設定
		if u.get('TargetCond.Multiplisity.Set_or_Calc') == 'Set':
			self.multi = u.get('TargetCond.Multiplisity.Set.Multiplicity')
			self.density = 0
		elif u.get('TargetCond.Multiplisity.Set_or_Calc') == 'Calc':
			self.multi = 0
			self.density = u.get('TargetCond.Multiplisity.Calc.TargetDensity')
		## 収縮に関する設定
		if u.get('TargetCond.Shrinkage.Shrink') == 'Yes':
			if u.get('TargetCond.Shrinkage.Yes.Control') == 'Density':
				self.density = u.get('TargetCond.Shrinkage.Yes.Density.target_density')
				self.shrinkage = 0.
			elif u.get('TargetCond.Shrinkage.Yes.Control') == 'Shrink':
				self.shrinkage = u.get('TargetCond.Shrinkage.Yes.Shrinkage.value')
				self.density = 0.
		elif u.get('TargetCond.Shrinkage.Shrink') == 'No':
			self.shrinkage = 1.
		#####
		self.entanglement = u.get('TargetCond.Entanglement.Type')
		if self.entanglement == 'Entangled':
			self.step_rfc = u.get('TargetCond.Entanglement.Entangled.Step_rfc[]')
			self.step_rfc_time = u.get('TargetCond.Entanglement.Entangled.Time')
			self.expand = 1.0
			self.step_press = []
			self.press_time = []
		elif self.entanglement == 'NO_Entangled':
			self.step_rfc = []
			self.step_rfc_time = []
			self.expand = u.get('TargetCond.Entanglement.NO_Entangled.ExpansionRatio')
			self.step_press = np.round(np.array(u.get('TargetCond.Entanglement.NO_Entangled.StepPress[]')), 5)
			self.press_time = u.get('TargetCond.Entanglement.NO_Entangled.Time')
		##########
		## シミュレーションの条件
		self.equilib_repeat = u.get('SimulationCond.Equilib_Condition.Repeat')
		self.equilib_time = u.get('SimulationCond.Equilib_Condition.Time')
		#####
		self.greenkubo = u.get('SimulationCond.GreenKubo.Calc')
		self.greenkubo_repeat = 0
		self.greenkubo_time = []
		if self.greenkubo == 'Yes':
			self.greenkubo_repeat = u.get('SimulationCond.GreenKubo.Yes.Repeat')
			self.greenkubo_time = u.get('SimulationCond.GreenKubo.Yes.Time')
		#####
		self.l_bond = u.get('SimulationCond.l_bond')
		#####
		if self.n_segments <= 10:
			self.c_n = 1.5
		elif self.n_segments <= 20:
			self.c_n = 1.65
		elif self.n_segments <= 40:
			self.c_n = 1.7
		else:
			self.c_n = 1.75
		#########################################################################################
		
		rnd_cond = [self.restart, self.cond_top, n_hist]

		return basic_cond, rnd_cond

	############################################
	#-----ネットワークポリマーの諸量を計算
	def calc_conditions(self):
		## 計算システムの諸量を計算して、出力
		n_chains, n_beads_unit = self.set_length()
		target_cond = self.init_calc(n_chains, n_beads_unit)

		return target_cond

	#####################
	#
	def set_length(self):
		self.e2e = self.l_bond*((self.n_segments + 1)*self.c_n)**0.5					# 理想鎖状態での末端間距離

		if self.nw_model == "Regular":
			if self.strand == "3_Chain_S":
				n_chains = 12						        					# サブチェインの本数
				n_beads_unit = 8 + self.n_segments*(1 + self.n_sc)*n_chains		# ユニットセル当たりの粒子数
				self.org_unitcell = (2*2**0.5)*self.e2e				        			# 理想鎖状態でのユニットセル長
			elif self.strand == "3_Chain_D":
				n_chains = 24						       
				n_beads_unit = 16 + self.n_segments*(1 + self.n_sc)*n_chains	
				self.org_unitcell = (2*2**0.5)*self.e2e		
			elif self.strand == "4_Chain":
				n_chains = 16						      
				n_beads_unit = 8 + self.n_segments*(1 + self.n_sc)*n_chains	
				self.org_unitcell = (4*3**0.5)*self.e2e/3			 
			elif self.strand == "6_Chain":
				n_chains = 3						  
				n_beads_unit = 1 + self.n_segments*(1 + self.n_sc)*n_chains		
				self.org_unitcell = self.e2e						      
			elif self.strand == "8_Chain":
				n_chains = 8						   
				n_beads_unit = 2 + self.n_segments*(1 + self.n_sc)*n_chains   
				self.org_unitcell = (2*3**0.5)*self.e2e/3	

		elif self.nw_model == "Random":
			n_chains = self.n_strand
			n_beads_unit = 2 + self.n_segments*(1 + self.n_sc)*n_chains
			self.org_unitcell = (2*3**0.5)*self.e2e/3	

		return n_chains, n_beads_unit

	###############################################################
	def init_calc(self, n_chains, n_beads_unit):
		calc_flag = 0
		err_dens = 0
		if self.multi == 0:
			calc_flag = 1
			if self.shrinkage == 1.:
				calcd_multi = round(self.density*self.org_unitcell**3/n_beads_unit)	# 密度を設定値とした場合に必要な多重度
				calcd_density = n_beads_unit*calcd_multi/self.org_unitcell**3		# 上記の多重度での密度
				err_dens = round((calcd_density/self.density - 1)*100, 2) 		# 設定密度との誤差(%)
				single_net_atom = int(n_beads_unit*self.n_cell**3.)	    		# 一つのネットワーク中の粒子数
				self.total_net_atom = int(calcd_multi*single_net_atom)    			# 全システム中のネットワーク粒子数
				mod_unitcell = (calcd_multi*n_beads_unit/self.density)**(1/3)					# その際のシステムサイズ
				self.shrinkage = mod_unitcell/self.org_unitcell								# 収縮比
				self.system = mod_unitcell*self.n_cell
				# vol = self.system**3									    			# システム体積
				print(self.n_cell)
				self.multi = calcd_multi
				mod_e2e = self.shrinkage*self.e2e											# 収縮後の末端間距離
				unit_cell = mod_unitcell
				self.nu = n_chains*self.multi/unit_cell**3
			elif self.shrinkage != 1.:
				sys.exit(u"\n############################################## \n多重度を自動計算にした場合、収縮条件は選択できません\n条件設定を見直してください。\n##############################################\n")
		elif self.multi != 0:
			if self.shrinkage == 1.:
				err_dens = 0.
				self.system = self.org_unitcell*self.n_cell						# e2e から決めたシステムサイズ
				self.density = n_beads_unit*self.multi/self.org_unitcell**3	# 多重度での密度
				single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
				self.total_net_atom = int(self.multi*single_net_atom)    	# 全システム中のネットワーク粒子数
				# vol = self.system**3									    	# システム体積
				mod_e2e = self.e2e											# 収縮後の末端間距離
				unit_cell = self.org_unitcell
				self.nu = n_chains*self.multi/unit_cell**3
			elif self.shrinkage != 1.:
				if self.density == 0.:
					err_dens = 0.
					mod_unitcell = self.org_unitcell*self.shrinkage
					self.system = mod_unitcell*self.n_cell						# e2e から決めたシステムサイズ
					self.density = n_beads_unit*self.multi/mod_unitcell**3	# 多重度での密度
					single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
					self.total_net_atom = int(self.multi*single_net_atom)    	# 全システム中のネットワーク粒子数
					# vol = system**3									    	# システム体積
					mod_e2e = self.shrinkage*self.e2e		
					unit_cell = mod_unitcell									# 収縮後の末端間距離
					self.nu = n_chains*self.multi/unit_cell**3
				elif self.density > 0.:
					err_dens = 0.
					mod_unitcell = (n_beads_unit*self.multi/self.density)**(1/3)
					self.shrinkage = mod_unitcell/self.org_unitcell
					single_net_atom = int(n_beads_unit*self.n_cell**3.)	    # 一つのネットワーク中の粒子数
					self.total_net_atom = int(self.multi*single_net_atom)    	# 全システム中のネットワーク粒子数
					self.system = mod_unitcell*self.n_cell						# e2e から決めたシステムサイズ
					# vol = self.system**3									    	# システム体積
					mod_e2e = self.shrinkage*self.e2e											# 収縮後の末端間距離
					unit_cell = mod_unitcell
					self.nu = n_chains*self.multi/unit_cell**3
		else:
			sys.exit("Something Wrong!!")
		#
		text = "#########################################" + "\n"
		text += "計算に使用するコア数\t\t" + str(self.core ) + "\n"
		text += "#########################################" + "\n"
		text += "ネットワークトポロジー\t\t" + str(self.nw_model) + "\n"
		text += "ネットワークモデル\t\t" + str(self.strand) + "\n"
		if self.nw_model == "Random":
			if self.calc == 'Read':
				text += "\t** 過去の計算を読み込み **\n"
				text += "Directory:" + str(self.restart) + "\n"
			elif self.calc == 'Calc':
				text += "\t** ランダム構造を計算 **\n"
				text += "ランダム構造の計算条件\t" + str(self.cond_top) + "\n"
		text += "#########################################" + "\n"
		text += "ストランド中のセグメント数:\t" + str(self.n_segments) + "\n"
		text += "特性比:\t\t\t\t" + str(round(self.c_n, 2)) + "\n"
		text += "初期の末端間距離:\t\t" + str(round(self.e2e, 4)) + "\n"
		text += "当初の単位ユニット:\t\t" + str(round(self.org_unitcell, 4)) + "\n"
		text += "一辺当たりの単位ユニット数:\t" + str(self.n_cell) + "\n"
		# text += "当初のシステムサイズ:\t\t" + str(round(org_system, 4)) + "\n"
		text += "#########################################" + "\n"
		if calc_flag == 1:
			text += "設定密度:\t\t\t" + str(self.density) + "\n"
			text += "算出された多重度:\t\t" + str(self.multi) + "\n"
			text += "上記の多重度での密度:\t\t" + str(round(calcd_density, 4)) + "\n"
			text += "密度の誤差:\t\t\t" + str(round(err_dens, 3)) + "%" + "\n"
			text += "収縮比:\t\t\t\t" + str(round(self.shrinkage, 4)) + "\n"
		else:
			text += "多重度:\t\t\t\t" + str(self.multi) + "\n"
			text += "密度:\t\t\t\t" + str(round(self.density, 4)) + "\n"
			text += "収縮比:\t\t\t\t" + str(round(self.shrinkage, 4)) + "\n"
		text += "NW の全セグメント数:\t\t" + str(self.total_net_atom) + "\n"
		text += "システムサイズ:\t\t\t" + str(round(self.system, 4)) + "\n"
		text += "#########################################" + "\n"
		text += "絡み合いの有無:\t\t\t" + str(self.entanglement) + "\n"
		if self.entanglement == 'Entangled':
			text += "Slow Push Off 条件: " + ', '.join(map(str, self.step_rfc)) + "\n"
			text += "Slow Push Off 時間条件: " + str(self.step_rfc_time) + "\n"
		else:
			text += "NPT 計算時の初期膨張率:\t\t" + str(self.expand) + "\n"
			text += "ステップ圧力:\t" + ', '.join(map(str, self.step_press)) + "\n"
			text += "圧力時間条件:\t\t" + str(self.press_time) + "\n"
		text += "#########################################" + "\n"
		text += "平衡化計算繰り返し:\t\t" + str(self.equilib_repeat) + "\n"
		text += "平衡化時間条件:\t\t" + str(self.equilib_time ) + "\n"
		if self.greenkubo == 'Yes':
			text += "応力緩和計算繰り返し:\t\t" + str(self.greenkubo_repeat) + "\n"
			text += "応力緩和時間条件:\t" + str(self.greenkubo_time) + "\n"
		text += "#########################################" + "\n"
		text += "ストランドの数密度:\t\t" + str(round(self.nu, 5)) + "\n"
		text += "#########################################" + "\n"
		print(text)

		if abs(err_dens) > 1:
			print(u"############################################## \n圧縮後の密度が、設定したい密度と 1% 以上違います。\nそれでも計算しますか？\n##############################################\n")

		target_cond = [self.system, unit_cell, self.total_net_atom]

		return target_cond

	################################################################################
	# 計算用のディレクトリーを作成
	def make_dir(self):
		target_name = self.nw_model + "_" + self.entanglement + "_" + self.strand + '_N_' + str(self.n_segments) + "_Cells_" + str(self.n_cell) + "_Multi_" + str(self.multi)
		os.makedirs(target_name, exist_ok = True)
		# with open(os.path.join(target_name, "calc.dat"), "w") as f:
		# 	f.write("# segments\tbond_length\tCN\tfunc\tnu\tNW_type\n" + str(self.n_segments) + '\t' + str(self.l_bond) + '\t' + str(self.c_n) + "\t" + str(round(self.nu, 5)) + '\t' + self.nw_model)
		self.make_cond_udf(target_name)
		return target_name
	
	###########################################
	# make new udf when not found.
	def make_cond_udf(self, target_name):
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
		u.put(self.ver_cognac, 'CalcCond.Cognac_ver')
		u.put(self.core, 'CalcCond.Cores')
		##################
		u.put(self.nw_model, 'TargetCond.Model.TargetModel')
		u.put('none', 'TargetCond.Model.RandomData')
		u.put([], 'TargetCond.Model.SelectedValue[]')

		u.put(self.strand, 'TargetCond.NetWork.Strand')
		u.put(self.n_strand, 'TargetCond.NetWork.N_Strands')
		u.put(self.n_segments, 'TargetCond.NetWork.N_Segments')
		u.put(self.n_sc, 'TargetCond.NetWork.N_Subchain')
		u.put(self.n_cell, 'TargetCond.NetWork.N_UnitCells')

		u.put(self.c_n, 'TargetCond.Strand.Characteristic_Ratio')
		u.put(self.e2e, 'TargetCond.Strand.R')
		u.put(self.org_unitcell, 'TargetCond.Strand.Initial_Unit_Cell')

		u.put(self.multi, 'TargetCond.Multiplisity.Multiplicity')

		u.put(self.shrinkage, 'TargetCond.Shrinkage.Shrinkage')
		u.put(self.density, 'TargetCond.Shrinkage.Density')

		u.put(self.entanglement, 'TargetCond.Entanglement.Type')

		u.put(self.total_net_atom, 'TargetCond.System.Total_Segments')
		u.put(self.system, 'TargetCond.System.SystemSize')
		u.put(self.nu, 'TargetCond.System.Nu')
		###################
		u.put(self.equilib_repeat, 'SimulationCond.Equilib_Condition.repeat')
		u.put(self.equilib_time[0], 'SimulationCond.Equilib_Condition.Time.delta_T')
		u.put(self.equilib_time[1], 'SimulationCond.Equilib_Condition.Time.Total_Steps')
		u.put(self.equilib_time[2], 'SimulationCond.Equilib_Condition.Time.Output_Interval_Steps')

		u.put(self.l_bond, 'SimulationCond.l_bond')
		##################
		u.write()

		return
