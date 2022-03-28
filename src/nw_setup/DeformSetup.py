#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *
import sys
import os
import shutil
# import numpy as np
import platform
# import pandas as pd
################################################################################
def main():
	print("このスクリプトを直接読んでも、初期状態が入っていません。")
	exit(1)
################################################################################
class Setup:
	def __init__(self, path_mod, Ver_Cognac, calc_dir, f_data, core, stress_eval, rate_list, time_div, deform_max, res, def_mode):
		self.path_mod = path_mod
		self.Ver = Ver_Cognac
		self.calc_dir = calc_dir
		self.f_data = f_data
		self.core = core
		self.stress_eval = stress_eval
		#
		self.rate_list = rate_list
		self.time_div = time_div
		self.deform_max = deform_max
		self.res = res
		self.def_mode = def_mode
		#
		self.base_udf = "base_uin.udf"

	############################################################################
	##### Main #####
	# 
	def make_all(self):
		# 平衡計算したUDFファイルとその場所を選択
		t_udf, data_dir = self.file_select()
		#
		# func, nu, structure = self.read_data(data_dir)
		#
		base = self.make_base_udf(t_udf)
		#
		self.make_batch(base)
		#
		self.make_script()
		return

	############################################################################
	##### Function #####
	# 平衡計算したUDFファイルとその場所を選択
	def file_select(self):
		param = sys.argv
		if len(param) == 1:
			print("usage: python", param[0], "Honya_out.udf")
			exit(1)
		elif not os.access(param[1],os.R_OK):
			print(param[1], "not exists.")
			exit(1)
		else:
			target = param[1]
			data_dir = os.path.dirname(target)
		return target, data_dir

	# 
	def read_data(self, data_dir):
		t_data = os.path.join(data_dir, self.f_data)
		if not os.path.exists(t_data):
			exit("network.dat is not exist!")
		else:
			with open(t_data, 'r') as f:
				calc_cond = f.readlines()[-1].strip().split('\t')
			func = calc_cond[3]
			nu = calc_cond[4]
			structure = calc_cond[5]
			return func, nu, structure

	# 
	def make_base_udf(self, t_udf):
		if os.path.exists(self.calc_dir):
			print("Use existing dir of ", self.calc_dir)
		else:
			print("Make new dir of ", self.calc_dir)
			os.makedirs(self.calc_dir)
		#
		base = os.path.join(self.calc_dir, self.base_udf)
		print("Readin file = ", t_udf)
		u = UDFManager(t_udf)
		u.jump(1)
		u.eraseRecord(record_pos=-999,record_num=-999)
		u.write(base)
		return base

	# ファイル名を設定し、バッチファイルを作成
	def make_batch(self, base):
		batch = "#!/bin/bash\n"
		#
		for rate in self.rate_list:
			# UDFファイル名を設定
			rate_str = "{0:4.0e}".format(rate)
			if self.def_mode == 'elong':
				uin = 'Elong_rate_' + rate_str + "_uin.udf"
			elif self.def_mode == 'shear':
				uin = 'Shear_rate_' + rate_str + "_uin.udf"
			# 
			batch = self.make_title(batch, "Calculating rate=" + rate_str)
			batch += self.Ver + ' -I ' + uin + ' -O ' + uin.replace("uin", "out") + ' -n ' + str(self.core) +' \n'
			batch += 'python ' + self.stress_eval + '\n'
			udf_in =  os.path.join(self.calc_dir, uin)
			shutil.copy(base, udf_in)
			self.mod_udf(udf_in, rate)
		# バッチファイルを作成
		if self.def_mode == 'elong':
			f_batch = os.path.join(self.calc_dir, '_elong_bat.bat')
		elif self.def_mode == 'shear':
			f_batch = os.path.join(self.calc_dir, '_shear_bat.bat')
		
		with open(f_batch, 'w') as f:
			f.write(batch)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)
		return


	#-----
	def mod_udf(self, udf_in, rate):
		if self.def_mode == 'elong':
			time_1_step = int(self.res/self.time_div/rate)
			time_total = time_1_step*(self.deform_max - 1)/self.res
		elif self.def_mode == 'shear':
			deform_time = self.deform_max/rate
			time_total = round(deform_time/self.time_div)
			time_1_step = int(self.res/self.time_div/rate)
		#
		u = UDFManager(udf_in)
		# goto global data
		u.jump(-1)

		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000.,			p + 'Max_Force')
		u.put(self.time_div,	p + 'Time.delta_T')
		u.put(time_total,		p + 'Time.Total_Steps')
		u.put(time_1_step,		p + 'Time.Output_Interval_Steps')
		u.put(1.0,				p + 'Temperature.Temperature')
		u.put(0, 				p + 'Temperature.Interval_of_Scale_Temp')
		u.put(0,				p + 'Pressure_Stress.Pressure')

		# Deformation
		if self.def_mode == 'elong':
			p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
			u.put('Cell_Deformation', 		p + 'Method')
			u.put('Simple_Elongation', 		p + 'Cell_Deformation.Method')
			u.put('Initial_Strain_Rate', 	p + 'Cell_Deformation.Simple_Elongation.Input_Method')
			u.put(rate,	 		p + 'Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate')
			u.put(0.5, 						p + 'Cell_Deformation.Simple_Elongation.Poisson_Ratio')
			u.put('z', 						p + 'Cell_Deformation.Simple_Elongation.Axis')
			u.put(1, 						p + 'Cell_Deformation.Interval_of_Deform')
			u.put(0, 						p + 'Cell_Deformation.Deform_Atom')
		elif self.def_mode == 'shear':
			p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
			u.put('Lees_Edwards', 		p + 'Method')
			u.put('Steady', 		p + 'Lees_Edwards.Method')
			u.put(rate, 	p + 'Lees_Edwards.Steady.Shear_Rate')
		
		# Output_Flags
		u.put([1, 1, 1], 'Simulation_Conditions.Output_Flags.Structure')

		# Read_Set_of_Molecules
		p = 'Initial_Structure.Read_Set_of_Molecules'
		u.put(['', -1], p)

		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 		p + 'Method')
		u.put(['', -1, 1, 1], 	p + 'Restart')

		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		u.put(0, p + 'Relaxation')
		#--- Write UDF ---
		u.write(udf_in)
		return
	
	###########################
	# ターミナルのタイトルを設定
	def make_title(self, batch, title):
		if platform.system() == "Windows":
			batch += "title " + title + "\n"
		elif platform.system() == "Linux":
			batch += r'echo -ne "\033]0; ' + title + ' \007"' + '\n'
		return batch

	#######################
	# 必要なスクリプトを作成
	def make_script(self):
		script = self.script_content()
		with open(os.path.join(self.calc_dir, self.stress_eval),'w') as f:
			f.write(script)
		return

	# スクリプトの中身
	def script_content(self):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n################################\n'
		script += 'import os \nimport sys \n'
		script += 'path_mod = "' + self.path_mod + '"\n'
		script += 'sys.path.append(path_mod)\n'
		script += 'from network_evaluation import Read_Shear as rs\n################################\n'
		# script += 'func = ' + str(func) + '\nnu = ' + str(nu) + '\nstructure = "' + structure + '"\n'
		script += '################################\n'
		script += 't_udf_list = rs.file_listing()\n'
		script += 'stress_data = rs.calc_stress_all(t_udf_list)\n'
		script += 'target_list = rs.save_data(stress_data, t_udf_list)\n'
		script += 'rs.plot(target_list)\n'
		return script
		
################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()
