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
	def __init__(self, path_mod, ver_Cognac, f_data, core, evaluate, rate, lamd_list, cycle_lamd, resol, calc_mode, cycle):
		self.path_mod = path_mod
		self.ver = ver_Cognac
		self.f_data = f_data
		self.core = core
		self.evaluate = evaluate
		#
		self.rate = rate
		self.time_div = 0.01
		self.lamd_list = lamd_list
		self.cycle_lamd = cycle_lamd
		self.resol = resol
		self.calc_mode = calc_mode
		self.cycle = cycle
		#
		self.base_udf = "base_uin.udf"

	############################################################################
	##### Main #####
	# 
	def make_all(self):
		# 平衡計算したUDFファイルとその場所を選択
		t_udf, data_dir = self.file_select()
		#
		func, nu, structure = self.read_data(data_dir)
		#
		calc_dir = self.make_batch(t_udf)
		#
		self.make_script(calc_dir, func, nu, structure)
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
			exit("calc.dat is not exist!")
		else:
			with open(t_data, 'r') as f:
				calc_cond = f.readlines()[-1].strip().split('\t')
			func = calc_cond[3]
			nu = calc_cond[4]
			structure = calc_cond[5]
			return func, nu, structure

	# ファイル名を設定し、バッチファイルを作成
	def make_batch(self, t_udf):
		batch = "#!/bin/bash\n"
		#
		if self.calc_mode == 'cycle':
			# base_udfとcalc_dirを作成し、system_sizeを戻す。
			base, calc_dir = self.make_base_udf(t_udf)
			for n_cycle in range(self.cycle):
				lamd = self.cycle_lamd
				lamd_str = str(lamd).replace('.', '_')
				# batch fileに追加
				batch = self.make_title(batch, "Calculating Cycle: Lambda=" + lamd_str)
				for direction in ['forward', 'backward']:
					uin = 'Cycle_' + str(n_cycle) + '_' + lamd_str + '_' + direction + '_uin.udf'
					uout = uin.replace("uin", "out")
					batch += self.ver + ' -I ' + uin + ' -O ' + uout + ' -n ' + str(self.core) +' \n'
					batch += 'python ' + self.evaluate + ' ' + uout + '\n'
					shutil.copy(base, os.path.join(calc_dir, uin))
					if n_cycle == 0 and direction == 'forward':
						self.mod_udf(os.path.join(calc_dir, uin), '', direction, lamd)
					else:
						self.mod_udf(os.path.join(calc_dir, uin), pre, direction, lamd)
					pre = uout

		elif self.calc_mode == 'step':
			# base_udfとcalc_dirを作成し、system_sizeを戻す。
			base, calc_dir = self.make_base_udf(t_udf)
			for i, lamd in enumerate(self.lamd_list):
				lamd_str = str(lamd).replace('.', '_')
				# batch fileに追加
				batch = self.make_title(batch, "Calculating Step: Lambda=" + lamd_str)
				for direction in ['forward', 'backward']:
					uin = 'Step_' + lamd_str + '_' + direction + '_uin.udf'
					uout = uin.replace("uin", "out")
					batch += self.ver + ' -I ' + uin + ' -O ' + uout + ' -n ' + str(self.core) +' \n'
					batch += 'python ' + self.evaluate + ' ' + uout + '\n'
					shutil.copy(base, os.path.join(calc_dir, uin))
					if i == 0 and direction == 'forward':
						self.mod_udf(os.path.join(calc_dir, uin), '', direction, lamd)
					else:
						self.mod_udf(os.path.join(calc_dir, uin), pre, direction, lamd)
					pre = uout

		# バッチファイルを作成
		f_batch = os.path.join(calc_dir, '_hysteresis.bat')
		with open(f_batch, 'w') as f:
			f.write(batch)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)
		return calc_dir

	def make_base_udf(self, t_udf):
		if self.calc_mode == 'cycle':
			lamd_str = str(self.cycle_lamd).replace('.', '_')
			rate_str = str("{0:4.0e}".format(self.rate)).replace('.', '_')
			calc_dir = 'Hysteresis_' + self.calc_mode + '_rate_' + rate_str + '_lambda_' + lamd_str
		if os.path.exists(calc_dir):
			print("Use existing dir of ", calc_dir)
		else:
			print("Make new dir of ", calc_dir)
			os.makedirs(calc_dir)
		#
		base = os.path.join(calc_dir, self.base_udf)
		print("Readin file = ", t_udf)
		udf = UDFManager(t_udf)
		self.mod_base_udf(udf)
		udf.write(base)
		return base, calc_dir

	def mod_base_udf(self, udf):
		udf.jump(1)
		system = udf.get('Structure.Unit_Cell.Cell_Size.c')
		udf.eraseRecord(record_pos=-999,record_num=-999)
		# goto global data
		udf.jump(-1)
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		udf.put(100000.,		p + 'Max_Force')
		udf.put(self.time_div,	p + 'Time.delta_T')
		udf.put(10,				p + 'Time.Total_Steps')
		udf.put(1,				p + 'Time.Output_Interval_Steps')
		udf.put(1.0,			p + 'Temperature.Temperature')
		udf.put(0, 				p + 'Temperature.Interval_of_Scale_Temp')
		udf.put(0,				p + 'Pressure_Stress.Pressure')
		# Deformation
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		udf.put('Cell_Deformation', 	p + 'Method')
		udf.put('Simple_Elongation', 	p + 'Cell_Deformation.Method')
		udf.put('Deformation_Speed', 	p + 'Cell_Deformation.Simple_Elongation.Input_Method')
		udf.put(system*self.rate,	 	p + 'Cell_Deformation.Simple_Elongation.Deformation_Speed.Speed')
		udf.put(0.5, 					p + 'Cell_Deformation.Simple_Elongation.Poisson_Ratio')
		udf.put('z', 					p + 'Cell_Deformation.Simple_Elongation.Axis')
		udf.put(1, 						p + 'Cell_Deformation.Interval_of_Deform')
		udf.put(0, 						p + 'Cell_Deformation.Deform_Atom')
		# # Read_Set_of_Molecules
		# p = 'Initial_Structure.Read_Set_of_Molecules'
		# udf.put([self.base_udf, -1], p)
		# Generate_Method
		# p = 'Initial_Structure.Generate_Method.'
		# udf.put('Restart', 			p + 'Method')
		# udf.put(['', -1, 1, 0], 	p + 'Restart')
		# Relaxation
		p = 'Initial_Structure.Relaxation.'
		udf.put(0, p + 'Relaxation')
		return
	
	def mod_udf(self, uin, pre, direction, lamd):
		time_1_step = round(self.resol/self.time_div/self.rate, -2)
		time_total = round(time_1_step*(lamd - 1)/self.resol, -2)
		#
		udf = UDFManager(uin)
		udf.jump(-1)
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		udf.put(time_total,		p + 'Time.Total_Steps')
		udf.put(time_1_step,	p + 'Time.Output_Interval_Steps')

		if direction == 'backward':
			# Deformation
			p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
			srate = udf.get(p + 'Cell_Deformation.Simple_Elongation.Deformation_Speed.Speed')
			udf.put(-1*srate,	 p + 'Cell_Deformation.Simple_Elongation.Deformation_Speed.Speed')
		# Generate_Method
		p = 'Initial_Structure.Generate_Method.'
		udf.put('Restart', 			p + 'Method')
		udf.put([pre, -1, 1, 0], 	p + 'Restart')
		#--- Write UDF ---
		udf.write(uin)
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
	def make_script(self, calc_dir, func, nu, structure):
		script = self.script_content(func, nu, structure)
		with open(os.path.join(calc_dir, self.evaluate), 'w') as f:
			f.write(script)
		return

	# スクリプトの中身
	def script_content(self, func, nu, structure):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n################################\n'
		script += 'import os \nimport sys \n'
		script += 'path_mod = "' + self.path_mod + '"\n'
		script += 'sys.path.append(path_mod)\n'
		script += 'from network_evaluation import Read_Stress2 as rs\n################################\n'
		script += 'func = ' + str(func) + '\nnu = ' + str(nu) + '\nstructure = "' + structure + '"\n'
		script += '################################\n'
		script += 'rs.calc_all(func, nu, structure)\n'
		return script
		
################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()
