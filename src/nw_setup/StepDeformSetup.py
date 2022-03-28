#!/usr/bin/env python
# -*- coding: utf-8 -*-
##### Import #####
from UDFManager import *
import sys
import os
import math
import shutil
import numpy as np
import platform
import pandas as pd
################################################################################
def main():
	print("このスクリプトを直接読んでも、初期状態が入っていません。")
	exit(1)
################################################################################
class Setup:
	def __init__(self, py_mod, ver_cognac, f_data, core, step_cond, relax, repeat):
		self.py_mod = py_mod
		self.ver = ver_cognac
		self.f_data = f_data
		self.core = ' -n ' + str(core) +' \n'
		#
		self.mode = step_cond[0]
		self.step_rate = step_cond[1]
		self.deform_max = step_cond[2]
		#
		self.relax_max = relax
		#
		self.repeat = repeat
		#
		if self.mode == 'elong':
			self.calc_dir = "Step_Elong_Lambda_" + str(self.deform_max).replace('.', '_') + '_rate_' + str(self.step_rate).replace('.', '_')
		elif self.mode == 'shear':
			self.calc_dir = "Step_Shear_Gamma_" + str(self.deform_max).replace('.', '_') + '_rate_' + str(self.step_rate).replace('.', '_')
		#
		self.base_udf = "base_uin.udf"

	############################################################################
	##### Main #####
	# 
	def calc_step(self):
		# 計算条件を出力し、確認
		conditions = self.calc_cond()
		# 平衡計算したUDFファイルとその場所を選択
		base = self.select_udf()
		# 
		self.setup_batch(conditions, base)

		return

	############################################################################
	# シミュレーションコンディションをプリント
	def calc_cond(self):
		def_time = 0.
		conditions = []
		# Step Deformation
		temp_scale = 1
		delta_t = 1e-3
		if self.mode == 'elong':
			def_rate = 5e-2
			if def_rate != self.step_rate:
				print(self.step_rate)
				print('おすすめの変形レートと異なります。')
				self.prompt()
				rate = self.step_rate
			else:
				rate = def_rate
			deform_time = (self.deform_max - 1.)/rate
		elif self.mode == 'shear':
			def_rate = 2e0
			if def_rate != self.step_rate:
				print('おすすめの変形レートと異なります。')
				self.prompt(self)
				rate = self.step_rate
			else:
				rate = def_rate
			deform_time = self.deform_max/rate
		t_steps = round(deform_time/delta_t)
		interval = 1
		for interval in [1, 2, 5, 10]:
			if round(t_steps/interval) <= 2000:
				break
		print(t_steps, interval, deform_time)
		conditions.append([delta_t, t_steps, interval, deform_time, temp_scale])
		# Quench
		temp_scale = 0
		quench_cond = []
		#
		if interval < 10:
			# 1st
			max_t = 0.1
			t_steps = round(max_t/delta_t)
			quench_cond.append([delta_t, t_steps, interval, max_t, temp_scale])
		# Quench_all
		delta_t = 1e-2
		cycles = 100
		for max_t in [1, 1e1, 1e2, 1e3, 1e4, 1e5, 1e6]:
			if max_t > self.relax_max:
				break
			interval = round(max_t/cycles/delta_t)
			t_steps = cycles*interval
			quench_cond.append([delta_t, t_steps, interval, max_t, temp_scale])
		#
		conditions.append(quench_cond)
		#
		deform_cond = conditions[0]
		print('\n##################################')	
		print('##### Step Deform Conditions #####')
		print('##################################')
		print(f'{"Deform Mode: ":^14}', self.mode)
		print(f'{"Max Deformation":^14}', '=', self.deform_max)
		print(f'{"Deform Rate":^14}', '=', self.step_rate)
		print('#################################')
		print(f"{'Sim. Period':^14}", '=', deform_cond[0]*deform_cond[1])
		print(f"{'delta_T':^14}", '=', deform_cond[0])
		print(f"{'interval':^14}", '=', deform_cond[2])
		print(f"{'Records':^14}", '=', int(deform_cond[1]/deform_cond[2]))
		if deform_cond[4] == 1:
			print(f'{"# Temp. Scalling is ON! #":^30}')
		print('#################################')
		print()
		#
		quench_cond = conditions[1]
		print('#################################')
		print('#####   Quench Conditions   #####')
		print('#################################\n')
		for i, cond in enumerate(quench_cond):
			# acc_time += cond[0]*cond[2]
			print(f'#### Quench untill: {cond[3]: .1e} ####')
			# print(f"{'Quench Period':^14}", '=', cond[0]*cond[1])
			print(f"{'delta_T':^14}", '=', cond[0])
			print(f"{'interval':^14}", '=', cond[2])
			print(f"{'Records':^14}", '=', round(cond[1]/cond[2]))

			if cond[4] == 1:
				print(f'{"# Temp. Scalling is ON! #":^30}')
			print()
		#
		self.prompt()
		return conditions

	###############
	# 計算条件の確認
	def prompt(self):
		dic={'y':True,'yes':True,'n':False,'no':False}
		while True:
			inp = input(u'計算を続行 ==> [Y]es >> ').lower()
			if inp in dic:
				inp = dic[inp]
				break
			print(u'##### \nもう一度入力してください。')
		if inp :
			print(u"計算を続行します")
		else:
			sys.exit(u"##### \n計算を中止します。")
		return

	#################################################################################
	# 平衡計算したUDFファイルとその場所を選択
	def select_udf(self):
		t_udf, data_dir = self.file_select()
		#
		func, nu, structure = self.read_data(data_dir)
		#
		base = self.make_base_udf(t_udf)
		return base

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
		with open(t_data, 'r') as f:
			calc_cond = f.readlines()[-1].strip().split('\t')
		if len(calc_cond) == 3:
			func = 0
			nu = 0
			structure = 'none'
		elif len(calc_cond) == 6:
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
		# u.jump(u.totalRecord() - 1)
		# ave_xy = u.get('Statistics_Data.Stress.Total.Total_Average.xy')
		u.eraseRecord(0, u.totalRecord() - 2)
		u.write(base)
		return base

	#######################################################################
	# ファイル名を設定し、バッチファイルを作成
	def setup_batch(self, conditions, base):
		batch_series = ''
		for axis in ['x', 'y', 'z']:
			for rep in range(self.repeat):
				repeat = 'rotate_' + axis + "_" + str(rep)
				if platform.system() == "Windows":
					batch_series += 'cd /d %~dp0\\' + repeat +'\n'
					batch_series += 'call _step_deform.bat\n'
				elif platform.system() == "Linux":
					batch_series += 'cd ./' + repeat +'\n'
					batch_series += './_step_deform.bat\n'
					batch_series += 'cd ../\n'
				target_dir = self.set_dir(repeat)
				self.make_batch(base, axis, conditions, repeat, target_dir)
				# self.make_script_ss(target_dir)
				self.make_script(target_dir)
			if platform.system() == "Windows":
				batch_series += 'cd /d %~dp0\n'
			#
		batch_series += "python calc_ave.py\n"

		f_batch = os.path.join(self.calc_dir, '_calc_all.bat')
		with open(f_batch, 'w') as f:
			f.write(batch_series)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)
		return

	def set_dir(self, repeat):
		target_dir = os.path.join(self.calc_dir, repeat)
		if os.path.exists(target_dir):
			print("Use existing dir of ", target_dir)
		else:
			print("Make new dir of ", target_dir)
			os.makedirs(target_dir)
		return target_dir

	def make_batch(self, base, axis, conditions, rep, target_dir):
		if platform.system() == "Windows":
			batch = ""
		elif platform.system() == "Linux":
			batch = "#!/bin/bash\n\n"
		#
		deform_str = str(self.deform_max).replace('.', '_')
		rate_str = str(f"{self.step_rate:.1e}").replace('.', '_')
		if self.mode == 'elong':
			uin = 'Step_Elong_Lambda_' + deform_str + '_rate_' + rate_str + "_uin.udf"
		elif self.mode == 'shear':
			uin = 'Step_Shear_Gamma_' + deform_str + '_rate_' + rate_str + "_uin.udf"
		#
		batch = self.make_title(batch, "Calculating " + rep + " Step rate=" + rate_str)
		read_udf, batch = self.make_step(uin, batch)
		self.step_deform(base, axis, uin, conditions[0], target_dir)
		pre = read_udf
		template = uin
		# Plot Stress during deformation with Temperature
		# batch += 'python calc_ss.py\n'
		#
		for i, cond in enumerate(conditions[1]):
			present_udf = 'Quench_untill_' + str(f'{cond[3]:.0e}') + "_uin.udf"
			batch = self.make_title(batch, "Calculating " + rep + " Quench_untill_" + str(f'{cond[3]:.0e}'))
			read_udf, batch = self.make_step(present_udf, batch)
			self.eq_udf(template, pre, present_udf, cond, target_dir)
			pre = read_udf
			template = present_udf
		# Plot G(t) during Quench
		batch += 'python calc_gt_all.py\n'		

		# バッチファイルを作成
		f_batch = os.path.join(target_dir, '_step_deform.bat')
		with open(f_batch, 'w') as f:
			f.write(batch)
			if platform.system() == "Linux":
				os.chmod(f_batch, 0o777)
		return

	###########################
	# ターミナルのタイトルを設定
	def make_title(self, batch, title):
		if platform.system() == "Windows":
			batch += "title " + title + "\n"
		elif platform.system() == "Linux":
			batch += r'echo -ne "\033]0; ' + title + ' \007"' + '\n'
		return batch

	# ファイル名の処理
	def make_step(self, present_udf, batch):
		out_udf = present_udf.replace("uin", "out")
		batch += self.ver + ' -I ' + present_udf + ' -O ' + out_udf + self.core
		read_udf = out_udf
		return read_udf, batch
		
	#-----
	def step_deform(self, base, axis, udf_in, time_temp, target_dir):
		u = UDFManager(base)
		u.jump(-1)

		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(100000.,		p + 'Max_Force')
		u.put(time_temp[0],	p + 'Time.delta_T')
		u.put(time_temp[1],	p + 'Time.Total_Steps')
		u.put(time_temp[2],	p + 'Time.Output_Interval_Steps')
		#
		u.put(1.0,			p + 'Temperature.Temperature')
		u.put(time_temp[4], p + 'Temperature.Interval_of_Scale_Temp')
		#
		u.put(0,			p + 'Pressure_Stress.Pressure')

		# Deformation
		if self.mode == 'elong':
			p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
			u.put('Cell_Deformation', 		p + 'Method')
			u.put('Simple_Elongation', 		p + 'Cell_Deformation.Method')
			u.put('Initial_Strain_Rate', 	p + 'Cell_Deformation.Simple_Elongation.Input_Method')
			u.put(self.step_rate,	 		p + 'Cell_Deformation.Simple_Elongation.Initial_Strain_Rate.Rate')
			u.put(0.5, 						p + 'Cell_Deformation.Simple_Elongation.Poisson_Ratio')
			u.put('z', 						p + 'Cell_Deformation.Simple_Elongation.Axis')
			u.put(1, 						p + 'Cell_Deformation.Interval_of_Deform')
			u.put(0, 						p + 'Cell_Deformation.Deform_Atom')
		elif self.mode == 'shear':
			p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
			u.put('Lees_Edwards', 		p + 'Method')
			u.put('Steady', 		p + 'Lees_Edwards.Method')
			u.put(self.step_rate, 	p + 'Lees_Edwards.Steady.Shear_Rate')
		
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

		# Rotate system
		self.rotate_position(u, axis)

		#--- Write UDF ---
		u.write(os.path.join(target_dir, udf_in))
		return
	
	################################################################################
	def eq_udf(self, template, read_udf, present_udf, time_temp, target_dir):
		u = UDFManager(os.path.join(target_dir, template))
		u.jump(-1)

		#--- Simulation_Conditions ---
		# Dynamics_Conditions
		p = 'Simulation_Conditions.Dynamics_Conditions.'
		u.put(time_temp[0], p+'Time.delta_T')
		u.put(time_temp[1], p+'Time.Total_Steps')
		u.put(time_temp[2], p+'Time.Output_Interval_Steps')
		u.put(time_temp[4], p + 'Temperature.Interval_of_Scale_Temp')

		# Deformation
		p = "Simulation_Conditions.Dynamics_Conditions.Deformation."
		u.put('None',	p + 'Method')

		#--- Initial_Structure ---
		p = 'Initial_Structure.Read_Set_of_Molecules'
		u.put([read_udf, -1], p)
		p = 'Initial_Structure.Generate_Method.'
		u.put('Restart', 			p+'Method')
		u.put([read_udf, -1, 1, 1], p+'Restart')

		#--- Write UDF ---
		u.write(os.path.join(target_dir, present_udf))
		return

	###############################################################
	# アトムのポジションを回転
	def rotate_position(self, u, axis):
		R = self.rotate(axis, np.pi/2.)
		u.jump(u.totalRecord() - 1)
		pos = u.get('Structure.Position.mol[].atom[]')
		for i, mol in enumerate(pos):
			for j, atom in enumerate(mol):
				tmp = list(np.array(R).dot(np.array(atom)))
				u.put(tmp, 'Structure.Position.mol[].atom[]', [i, j])
		return

	def rotate(self, axis, deg):
		if axis == 'x':
			R = [
				[1., 0., 0.],
				[0., np.cos(deg), -1*np.sin(deg)],
				[0., np.sin(deg), np.cos(deg)]
			]
		elif axis == 'y':
			R = [
				[np.cos(deg), 0., np.sin(deg)],
				[0., 1., 0.],
				[-1*np.sin(deg), 0., np.cos(deg)]
			]
		elif axis == 'z':
			R = [
				[np.cos(deg), -1*np.sin(deg), 0.],
				[np.sin(deg), np.cos(deg), 0.],
				[0., 0., 1.]
			]
		return R


	#######################
	# 必要なスクリプトを作成
	def make_script(self, target_dir):
		script = self.script_content_gt_all()
		with open(os.path.join(target_dir, 'calc_gt_all.py'),'w') as f:
			f.write(script)
		#
		script = self.script_content_ave()
		with open(os.path.join(self.calc_dir, 'calc_ave.py'),'w') as f:
			f.write(script)
		#
		return

	# スクリプトの中身
	def script_content_gt_all(self):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n'
		script += '################################\n'
		script += 'import sys \n'
		script += 'py_mod = "' + self.py_mod + '"\n'
		script += 'sys.path.append(py_mod)\n'
		script += 'from network_evaluation import ReadStepStress as rs\n'
		script += '################################\n'
		script += 'mode = "' + self.mode + '"\n'
		script += 'stress_step, strain_step, gt_step = rs.calc_step(mode)\n'
		script += '#\n'
		script += 'rs.calc_quench(stress_step, strain_step, gt_step, mode)\n'
		#
		return script

	# スクリプトの中身
	def script_content_ave(self):
		script = '#!/usr/bin/env python \n# -*- coding: utf-8 -*-\n'
		script += '################################\n'
		script += 'import sys \n'
		script += 'py_mod = "' + self.py_mod + '"\n'
		script += 'sys.path.append(py_mod)\n'
		script += 'from network_evaluation import ReadStepStress as rs\n'
		script += '################################\n'
		script += 'gt_ave_all, strain_ave = rs.average()\n'
		script += 'rs.irheo(gt_ave_all)'
		#
		return script

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()
