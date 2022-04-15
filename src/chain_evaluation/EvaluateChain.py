#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
# Import Modules
################################################################################
from UDFManager import UDFManager
import os
import sys
import math
import cmath
import numpy as np
import platform
import subprocess
import scipy.signal as signal
#
import CognacUtility as CU
from CognacBasicAnalysis import *
from CognacGeometryAnalysis import CognacGeometryAnalysis
################################################################################
################################################################################
def evaluate_nw():
	target = file_select()
	#
	calc_cond, chain_list = make_chain_list(target)
	# ポリマー鎖関連の特性情報を計算
	ec = EvaluateChain(calc_cond, chain_list, target)
	ec.eval_chain()
	return

##############################
# 対象となる udf ファイルを選択
def file_select():
	param = sys.argv
	if len(param) == 1:
		print("usage: python", param[0], "Honya_out.udf")
		exit(1)
	elif not os.access(param[1],os.R_OK):
		print(param[1], "not exists.")
		exit(1)
	else:
		target = param[1]
	return target

# 計算条件から、ホモポリマーとネットワークを判断し、chain_list を読み出す。
def make_chain_list(target):
	# 計算対象の条件を読み取る
	if not os.access('target_condition.udf', os.R_OK):
		print("'target_condition.udf' is not exists.")
		exit(1)
	else:
		cond_u = UDFManager('target_condition.udf')
		nw_type = cond_u.get('TargetCond.Model.TargetModel')
		func = cond_u.get('TargetCond.NetWork.N_Strands')
		n_seg = cond_u.get('TargetCond.NetWork.N_Segments')
		l_bond= cond_u.get('SimulationCond.l_bond')
		cn = cond_u.get('TargetCond.Strand.Characteristic_Ratio')
		nu = cond_u.get('TargetCond.System.Nu')
		calc_cond = [n_seg, l_bond, cn, func, nu, nw_type]
	#
	sc = SelectChain(target)
	if len(calc_cond) == 3:
		# ホモポリマーのリストを作成
		end_list, chain_list = sc.make_chain_list()
	elif len(calc_cond) == 6:
		# ネットワークストランドのリストを作成
		# ss = Init_Strand_Select(target)
		jp_list, jp_pair_list, chain_list = sc.make_strand_list()
	return calc_cond, chain_list

##############################################################################
class SelectChain:
	def __init__(self, target):
		self.uobj = UDFManager(target)

	# ##########################
	# # ホモポリマーのリストを作成
	# def make_chain_list(self):
	# 	atom_list = self.uobj.get("Set_of_Molecules.molecule[].atom[]")
	# 	self.uobj.jump(-1)
	# 	chain_list = []
	# 	end_list = []
	# 	for i, target_chain in enumerate(atom_list):
	# 		tmp = []
	# 		for j, atom in enumerate(target_chain):
	# 			tmp.append(j)
	# 		end_list.append([i, [tmp[0], tmp[-1]]])
	# 		chain_list.append([i, tmp])
	# 	return end_list, chain_list

	########################################
	# 架橋点およびストランドの構成アトムのリスト
	def make_strand_list(self):
		jp_list = self.make_jp_list()
		#
		jp_pair_list = []
		strand_list = []
		for target_jp in jp_list:
			jp_pair, strand = self.make_jp_pair(target_jp, jp_list)
			for i in jp_pair:
				jp_pair_list.append(i)
			if len(strand) > 0:
				for i in strand:
					strand_list.append(i)
		return jp_list, jp_pair_list, strand_list

	# 架橋点のリストを作成
	def make_jp_list(self):
		self.uobj.jump(-1)
		jp_list = []
		#
		mols = self.uobj.get("Set_of_Molecules.molecule[]")
		for i, mol in enumerate(mols):
			for j, atom in enumerate(mol[1]):
				tmp = []
				if atom[1] == 'JP_A' or atom[1] == 'JP_B':
					jp_list.append([i, j])
		return jp_list

	# 架橋点どうしのペアを作成
	def make_jp_pair(self, target_jp, jp_list):
		molecule = target_jp[0]
		start_jp = target_jp[1]
		jp_pair = []
		strand = []
		bonds = self.uobj.get("Set_of_Molecules.molecule[].bond[]")
		tmp_bonds = bonds[molecule]
		#
		for i, bond in enumerate(tmp_bonds):
			tmp = []
			if ((bond[1] == start_jp) or (bond[2] == start_jp)) and (i < len(tmp_bonds) - 1):
				if bond[1] == start_jp:
					adj = bond[2]
				else:
					adj = bond[1]
				tmp.append(start_jp)
				tmp.append(adj)
				tmp_id = i + 1
				while tmp_bonds[tmp_id][0] == "bond_Strand":
					if tmp_bonds[tmp_id][1] == adj:
						adj = tmp_bonds[tmp_id][2]
					elif tmp_bonds[tmp_id][2] == adj:
						adj = tmp_bonds[tmp_id][1]
					tmp.append(adj)
					tmp_id += 1
				#
				if tmp_bonds[tmp_id][1] == adj:
					end_jp = tmp_bonds[tmp_id][2]
				elif tmp_bonds[tmp_id][2] == adj:
					end_jp = tmp_bonds[tmp_id][1]
				if len(tmp)>2:
					tmp.append(end_jp)
					jp_pair.append([molecule, [start_jp, end_jp]])
					strand.append([molecule, tmp])
		return jp_pair, strand


###############################################################################
###############################################################################
class EvaluateChain:
	def __init__(self, calc_cond, chain_list, target):
		self.calc_cond = calc_cond
		self.n_seg = int(calc_cond[0])
		self.l_bond = float(calc_cond[1])
		#
		self.chain_list = chain_list
		#
		self.target = target
		self.target_name = target.split('.')[0]
		self.uobj = UDFManager(target)

	def eval_chain(self):
		bond_list = []
		angle_list = []
		Rx_list = []
		Ry_list = []
		Rz_list = []
		R_list = []
		#
		gr_list = []
		cn_list = []
		#
		if self.target.split('_')[0] == 'GK':
			self.calc_gk()
			#
			# gt2 = modify(data_list)

			#
			return
		else:
			rec_size = self.uobj.totalRecord()
			for rec in range(1, rec_size):
				print("Reading Rec=", rec, '/', rec_size)
				bond, angle, e2e_x, e2e_y, e2e_z, e2e, r2, gr, cn = self.read_chain(rec)
				bond_list.extend(bond)
				angle_list.extend(angle)
				Rx_list.extend(e2e_x)
				Ry_list.extend(e2e_y)
				Rz_list.extend(e2e_z)
				R_list.extend(e2e)
				#
				gr_list.append(gr)
				cn_list.append(cn)
			# 鎖に沿ったセグメント間距離の平均を計算
			cn_ave, cn_part = self.calc_cn(cn_list)
			#
			self.make_output(bond_list, angle_list, Rx_list, Ry_list, Rz_list, R_list, gr_list, cn_list, cn_ave, cn_part)
			return

	####################################################################################
	# Green Kubo での計算を処理
	def calc_gk(self):
		corr, corr_all = self.calc_corr()
		mm = MakeMulti(["Corr_stress", corr, ['Time', 'sigma']], self.target_name)
		mm.make_all()
		mm = MakeMulti(["Corr_stress_all", corr_all, ['Time', 'sigma', 'ave']], self.target_name)
		mm.make_all()
		# mm = MakeMulti(["Corr_stress_semi", corr, ['Time', 'sigma']], self.target_name)
		# mm.make_all()
		self.irheo(corr)
		# mm = MakeMulti(["Corr_stress_mod", corr_mod, ['Time', 'sigma']], self.target_name)
		# mm.make_all()
		return

	def calc_corr(self):
		self.uobj.jump(self.uobj.totalRecord() -1)
		#
		vol = self.uobj.get('Statistics_Data.Volume.Total_Average')
		corr_all = self.uobj.get('Correlation_Functions.Stress.Correlation[]')
		corr = []
		prev = 0.
		for data in corr_all:
			time = data[0]
			ave = vol*np.average(np.array(data[2:]))
			if data[1] > 0:
				g = data[1]
				prev = data[1]
			else:
				g = prev
		# g_mod = signal.savgol_filter(g, 3, 2)
			corr.append([time, g, ave])
		# corr_mod = np.stack([time, g_mod], 1)
		return corr, corr_all


	##################################
	# 
	def irheo(self, data_list):
		minmax = [1e-5, 1e2]
		div = 10
		#
		# mod = self.modify(data_list)
		# gt, mod_gt = self.modify_data(mod)
		#
		gw = self.calcgw(data_list, minmax, div)
		self.save_data(gw, 'gw.dat')
		#
		# self.save_data(data_list, 'modified.dat')
		# self.plotgtgw('modified.dat')
		# cmd = "corr2gw < modified.dat > gw.dat"
		# subprocess.call(cmd, shell=True)

		self.plotgtgw('gw.dat')
		#
		return

	def modify_data(self, data_list):
		fine_div = 100
		#
		glist = []
		timelist = []
		for data in data_list:
			time = data[0]
			g = data[1]
			if time == 0.0:
				timelist.append(time)
				glist.append(g)
			else:
				for i in range(1, fine_div + 1):
					timelist.append(pre_time + i*(time-pre_time)/fine_div)
					glist.append(pre_g + i*(g-pre_g)/fine_div)
			pre_time = time
			pre_g = g
			#
		mod_g = signal.savgol_filter(glist, 5, 3)
		#
		gt = np.stack([timelist, glist], 1)
		mod_gt = np.stack([timelist, mod_g], 1)
		#
		return gt, mod_gt

	def calcgw(self, gt, minmax, div):
		gw = []
		mag = math.log10(minmax[0])
		while mag < math.log10(minmax[1]):
			for i in range(div):
				omega = 10**(mag+i/div)
				gstar = self.gs(gt, omega)
				gw.append([omega, gstar.real, abs(gstar.imag)])
			mag += 1
		#
		return gw

	def gs(self, gt, omega):
		gstar = gt[0][1] + (1 - cmath.exp(-1j*omega*gt[1][0]))*(gt[1][1] - gt[0][1])/gt[1][0]/(1j*omega)
		for k in range(len(gt) - 2):
			gstar += (gt[k+2][1] - gt[k+1][1])*(cmath.exp(-1j*omega*gt[k+1][0]) - cmath.exp(-1j*omega*gt[k+2][0]))/(gt[k+2][0] - gt[k+1][0])/(1j*omega)
		#
		return gstar 

	#----- 計算結果をターゲットファイル名で保存
	def save_data(self, target, f_data):
		with open(f_data,'w') as f:
			for line in target:
				for data in line:
					f.write(str(data) + '\t')
				f.write('\n')
		return

	#----- 結果をプロット
	def plotgtgw(self, f_data):
		plt = self.make_gtgw(f_data)
		#
		if platform.system() == "Windows":
			subprocess.call([plt], shell=True)
		elif platform.system() == "Linux":
			subprocess.call(['gnuplot ' + plt], shell=True)
		return

	# 必要なスクリプトを作成
	def make_gtgw(self, f_data):
		script = self.gtgw_content(f_data)
		plt = f_data.replace('dat', 'plt')
		with open(plt, 'w') as f:
			f.write(script)
		return plt

	# スクリプトの中身
	def gtgw_content(self, f_data):
		out_png = f_data.replace('dat', 'png')
		script = 'set term pngcairo font "Arial,14"\n\n'
		script += 'set colorsequence classic\n\n'
		script += 'data = "' + f_data + '"\n\n'
		script += 'set output "' + out_png + '"\n\n'
		script += 'set key left\nset size square\n'
		script += '#set xrange [1:4]\n#set yrange [0:0.2]\n#set xtics 1\n#set ytics 0.1\n'

		if f_data == 'modified.dat' or f_data == 'ave_all_stress.dat':
			script += 'set logscale xy\n'
			script += 'set format x "10^{%L}" \nset format y "10^{%L}"\n'
			script += 'set xlabel "Time"\nset ylabel "Stress"\n'	
			script += 'plot	data u 1:2 axis x1y1 w l lw 2 lt 1 ti "Stress"'
		elif f_data == 'gw.dat' or f_data == 'freq_mod.dat':
			script += 'set xrange [:1e2]\nset yrange [1e-4:]\nset y2range [1e-1:1e1]\nset y2tics\n'
			script += 'set logscale xyy2\n'
			script += '# 斜辺の傾きが -2 の三角形の準備\n'
			script += 'a = 30; # グラフの中に入るように三角形の高さを調整\n'
			script += 'x1=5e-4; x2=1e-3;\n'
			script += 'y1=a*x1**(1);y2=a*x2**(1);\n'
			script += 'set object 1 polygon from x1,y1 to x2,y1 to x2,y2 to x1,y1 fs empty border\n\n'
			script += 'set format x "10^{%L}" \nset format y "10^{%L}"\nset format y2 "10^{%L}"\n'
			# script += 'set label 1 sprintf("{/Symbol l} = %.1f", deform) at graph 0.6, 0.9\n\n'
			script += 'set xlabel "Frequency"\nset ylabel "G' + "', G''" + '"\nset y2label "tan{/Symbol d}"\n\n'
			script += 'plot	'
			script += 'data u 1:2 w lp lt 1 ti "G' + "'" + '", \\\n'
			script += 'data u 1:3 w lp lt 2 ti "G' + "''" + '", \\\n'
			script += 'data u 1:($3/$2) axis x1y2 w lp lt 3 ti "tan{/Symbol d}"'
		script += '\n\nreset'

		return script

	def modify(self, data_list):
		a = 0.057
		tau = 190
		fitstart = 500
		mod_gt = []
		for data in data_list:
			time = float(data[0])
			g = float(data[1])
			if time < fitstart:
				# if g > 0:
				mod_gt.append([time, g])
			else:
				break
		time = fitstart
		while time < 1e5:
			tmp = a*np.exp(-time/tau)
			if tmp > 1e-10:
				mod_gt.append([time, tmp])
				time += 10**int(np.log10(time))/100
			else:
				break
		# save_data(mod_gt, 'mod_gt.dat')
		return mod_gt


	###################################################################
	# 鎖に沿ったセグメント間距離の平均を計算
	def calc_cn(self, cn_list):
		cn_ave = []
		cn_part = []
		#
		l_part = len(cn_list)//10
		# データの分割
		multi = 0
		part_cn = []
		tmp = []
		for i, part in enumerate(cn_list):
			if i < l_part*(multi + 1):
				tmp.append(part)
			else:
				part_cn.append(tmp)
				tmp = []
				tmp.append(part)
				multi += 1
		# 各パートごとに平均
		for part in part_cn:
			tmp = [ [i + 1, 0] for i in range(len(cn_list[0]))]
			count = 0
			cn_part_ave = []
			for data in part:
				for i, el in enumerate(data):
					tmp[i][1] += el[1]
				count += 1
			for data in tmp:
				cn_part_ave.append([data[0], data[1]/count])
			cn_part.append(cn_part_ave)
		# パートごとの平均をさらに平均
		tmp = [ [i + 1, 0] for i in range(len(cn_list[0]))]
		count = 0
		for data in cn_part:
			for i, el in enumerate(data):
				tmp[i][1] += el[1]
			count += 1
		for data in tmp:
			cn_ave.append([data[0], data[1]/count])
		return cn_ave, cn_part

	def make_output(self, bond_list, angle_list, Rx_list, Ry_list, Rz_list, R_list, gr_list, cn_list, cn_ave, cn_part):
		# 結果をヒストグラムで出力 
		hist_list = [
				["bond", bond_list, 200, "True", ['bond length', 'Freq.'], 'box'],
				["angle", angle_list, 200, "True", ['angle [deg]', 'Freq.'], 'box'],
				["Rx", Rx_list, 200, "True", ['|Rx|', 'Freq.'], self.calc_cond ],
				["Ry", Ry_list, 200, "True", ['|Ry|', 'Freq.'], self.calc_cond ],
				["Rz", Rz_list, 200, "True", ['|Rz|', 'Freq.'], self.calc_cond ],
				["R", R_list, 200, "True", ['|R|', 'Freq.'], self.calc_cond ]
				]
		for cond in hist_list:
			mh = MakeHist(cond, self.target_name)
			mh.make_hist_all()

		# マルチ形式での出力
		multi_list = [
				["gr", gr_list, ['Distance', 'g(r)']],
				["CN", cn_list, ['|i-j|', 'C_{|i-j|}']],
				["CN_part", cn_part, ['|i-j|', 'C_{|i-j|}']],
				["CN_ave", cn_ave, ['|i-j|', 'C_{|i-j|}']]
				]
		for cond in multi_list:
			mm = MakeMulti(cond, self.target_name)
			mm.make_all()
		return

	# ポリマー鎖関連の特性情報
	def read_chain(self, rec):
		# 初期化
		self.uobj.jump(rec)
		self.bound_setup()
		CU.setCell(tuple(self.uobj.get("Structure.Unit_Cell.Cell_Size")))
		# ステップの数に対応した空リストを作成
		r2_ij = [[] for i in range(len(self.chain_list[0][1]))]
		#
		e2e_x = []
		e2e_y = []
		e2e_z = []
		e2e_list = []
		r2_list = []
		bond_list = []
		cn = []
		#
		xp = [[] for i in range(len(self.chain_list[0][1]))]
		# 
		ba = CognacBasicAnalysis(self.target, rec)
		for chain in self.chain_list:
			mol = chain[0]
			c_len = len(chain[1])
			atom = self.uobj.get("Set_of_Molecules.molecule[].atom[]", [mol, chain[1][2]])[1]
			#		
			for step in range(1, c_len):
				for start in range(c_len - step):
					if len(self.calc_cond) == 3: # ポリマー鎖の場合
						e2e_vec = ba.vector([mol, chain[1][start]], [mol, chain[1][start + step]])
					elif len(self.calc_cond) == 6: # ストランドの場合
						end1 = tuple(self.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start]]))
						end2 = tuple(self.uobj.get("Structure.Position.mol[].atom[]", [mol, chain[1][start + step]]))
						e2e_vec = CU.distanceWithBoundary(end1, end2)
					e2e_dist = np.linalg.norm(np.array(e2e_vec))
					r2 = e2e_dist**2
					r2_ij[step].append(r2)
					if step == 1:
						bond_list.append(e2e_dist)
					if step == c_len -1:
						e2e_x.append(e2e_vec[0])
						e2e_y.append(e2e_vec[1])
						e2e_z.append(e2e_vec[2])
						#
						e2e_list.append(e2e_dist)
						r2_list.append(r2)
			#
			# for p in range(c_len):
			# 	xp[p].append(np.linalg.norm(ba.Xp(mol, p)))
		#
		# xp_list = []
		# for i in range(c_len):
		# 	xp_list.append([i+1, np.average(np.array(xp[i]))])
		# print(xp_list)

		# gr
		cg = CognacGeometryAnalysis(self.target, rec)
		gr = cg.gr([atom])
		# cn
		for i in range(1, len(r2_ij)):
			cn.append([i, np.average(np.array(r2_ij[i]))/(i*self.l_bond**2)])
		# angle
		anglename = self.uobj.get("Molecular_Attributes.Angle_Potential[].Name")
		tmp = np.array(ba.angle(anglename[0]))
		angle_list = list(tmp[~np.isnan(tmp)])
		# print(cn[:3])
		# print(r2_list[:3])
		return bond_list, angle_list, e2e_x, e2e_y, e2e_z, e2e_list, r2_list, gr, cn
	
	# 周期境界条件の設定
	def bound_setup(self):
		axis = self.uobj.get("Simulation_Conditions.Boundary_Conditions")
		boundarylist = [0,0,0]
		#
		for i in range(0,3):
			if axis[i] == "NONE" :
				boundarylist[i] = 0
			elif axis[i] == "PERIODIC" :
				boundarylist[i] = 1
			elif axis[i] == "REFLECTIVE1" :
				boundarylist[i] = 2
			elif axis[i] == "REFLECTIVE2" :
				boundarylist[i] = 3
		CU.setBoundary(tuple(boundarylist))
		return


##############################################################################################
##############################################################################################
class MakeHist:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, bins, normalize, Legend, option]
		self.list = cond_list[1]
		self.bins = cond_list[2]

		self.dir = os.path.join(target_name, cond_list[0])

		self.base = cond_list[0]
		self.norm = cond_list[3]
		#
		self.f_dat = cond_list[0] + "_hist.dat"
		self.f_plt = cond_list[0] + ".plt"
		self.f_png = cond_list[0] + ".png"
		self.leg = cond_list[4]
		self.option = cond_list[5]

	# ヒストグラムのグラフの作成
	def make_hist_all(self):
		# ヒストグラムのデータ作成
		bin_width, hist_data = self.make_hist_data()
		# ヒストグラムのデータを書き出し 
		self.write_data(hist_data, bin_width)
		# グラフを作成
		self.make_graph(bin_width)
		return

	# ヒストグラムのデータ作成
	def make_hist_data(self):
		# ヒストグラムを作成
		weights = np.ones(len(self.list))/float(len(self.list))
		if self.norm:
			val, x = np.histogram(self.list, bins=self.bins, weights= weights)
		else:
			val, x = np.histogram(self.list, bins=self.bins)
		# グラフ用にデータを変更
		bin_width = (x[1]-x[0])
		mod_x = (x + bin_width/2)[:-1]
		hist_data = np.stack([mod_x, val], axis = 1)
		return bin_width, hist_data

	# ヒストグラムのデータを書き出し 
	def write_data(self, hist_data, bin_width):
		os.makedirs(self.dir, exist_ok=True)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# Histgram data:\n\n")
			for line in hist_data:
				f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

	# グラフを作成
	def make_graph(self, bin_width):
		self.make_script(bin_width)
		cwd = os.getcwd()
		os.chdir(self.dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return

	# 必要なスクリプトを作成
	def make_script(self, bin_width):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content(bin_width)
			f.write(script)
		return

	# スクリプトの中身
	def script_content(self, bin_width):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		#
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		#
		script += '#\nset size square\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		#
		if self.base == "Rx" or self.base == "Ry" or self.base == "Rz":
			if type(self.option) == list:
				n_seg = self.option[0]
				bond = self.option[1]
				cn = self.option[2]
				func = self.option[3]
				nu = self.option[4]
				nw_type = self.option[5]
			script += 'N = ' + str(n_seg) + '\n'
			script += 'bond = ' + str(bond) + '\n'
			script += 'CN = ' + str(cn) + '\n'
			script += 'func = ' + str(func) + '\n\n'
			#
			script += 'R1 = bond*(CN*(N+1))**0.5\n'
			script += 'C=0.1\n'
			script += 'delta=R1\n'
			script += 'frc = 1.0\n\n'
			#
			if nw_type == 'Regular' and func == 3:
				script += 'set xrange [0:]\n#set yrange [0:100]\n'
				script += 'Pos = R1/2**0.5\n\n'
				script += 'f(x) = C*(1./2.)*(1./(frc*delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*(frc*delta)**2)) + exp(-1.*((x+Pos)**2)/(2.*(frc*delta)**2)))\n\n'
				script += 'fit f(x) data via C, frc\n\n'
				script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
			elif nw_type == 'Regular' and func == 4:
				script += 'set xrange [0:]\n#set yrange [0:100]\n'
				script += 'Pos = R1/3**0.5\ndelta = Pos*(1. - 2./func)**0.5\n\n'
				script += 'f(x) = C*(1./2.)*(1./(frc*delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*(frc*delta)**2)) + exp(-1.*((x+Pos)**2)/(2.*(frc*delta)**2)))\n\n'
				script += 'fit f(x) data via C, frc\n\n'
				script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
			elif nw_type == 'Regular' and func == 6:
				script += 'set xrange [0:]\n#set yrange [0:100]\n'
				script += 'Pos = R1\ndelta = Pos*(1. - 2./func)**0.5\n\n'
				script += 'f(x) = C*(1./2.)*(1./(frc*delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*(frc*delta)**2)) + exp(-1.*((x+Pos)**2)/(2.*(frc*delta)**2)))\n\n'
				script += 'fit f(x) data via C, frc\n\n'
				script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
			elif nw_type == 'Regular' and func == 8:
				script += 'set xrange [0:]\n#set yrange [0:100]\n'
				script += 'Pos = R1/3**0.5\ndelta = Pos*(1. - 2./func)**0.5\n\n'
				script += 'f(x) = C*(1./2.)*(1./(frc*delta*(3.142*2.)**0.5))*(exp(-1.*((x-Pos)**2)/(2.*(frc*delta)**2)) + exp(-1.*((x+Pos)**2)/(2.*(frc*delta)**2)))\n\n'
				script += 'fit f(x) data via C, frc\n\n'
				script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
			else:
				script += '#set xrange [0:]\n#set yrange [0:100]\nfrc=1.0\n\n'
				script += 'f(x) = C*exp(-1.*x**2./(2.*(frc*R1)**2.))/(2.*pi*(frc*R1)**2.)**(1/2)\n\n'
				script += 'fit f(x) data via C, frc\n\n'
				script += '#\nset label 1 sprintf("frc =%.3f", frc) at graph 0.7, 0.8\n\n'
			#
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			script += ', \\\n f(x)'

		if self.base == "R":
			if type(self.option) == list:
				n_seg = self.option[0]
				bond = self.option[1]
				cn = self.option[2]
				func = self.option[3]
			script += 'N = ' + str(n_seg) + '\n'
			script += 'bond = ' + str(bond) + '\n'
			script += 'CN = ' + str(cn) + '\n'
			script += 'f = ' + str(func) + '\n'
			script += 'R1 = bond*(CN*(N+1))**0.5\n'
			script += 'C=0.1\nfrc=1.0\n\n'
			script += 'f(x) = C*exp(-1.*(x-frc*R1)**2./(2.*sigma**2.))/(2.*pi*sigma**2.)**(1/2)\n\n'
			script += '#f(x) = C*4.*pi*x**2.*(3./(2.*pi*(frc*R1)**2.))**(3./2.)*exp(-3.*x**2./(2.*(frc*R1)**2.))\n'	
			script += 'fit f(x) data via frc, C, sigma\n\n'
			script += '#\nset label 1 sprintf("frc.=%.3f", frc) at graph 0.7, 0.8\n'
			script += 'set label 1 sprintf("sigma=%.3f", sigma) at graph 0.7, 0.8\n\n'
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			script += ', \\\n f(x)'
		#
		if self.base == "angle":
			if self.option != "box":
				script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w l noti'
			else:
				script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
				script += 'plot data u 1:($2/(3.142*sin(3.142*$1/180))) w boxes noti'
		#
		if self.base == "bond":
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
		#
		elif self.option == "box":
			script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
			script += '#\nplot data w boxes noti'
			
		return script


#############################################################################################
class MakeMulti:
	def __init__(self, cond_list, target_name):
		# cond_list = [base_name, data_list, Legend]
		self.list = cond_list[1]

		self.dir = os.path.join(target_name, cond_list[0])

		self.base = cond_list[0]
		self.repeat = len(cond_list[1])
		#
		self.f_dat = cond_list[0] + ".dat"
		self.f_plt = cond_list[0] + ".plt"
		self.f_png = cond_list[0] + ".png"
		self.leg = cond_list[2]
	
	############################################################
	# マルチリストのグラフの作成
	def make_all(self):
		# データを書き出し 
		self.write_data()
		# グラフを作成
		self.make_graph()
		return

	##############################
	# データを書き出し 
	def write_data(self):
		os.makedirs(self.dir, exist_ok=True)
		with open(os.path.join(self.dir, self.f_dat), 'w') as f:
			f.write("# data:\n")
			if self.base == 'CN_ave' or self.base == 'Corr_stress' or self.base == 'Corr_stress_semi' or self.base == 'Corr_stress_mod' or self.base == 'Corr_stress_all':
				for line in self.list:
					for data in line:
						f.write(str(data) + '\t')
					f.write('\n')
			else:
				for i, data in enumerate(self.list):
					f.write("\n\n# " + str(i) +":\n\n")
					for line in data:
						# print(line)
						f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

	###########################################################
	# グラフを作成
	def make_graph(self):
		self.make_script()
		cwd = os.getcwd()
		os.chdir(self.dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return

	#######################
	# 必要なスクリプトを作成
	def make_script(self):
		with open(os.path.join(self.dir, self.f_plt), 'w') as f:
			script = self.script_content()
			f.write(script)
		return

	#################
	# スクリプトの中身
	def script_content(self):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		#
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		#
		script += '#\nset size square\n#set xrange [1:]\n#set yrange [1:]\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		#
		if self.base == "CN" or self.base == "CN_ave" or self.base == "CN_part":
			script += '#\nset xrange [1:]\nset yrange [1:]\n'
			script += 'set key bottom\n\n'
			script += 'ct = 0.274\n'
			script += "f(x) = (1+ct)/(1-ct) -(2*ct*(1-ct**x))/(1-ct)**2/x\n\n"
			script += 'plot '
			if self.base == "CN":
				for i in range(self.repeat):
					script += 'data ind ' + str(i) + ' w l lc ' + str(i) + ' noti, \\\n'
			elif self.base == "CN_part":
				for i in range(self.repeat):
					script += 'data ind ' + str(i) + ' w l lc ' + str(i) + ' ti "part:' + str(i) + '", \\\n'
			else:
				script += 'data w l ti "averaged", \\\n'
			script += 'f(x) w l lw 2 ti "FreeRotationalModel"'
		elif self.base == 'Corr_stress' or self.base == 'Corr_stress_mod':
			script += 'set logscale xy \n\nset format x "10^{%L}" \nset format y "10^{%L}"\n\n'
			script += 'plot '
			script += 'data w l ti "Stress" \\\n'
		elif self.base == 'Corr_stress_semi':
			script += 'set logscale y \n\n#set format x "10^{%L}" \nset format y "10^{%L}"\n\n'
			script += 'a = 1\ntau =1000\n\ns = 100\ne = 1000\n\n'
			script += 'f(x) = a*exp(-1*x/tau) \n'
			script += 'fit [s:e] f(x) data usi 1:2 via a,tau\n\n'
			script += 'set label 1 sprintf("Fitted \\nA = %.1e \\n{/Symbol t} = %.1e \\nFitting Region: %d to %d", a, tau, s, e) at graph 0.35, 0.75\n\n'
			script += 'plot '
			script += 'data w l ti "Stress", \\\n'
			script += '[s:e] f(x) noti'
		elif self.base == 'Corr_stress_all':
			script += 'set logscale xy \n\nset format x "10^{%L}" \nset format y "10^{%L}"\n\n'
			script += 'plot data u 1:2 w l ti "G_t", \\\n'
			script += 'data u 1:3 w l ti "xy", \\\n'
			script += 'data u 1:4 w l ti "yz", \\\n'
			script += 'data u 1:5 w l ti "zx", \\\n'
			script += 'data u 1:6 w l ti "xx-yy", \\\n'
			script += 'data u 1:7 w l ti "yy-zz"'
		else:
			script += 'plot '
			for i in range(self.repeat):
				script += 'data ind ' + str(i) + ' w l lc ' + str(i) + 'noti, \\\n'

		return script




#######################################################

