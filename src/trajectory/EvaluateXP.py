#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
# Import Modules
################################################################################
import os
import sys
import math
import cmath
import numpy as np
import platform
import subprocess
import scipy.signal as signal
#
from UDFManager import UDFManager
import CognacUtility as cu
from CognacBasicAnalysis import *
from CognacGeometryAnalysis import CognacGeometryAnalysis
from trajectory.XpCalc import XpCalc
#
import trajectory.variables as var
################################################################################
# MAIN
################################################################################
def calc_xp():
	# 対象となる udf ファイルを選択
	select_udf()
	#
	evaluate_xp()
	# 計算結果を出力
	make_output()
	return

##########################################
# 対象となる udf ファイルを選択
def select_udf():
	param = sys.argv
	if len(param) == 1:
		print("usage: python", param[0], "Honya_out.udf")
		exit(1)
	elif not os.access(param[1],os.R_OK):
		print(param[1], "not exists.")
		exit(1)
	var.target = param[1]
	var.target_name = var.target.split('.')[0]
	var.uobj = UDFManager(var.target)
	return

###############################################################################
# ポリマー鎖関連の特性情報を計算
###############################################################################
def evaluate_xp():
	molname = "polymerA"
	time()
	make_chainList(molname, 'Xp')
	p = 1
	cp = normalCoordinate(p)
	ndata = len(cp)

	for i in range(0, ndata):
		var.xp_data.append([cp[i][0],cp[i][1],cp[i][2][0],cp[i][2][1],cp[i][2][2] ])
	return


#----- utility functions -----
def time():
	startTime = 0.0        # output UDF always starts at time=0
	var.timeRecord = [0]
	var.timeList = [startTime]
	for i in range(0, var.uobj.totalRecord()):
		var.uobj.jump(i)
		time = var.uobj.get("Time")
		if time == startTime:
			var.timeRecord[0] = i
		else:
			var.timeRecord.append(i)
			var.timeList.append(time - startTime)
	return

def make_chainList(molname, suffix=None):
	'''list of molecules having the specifiled molname

	Args:
		molname (str): name of the molecules to be selected
		suffix (str): suffix to be added to the 'key'

	Returns:
		list [(molIndex, key), ...]
		where key = 'molname:molIndex:suffix'
	'''
	key = '{}:{}'
	if suffix != None:
		key += ':' + suffix
	
	for m, atoms in enumerate(var.uobj.get('Set_of_Molecules.molecule[]')):
		tmp = []
		if var.uobj.get("Set_of_Molecules.molecule[].Mol_Name",[m]) == molname:
			for i in range(len(atoms[1])):
				tmp.append(i)
			var.chainList.append((m, key.format(molname, m)))
			var.chainAtom.append([m, tmp])
	if len(var.chainList) == 0:
		print('no molecule with Mol_Name: ' + molname)
	return


def normalCoordinate(p=1):
	cu.clearVectorMap()
	for rec in var.timeRecord:
		var.uobj.jump(rec)
		pos_list = tuple(var.uobj.get("Structure.Position.mol[].atom[]"))
		chainatom_pos = make_chains(pos_list)
		for m, key in var.chainList:
			# pos = np.array(pos_list[m])
			pos = np.array(chainatom_pos[m])
			cu.pushVector(key, Xp(pos, p))

	CpAve = cu.vectorCorrelation()
	S = np.sum(CpAve[0])
	x0, y0, z0 = CpAve[0]
	results = [ (var.timeList[i], np.sum(CpAve[i])/S, (CpAve[i][0]/x0, CpAve[i][1]/y0, CpAve[i][2]/z0)) 
					for i in range(len(CpAve)) ]
	return results

def make_chains(pos_list):
	chainatom_pos = []
	for chain in var.chainAtom:
		mol = chain[0]
		pos = []
		for atom in range(len(var.chainAtom[mol][1])):
			segment = pos_list[mol][atom]
			pos.append(segment)
		chainatom_pos.append(pos)
	return chainatom_pos

def Xp(pos, p=1):
	N = len(pos)
	k = np.pi*p/N
	xp = np.zeros(3)
	for n in range(N):
		xp += np.cos(k*(n+0.5))*pos[n]
	return tuple(xp/N)



###############################################################################
# 計算結果を出力
###############################################################################
def make_output():
	# マルチ形式での出力
	multi_list = [
			["Xp", var.xp_data, ['time', 'AutoCorr.']],
			]
	for cond in multi_list:
		make_multi(cond)
	return








def atomMsd(self, molname, atomIndex=None,
			start_record="first", end_record="last", cancel_trans=False):
	m_aList = self._mol_atomList(molname, atomIndex)
	if len(m_aList) == 0:
		return 0

	cu.clearVectorMap()
	for rec in self._recList(start_record, end_record):
		self.jump(rec)
		if cancel_trans:
			Gx, Gy, Gz = self.system_centerofmass()
		for m, a, key in m_aList:
			x, y, z = self.position([m,a])
			if cancel_trans:
				x -= Gx
				y -= Gy
				z -= Gz
			cu.pushVector(key, (x,y,z))

	msdAve = cu.msd()

	return [ (self.timeList[i], np.sum(msdAve[i]), tuple(msdAve[i]))
				for i in range(1,len(msdAve)) ]



def _atomList(self, atomname):
	'''list of atoms having the specified atomname

	Args:
		atomname (str): name of the atoms to be selected

	Returns:
		list [ (molIndex, atomIndex, key), ...],
		where key = 'atomname:molIndex:atomIndex'
	'''
	udfpath = 'Set_of_Molecules.molecule[].atom[]'
	atomList = []
	for m in range(self.totalmol):
		for a in range(self.size(udfpath, [m])):
			if self.get(udfpath+'.Atom_Name', [m,a]) == atomname:
				atomList.append((m, a, '{}:{}:{}'.format(atomname, m, a)))
	if len(atomList) == 0:
		print('no atom with Atom_Name: ' + atomname)
	return atomList























def read_chain2(rec):
	val.uobj.jump(rec)
	val.systemsize = val.uobj.get('Structure.Unit_Cell.Cell_Size.a')
	# bound_setup()
	# CU.setCell(tuple(val.uobj.get("Structure.Unit_Cell.Cell_Size")))
	ba = CognacBasicAnalysis(val.target, rec)
	mols = val.uobj.get('Structure.Position.mol[].atom[]')
	# ステップの数に対応した空リストを作成
	val.n_bonds = len(mols[0])
	r2_ij = [[] for i in range(val.n_bonds)]
	xp = [[] for i in range(val.n_bonds)]
	for chain in mols:
		for step in range(1, val.n_bonds):
			for start in range(val.n_bonds - step):
				end1 = tuple(chain[start])
				end2 = tuple(chain[start + step])
				# e2e_vec = CU.distanceWithBoundary(end1, end2)
				e2e_vec = np.array(end2) - np.array(end1)
				e2e_dist = np.linalg.norm(np.array(e2e_vec))
				
				r2 = e2e_dist**2
				r2_ij[step].append(r2)
				if step == 1:
					val.bond_list.append(e2e_dist)
				if step == val.n_bonds -1:
					val.Rx_list.append(e2e_vec[0])
					val.Ry_list.append(e2e_vec[1])
					val.Rz_list.append(e2e_vec[2])
					#
					val.R_list.append(e2e_dist)
		#


	# cn
	cn = []
	val.l_bond = np.average(np.array(val.bond_list))
	for i in range(1, len(r2_ij)):
		cn.append([i, np.average(np.array(r2_ij[i]))/(i*val.l_bond**2)])
	val.cn_list.append(cn)
	# angle
	anglename = val.uobj.get("Molecular_Attributes.Angle_Potential[].Name")
	tmp = np.array(ba.angle(anglename[0]))
	val.angle_list.extend(list(tmp[~np.isnan(tmp)]))


	return













##########################
# マルチリストのグラフの作成
def make_multi(cond_list):
	var.base_name = cond_list[0]
	var.data_list = cond_list[1]
	var.leg = cond_list[2]
	var.target_dir = os.path.join(var.target_name, var.base_name)
	var.f_dat = var.base_name + "_hist.dat"
	var.f_plt = var.base_name + ".plt"
	var.f_png = var.base_name + ".png"

	# データを書き出し 
	write_multi_data()
	# グラフを作成
	make_multi_graph()
	return

# データを書き出し 
def write_multi_data():
	os.makedirs(var.target_dir, exist_ok=True)
	with open(os.path.join(var.target_dir, var.f_dat), 'w') as f:
		if var.base_name == 'Xp':
			f.write("# time\tCp\tCp_x\tCp_y\tCp_z\n\n")
			for line in var.data_list:
				for data in line:
					f.write(str(data) + '\t')
				f.write('\n')
		else:
			for i, data in enumerate(var.data_list):
				f.write("\n\n# " + str(i) +":\n\n")
				for line in data:
					f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
	return

# グラフを作成
def make_multi_graph():
	make_multi_script()
	cwd = os.getcwd()
	os.chdir(var.target_dir)
	if platform.system() == "Windows":
		subprocess.call(var.f_plt, shell=True)
	elif platform.system() == "Linux":
		subprocess.call('gnuplot ' + var.f_plt, shell=True)
	os.chdir(cwd)
	return

# 必要なスクリプトを作成
def make_multi_script():
	with open(os.path.join(var.target_dir, var.f_plt), 'w') as f:
		script = multi_script_content()
		f.write(script)
	return

# スクリプトの中身
def multi_script_content():
	repeat = len(var.data_list )
	#
	script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
	script += '# \ndata = "' + var.f_dat + '" \nset output "' + var.f_png + ' "\n'
	script += '#\nset size square\nset xrange [0:]\nset yrange [0:]\n'
	script += '#\nset xlabel "' + var.leg[0] + '"\nset ylabel "' + var.leg[1] + '"\n\n'
	#
	if var.base_name == "Xp":
		script += 'plot data u 1:2 w l ti "Xp-ave.", \\\n'
		script += 'data u 1:3 w l ti "Xp-x", \\\n'
		script += 'data u 1:4 w l ti "Xp-y", \\\n'
		script += 'data u 1:5 w l ti "Xp-z"\n\nreset'

	elif var.base_name == 'Corr_stress' or var.base_name == 'Corr_stress_mod':
		script += 'set logscale xy \n\nset format x "10^{%L}" \nset format y "10^{%L}"\n\n'
		script += 'plot '
		script += 'data w l ti "Stress" \\\n'
	elif var.base_name == 'Corr_stress_semi':
		script += 'set logscale y \n\n#set format x "10^{%L}" \nset format y "10^{%L}"\n\n'
		script += 'a = 1\ntau =1000\n\ns = 100\ne = 1000\n\n'
		script += 'f(x) = a*exp(-1*x/tau) \n'
		script += 'fit [s:e] f(x) data usi 1:2 via a,tau\n\n'
		script += 'set label 1 sprintf("Fitted \\nA = %.1e \\n{/Symbol t} = %.1e \\nFitting Region: %d to %d", a, tau, s, e) at graph 0.35, 0.75\n\n'
		script += 'plot '
		script += 'data w l ti "Stress", \\\n'
		script += '[s:e] f(x) noti'
	elif var.base_name == 'Corr_stress_all':
		script += 'set logscale xy \n\nset format x "10^{%L}" \nset format y "10^{%L}"\n\n'
		script += 'plot data u 1:2 w l ti "G_t", \\\n'
		script += 'data u 1:3 w l ti "xy", \\\n'
		script += 'data u 1:4 w l ti "yz", \\\n'
		script += 'data u 1:5 w l ti "zx", \\\n'
		script += 'data u 1:6 w l ti "xx-yy", \\\n'
		script += 'data u 1:7 w l ti "yy-zz"'
	elif var.base_name == 'Sq':
		script += 'plot data u 1:2 w l ti "S(q)"'
	elif var.base_name == 'Guinier':
		script += 'set logscale y \n\nset format y "10^{%L}"\n\n'
		script += 'plot data u 1:($2**2) w l ti "G_t", \\\n'
	else:
		script += 'plot '
		for i in range(repeat):
			script += 'data ind ' + str(i) + ' w l lc ' + str(i) + 'noti, \\\n'

	return script



































def msd():
	samples = 10
	records = 7

	target = np.array([range(i, i+records) for i in np.arange(samples)])

	base = np.zeros(records)

	tlist = [[[-1. if i == number else 1.0 if i == number + step else x for i, x in enumerate(base)] if number < records - step else base for number in range(records - 1)] for step in range(1, records)]
	# tlist = []
	# for step in range(1, records):
	# 	tmp2 = []
	# 	for number in range(records-1):
	# 		tmp = []
	# 		for i, elm in enumerate(base):
	# 			if i == number:
	# 				tmp.append(-1.)
	# 			elif i == number+step:
	# 				tmp.append(1.)
	# 			else:
	# 				tmp.append(elm)
	# 		if number < records-step:
	# 			tmp2.append(tmp)
	# 		else:
	# 			tmp2.append(base)
	# 	tlist.append(tmp2)


	modar = np.array(tlist)


	norm_ar = np.array([1./x for x in reversed(range(1, records))])
	print(norm_ar)

	abs_d = np.abs(np.matmul(modar, np.transpose(target)))
	sum_data = np.sum(np.average(abs_d, axis = 2), axis = 1)
	ave_abs = np.multiply(sum_data, norm_ar)
	print(ave_abs)

	sqred_d = np.square(np.matmul(modar, np.transpose(target)))
	sum_sq_data = np.sum(np.average(sqred_d, axis = 2), axis = 1)
	ave_sq = np.multiply(sum_sq_data, norm_ar)
	print(ave_sq)












####################################################################################
# Green Kubo での計算を処理
def calc_gk():
	corr, corr_all = calc_corr()
	cond_list = ["Corr_stress", corr, ['Time', 'sigma']]
	make_multi(cond_list)
	cond_list = ["Corr_stress_all", corr_all, ['Time', 'sigma', 'ave']]
	make_multi(cond_list)

	# irheo(corr)
	return

def calc_corr():
	val.uobj.jump(val.uobj.totalRecord() - 1)
	#
	vol = val.uobj.get('Statistics_Data.Volume.Total_Average')
	corr_all = val.uobj.get('Correlation_Functions.Stress.Correlation[]')
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
		corr.append([time, g, ave])
	# g_mod = signal.savgol_filter(g, 3, 2)
	# corr_mod = np.stack([time, g_mod], 1)
	return corr, corr_all


##################################
# 
def irheo(data_list):
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












def test():

	samples = 10
	records = 7

	target = np.array([range(i, i+records) for i in np.arange(samples)])

	base = np.zeros(records)

	tlist = [[[-1. if i == number else 1.0 if i == number + step else x for i, x in enumerate(base)] if number < records - step else base for number in range(records - 1)] for step in range(1, records)]
	# tlist = []
	# for step in range(1, records):
	# 	tmp2 = []
	# 	for number in range(records-1):
	# 		tmp = []
	# 		for i, elm in enumerate(base):
	# 			if i == number:
	# 				tmp.append(-1.)
	# 			elif i == number+step:
	# 				tmp.append(1.)
	# 			else:
	# 				tmp.append(elm)
	# 		if number < records-step:
	# 			tmp2.append(tmp)
	# 		else:
	# 			tmp2.append(base)
	# 	tlist.append(tmp2)


	modar = np.array(tlist)


	norm_ar = np.array([1./x for x in reversed(range(1, records))])
	print(norm_ar)

	abs_d = np.abs(np.matmul(modar, np.transpose(target)))
	sum_data = np.sum(np.average(abs_d, axis = 2), axis = 1)
	ave_abs = np.multiply(sum_data, norm_ar)
	print(ave_abs)

	sqred_d = np.square(np.matmul(modar, np.transpose(target)))
	sum_sq_data = np.sum(np.average(sqred_d, axis = 2), axis = 1)
	ave_sq = np.multiply(sum_sq_data, norm_ar)
	print(ave_sq)

def evaluate_all():
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




def array_test():
	t=5
	tlist = []
	base = np.zeros(t)
	for i, elm in enumerate(base):
		if i == 0:
			tlist.append(-1.)
		elif i == 1:
			tlist.append(1.)
		else:
			tlist.append(elm)

	tlist2 = [-1. if i ==0 else x for i, x in enumerate(base)]

	print(tlist)
	print(tlist2)
	# test = np.array()

	return