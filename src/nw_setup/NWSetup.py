#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
# from ast import arg
import re
import numpy as np
import copy
import random
import platform
import subprocess
import sys
import os
import pickle
import multiprocessing as mp
from UDFManager import UDFManager
################################################################################
#
def select_set(nw_cond, sim_cond2, rnd_cond, target_dir):
	nw_model = nw_cond[0]
	# ネットワークを設定
	if nw_model == "Regular":
		calcd_data_dic = regnw_setup(nw_cond, sim_cond2)
	elif nw_model == "Random":
		calcd_data_dic = rndnw_setup(nw_cond, sim_cond2, rnd_cond, target_dir)

	return calcd_data_dic

#######################################
def regnw_setup(nw_cond, sim_cond2):
	calcd_data_dic = calc_all(nw_cond, sim_cond2)
	return calcd_data_dic

def rndnw_setup(nw_cond, sim_cond2, rnd_cond, target_dir):
	multi = sim_cond2[0]

	restart = rnd_cond[0] 
	cond_top = rnd_cond[1]
	hist_bin = rnd_cond[2]
	base_top_list = make_8chain_dic(nw_cond)

	if restart == '':
		# トポロジーの異なるネットワークを探索して、任意の多重度のネットワークポリマーの代数的連結性の分布関数を策定
		candidate_list, random_dir = top_search(nw_cond, cond_top, base_top_list)
	else:
		# 設定したリスタートファイルを読み込んで、リストを作成
		candidate_list, random_dir = top_select(restart)

	top_dic_list = nw_search(candidate_list, random_dir, hist_bin, multi, target_dir)

	###########################################
	# ターゲットとなるネットワーク全体の辞書を設定。
	calcd_data_dic = make_data_dic(nw_cond, base_top_list, top_dic_list)
	return calcd_data_dic


################################################
## REGULAR NW SETUP
################################################
def calc_all(nw_cond, sim_cond2):
	# 架橋点 JP を設定
	jp_xyz, strand_se_xyz = calc_jp_strands(nw_cond)
	#
	calcd_data_dic = set_atom(sim_cond2, nw_cond, jp_xyz, strand_se_xyz)
	return calcd_data_dic
##########################################
# JPおよびサブチェインの始点と終点のXYZを設定
def calc_jp_strands(nw_cond):
	strand = nw_cond[1]
	# jp_xyz は、JPの座標のリスト
	# strand_se_xyz は、サブチェインの出発点と終点のリスト
	if strand == "3_Chain_S":
		# JPを設定
		jp_xyz = [
		[
		[0, 0, 0],
		[0, 0.25, 0.25],
		[0.25, 0.25, 0.5],
		[0.25, 0, 0.75],
		[0.5, 0.5, 0.5],
		[0.5, 0.75, 0.75],
		[0.75, 0.5, 0.25],
		[0.75, 0.75, 0]
		]
		]
		# サブチェインの出発点と終点を設定
		strand_se_xyz = [
		[
		[[0, 0, 0], [0, 0.25, 0.25]],
		[[0, 0.25, 0.25], [0.25, 0.25, 0.5]],
		[[0.25, 0.25, 0.5], [0.25, 0, 0.75]],
		[[0.25, 0.25, 0.5], [0.5, 0.5, 0.5]],
		[[0.5, 0.5, 0.5], [0.5, 0.75, 0.75]],
		[[0.5, 0.5, 0.5], [0.75, 0.5, 0.25]],
		[[0.75, 0.5, 0.25], [0.75, 0.75, 0]],
		[[0.75, 0.5, 0.25], [1, 0.25, 0.25]],
		[[0.25, 0, 0.75], [0, 0, 1]],
		[[0.5, 0.75, 0.75], [0.25, 1, 0.75]],
		[[0.5, 0.75, 0.75], [0.75, 0.75, 1]],
		[[0.75, 0.75, 0], [1, 1, 0]]
		]
		]

	elif strand == "3_Chain_D":
		# JPを設定
		jp_xyz = [
		[
		[0, 0, 0],
		[0, 0.25, 0.25],
		[0.25, 0.25, 0.5],
		[0.25, 0, 0.75],
		[0.5, 0.5, 0.5],
		[0.5, 0.75, 0.75],
		[0.75, 0.5, 0.25],
		[0.75, 0.75, 0]
		],
		[	#ここから二つ目
		[0, 0.5, 0.75],
		[0, 0.75, 0.5],
		[0.25, 0.75, 0.25],
		[0.25, 0.5, 0],
		[0.5, 0.25, 0],
		[0.5, 0, 0.25],
		[0.75, 0, 0.5],
		[0.75, 0.25, 0.75]
		]
		]
		# サブチェインの出発点と終点を設定
		strand_se_xyz = [
		[
		[[0, 0, 0], [0, 0.25, 0.25]],
		[[0, 0.25, 0.25], [0.25, 0.25, 0.5]],
		[[0.25, 0.25, 0.5], [0.25, 0, 0.75]],
		[[0.25, 0.25, 0.5], [0.5, 0.5, 0.5]],
		[[0.5, 0.5, 0.5], [0.5, 0.75, 0.75]],
		[[0.5, 0.5, 0.5], [0.75, 0.5, 0.25]],
		[[0.75, 0.5, 0.25], [0.75, 0.75, 0]],
		[[0.75, 0.5, 0.25], [1, 0.25, 0.25]],
		[[0.25, 0, 0.75], [0, 0, 1]],
		[[0.5, 0.75, 0.75], [0.25, 1, 0.75]],
		[[0.5, 0.75, 0.75], [0.75, 0.75, 1]],
		[[0.75, 0.75, 0], [1, 1, 0]]
		],
		[
		[[0, 0.5, 0.75], [0, 0.75, 0.5]],
		[[0, 0.75, 0.5], [0.25, 0.75, 0.25]],
		[[0.25, 0.75, 0.25], [0.25, 0.5, 0]],
		[[0.25, 0.5, 0], [0.5, 0.25, 0]],
		[[0.5, 0.25, 0], [0.5, 0, 0.25]],
		[[0.5, 0, 0.25], [0.75, 0, 0.5]],
		[[0.75, 0, 0.5], [0.75, 0.25, 0.75]],
		[[0, 0.5, 0.75], [0.25, 0.5, 1]],
		[[0.25, 0.75, 0.25], [0.5, 1, 0.25]],
		[[0.75, 0.25, 0.75], [0.5, 0.25, 1]],
		[[0.75, 0.25, 0.75], [1, 0.5, 0.75]],
		[[0.75, 1, 0.5], [1, 0.75, 0.5]]
		]
		]

	elif strand == "4_Chain":
		# JPを設定
		jp_xyz = [
		[
		[0.0, 0.0, 0.0],
		[0.0, 0.5, 0.5],
		[0.5, 0.0, 0.5],
		[0.5, 0.5, 0.0],
		[0.25, 0.25, 0.25],
		[0.25, 0.75, 0.75],
		[0.75, 0.25, 0.75],
		[0.75, 0.75, 0.25]
		]
		]
		# サブチェインの出発点と終点を設定
		strand_se_xyz = [
		[
		[[0.25, 0.25, 0.25], [0.0, 0.0, 0.0]],	# No.1
		[[0.25, 0.25, 0.25], [0.0, 0.5, 0.5]],
		[[0.25, 0.25, 0.25], [0.5, 0.0, 0.5]],
		[[0.25, 0.25, 0.25], [0.5, 0.5, 0.0]],
		[[0.25, 0.75, 0.75], [0.0, 0.5, 0.5]],	# No.2
		[[0.25, 0.75, 0.75], [0.0, 1.0, 1.0]],
		[[0.25, 0.75, 0.75], [0.5, 0.5, 1.0]],
		[[0.25, 0.75, 0.75], [0.5, 1.0, 0.5]],
		[[0.75, 0.25, 0.75], [0.5, 0.0, 0.5]],	# No.3
		[[0.75, 0.25, 0.75], [0.5, 0.5, 1.0]],
		[[0.75, 0.25, 0.75], [1.0, 0.0, 1.0]],
		[[0.75, 0.25, 0.75], [1.0, 0.5, 0.5]],
		[[0.75, 0.75, 0.25], [0.5, 0.5, 0.0]],	# No.4
		[[0.75, 0.75, 0.25], [0.5, 1.0, 0.5]],
		[[0.75, 0.75, 0.25], [1.0, 0.5, 0.5]],
		[[0.75, 0.75, 0.25], [1.0, 1.0, 0.0]]
		]
		]

	elif strand == "6_Chain":
		# JPを設定
		jp_xyz = [
		[
		[0.,0.,0.]
		]
		]
		# サブチェインの出発点と終点を設定
		strand_se_xyz = [
		[
		[[0., 0., 0.], [1, 0, 0]],
		[[0., 0., 0.], [0, 1, 0]],
		[[0., 0., 0.], [0, 0, 1]]
		]
		]

	elif strand == "8_Chain":
		# JPを設定
		jp_xyz = [
		[
		[0.,0.,0.],
		[0.5,0.5,0.5]
		]
		]
		# サブチェインの出発点と終点を設定
		strand_se_xyz = [
		[
		[[0.5, 0.5, 0.5], [0, 0, 0]],
		[[0.5, 0.5, 0.5], [1, 0, 0]],
		[[0.5, 0.5, 0.5], [0, 1, 0]],
		[[0.5, 0.5, 0.5], [1, 1, 0]],
		[[0.5, 0.5, 0.5], [0, 0, 1]],
		[[0.5, 0.5, 0.5], [1, 0, 1]],
		[[0.5, 0.5, 0.5], [0, 1, 1]],
		[[0.5, 0.5, 0.5], [1, 1, 1]]
		]
		]

	return jp_xyz, strand_se_xyz

#########################################################
def set_atom(sim_cond2, nw_cond, jp_xyz, strand_se_xyz):
	multi = sim_cond2[0]

	n_segments = nw_cond[3]
	n_cell = nw_cond[4]
	n_sc = nw_cond[5]

	calcd_data_dic={}

	for mul in (range(multi)):
		# atom_all = []
		# pos_all = {}
		for mol, jp in enumerate(jp_xyz):
			atom_all = []
			pos_all = {}
			# システム全体にわたるジャンクションポイントのxyzとIDの辞書を作成
			jp_id_dic, jp_xyz_dic, atom_jp = set_jp_id_reg(jp, n_cell)
			atom_all.extend(atom_jp)
			pos_all.update(jp_xyz_dic)
			# サブチェイン中の各アトムのxyzリストとボンドリストを作成
			strand_xyz, bond_all, atom_sc, angle_all = set_strands(jp_id_dic, strand_se_xyz[mol], n_cell, n_segments, n_sc)
			#
			atom_all.extend(atom_sc)
			pos_all.update(strand_xyz)
			#
			calcd_data_dic[mul+mol] = {"atom_all":atom_all, "bond_all":bond_all, "pos_all":pos_all, "angle_all":angle_all}
	return calcd_data_dic

###################################################
# システム全体にわたるJPのxyzとIDの辞書を作成
def set_jp_id_reg(jp_xyz, n_cell):
	jp_id_dic = {}
	jp_xyz_dic = {}
	atom_jp = []
	jp_id = 0
	for z in range(n_cell):
		for y in range(n_cell):
			for x in range(n_cell):
				base_xyz = np.array([x,y,z])
				for jp in jp_xyz:
					jp_id_dic[tuple(np.array(jp) + base_xyz)] = (jp_id)
					jp_xyz_dic[(jp_id)] = tuple(np.array(jp) + base_xyz)
					atom_jp.append([jp_id, 0, 0])
					jp_id += 1
	return jp_id_dic, jp_xyz_dic, atom_jp

#########################################################
# サブチェイン中の各アトムのxyzリストとボンドリストを作成
def set_strands(jp_id_dic, strand_se_xyz, n_cell, n_segments, n_sc):
	strand_xyz = {}
	bond_all = {}
	atom_sc = []
	angle_all = []
	sub_id = len(jp_id_dic)
	bond_id = 0
	for z in range(n_cell):
		for y in range(n_cell):
			for x in range(n_cell):
				b_xyz = (x,y,z)
				for se_xyz in strand_se_xyz:
					tmp_xyz, tmp_bond, new_sub_id, new_bond_id, tmp_atom_sc, tmp_angle = calc_single_strand_reg(jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, n_segments, n_cell, n_sc)
					strand_xyz.update(tmp_xyz)
					bond_all.update(tmp_bond)
					atom_sc.extend(tmp_atom_sc)
					angle_all.append(tmp_angle)
					sub_id = new_sub_id
					bond_id = new_bond_id
	return strand_xyz, bond_all, atom_sc, angle_all

###############################################################
# 一本のサブチェイン中の各アトムのxyzリストとボンドリストを作成
def calc_single_strand_reg(jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, n_segments, n_cell, n_sc):
	tmp_xyz = {}
	tmp_bond = {}
	tmp_angle = []
	tmp_atom_sc = []
	bas_xyz = np.array(b_xyz)
	# サブチェインの末端間のベクトルを設定
	start_xyz = np.array(se_xyz[0]) + bas_xyz
	end_xyz = np.array(se_xyz[1]) + bas_xyz
	vec = end_xyz - start_xyz
	# ストランドの鎖長分のループ処理
	unit_len = 1./(n_segments + 1)
	ortho_vec = find_ortho_vec(vec)
	mod_o_vec = np.linalg.norm(vec)*ortho_vec
	# 始点のアトムのIDを設定
	mod_xyz = list(start_xyz)[:]
	for dim in range(3):
		if mod_xyz[dim] == n_cell:
			mod_xyz[dim] = 0
	s_id = jp_id_dic[tuple(mod_xyz)]
	tmp_angle.append(s_id)
	# 終点のアトムのIDを周期境界条件で変更
	mod_xyz = list(end_xyz)[:]
	for dim in range(3):
		if mod_xyz[dim] == n_cell:
			mod_xyz[dim] = 0
	E_id = jp_id_dic[tuple(mod_xyz)]
	# サブチェインの鎖長分のループ処理
	for seg in range(n_segments):
		pos = tuple(start_xyz + vec*(seg + 1)/(n_segments + 1.))
		tmp_xyz[sub_id] = pos

		if seg == 0 or seg == n_segments - 1:
			tmp_atom_sc.append([sub_id, 1, 1])
		else:
			tmp_atom_sc.append([sub_id, 2, 2])
		e_id = sub_id
		#
		if seg == 0:
			bond = 0
		else:
			bond = 1
		tmp_bond[bond_id] = tuple([bond, [s_id, e_id]])
		bond_id += 1
		tmp_angle.append(e_id)
		s_id = e_id
		sub_id += 1
		
		if n_sc != 0:
			sc_s_id = s_id
			for i in range(n_sc):
				tmp_xyz[sub_id] = tuple(np.array(pos)  +  (i + 1)*mod_o_vec*unit_len)
				tmp_atom_sc.append([sub_id, 2, 1])
				sc_e_id = sub_id
				#
				bond = 2
				tmp_bond[bond_id] = tuple([bond, [sc_s_id, sc_e_id]])
				sc_s_id = sc_e_id
				sub_id += 1
				bond_id += 1
		
	e_id = E_id
	bond = 0
	tmp_bond[bond_id] = tuple([bond, [s_id, e_id]])
	tmp_angle.append(e_id)
	bond_id += 1
	return tmp_xyz, tmp_bond, sub_id, bond_id, tmp_atom_sc, tmp_angle

######
def find_ortho_vec(vec2):
	vec = np.array(vec2).reshape(-1,1)
	# 線形独立である新たな三次元ベクトルを見つける。
	rank = 0
	while rank != 2:
		a = np.array(np.random.rand(3)).reshape(-1,1)
		target = np.hstack((vec, a))
		rank = np.linalg.matrix_rank(target)
	# QR分解により
	q, r = np.linalg.qr( target )
	ortho_vec = q[:,1]
	return ortho_vec


################################################################################
## RANDOM NW SETUP
################################################################################

# 基本構造として、8Chainモデルを設定
def make_8chain_dic(nw_cond):
	nw_model = nw_cond[0]
	strand = nw_cond[1]
	n_strand = nw_cond[2]
	n_segments = nw_cond[3]
	n_cell = nw_cond[4]
	n_sc = nw_cond[5]
	# ユニットセルでの、jp およびサブチェインの始点と終点のXYZを設定
	jp_xyz = [
			[0.,0.,0.], 
			[0.5, 0.5, 0.5]
			]
	strand_se_xyz = [
				[[0.5, 0.5, 0.5], [0, 0, 0]],
				[[0.5, 0.5, 0.5], [1, 0, 0]],
				[[0.5, 0.5, 0.5], [0, 1, 0]],
				[[0.5, 0.5, 0.5], [1, 1, 0]],
				[[0.5, 0.5, 0.5], [0, 0, 1]],
				[[0.5, 0.5, 0.5], [1, 0, 1]],
				[[0.5, 0.5, 0.5], [0, 1, 1]],
				[[0.5, 0.5, 0.5], [1, 1, 1]]
				]
	start_jp = [0.5, 0.5, 0.5]
	# システム全体にわたるピボットのxyzとIDの辞書を作成
	jp_id_dic, jp_xyz_dic, atom_jp, n_jp = set_jp_id_rnd(jp_xyz, n_cell)
	# ストランドの結合状態を記述
	init_8ch_dic, vector_dic = set_strands_8(jp_id_dic, strand_se_xyz, n_cell, start_jp)
	#
	base_top_list = [init_8ch_dic, jp_xyz_dic, atom_jp, n_jp, vector_dic]

	return base_top_list 

##########################################
# システム全体にわたるjpのxyzとIDの辞書を作成
def set_jp_id_rnd(jp_xyz, n_cell):
	jp_id = 0
	jp_id_dic = {}
	jp_xyz_dic = {}
	atom_jp = []
	for z in range(n_cell):
		for y in range(n_cell):
			for x in range(n_cell):
				base_xyz = np.array([x,y,z])
				for jp in jp_xyz:
					jp_id_dic[tuple(np.array(jp) + base_xyz)] = (jp_id)
					jp_xyz_dic[jp_id] = tuple(np.array(jp) + base_xyz)
					atom_jp.append([jp_id, 0, 0])
					jp_id += 1
	n_jp = jp_id
	return jp_id_dic, jp_xyz_dic, atom_jp, n_jp

#####################################################
# ストランドの結合状態を記述
def set_strands_8(jp_id_dic, strand_se_xyz, n_cell, start_jp):
	init_8ch_dic = {}
	vector_dic = {}
	str_id = 0
	for z in range(n_cell):
		for y in range(n_cell):
			for x in range(n_cell):
				base_xyz = np.array([x,y,z])
				for xyz in strand_se_xyz:
					# サブチェインの末端間のベクトルを設定
					start_xyz = np.array(start_jp) + base_xyz
					end_xyz = np.array(xyz[1]) + base_xyz
					vector = end_xyz - start_xyz
					# 始点のアトムのIDを設定
					start_id = jp_id_dic[tuple(start_xyz)]
					# 終点のアトムのIDを周期境界条件で変更
					mod_end_xyz = list(end_xyz)[:]
					for i in range(3):
						if mod_end_xyz[i] == n_cell:
							mod_end_xyz[i] = 0
					end_id = jp_id_dic[tuple(mod_end_xyz)]
					#
					init_8ch_dic[str_id] = tuple([start_id, end_id])
					vector_dic[str_id] = vector
					str_id+=1
	return init_8ch_dic, vector_dic


################################################################################
## トポロジーの異なるネットワークを探索して、代数的連結性の分布関数を策定
################################################################################
#########################################################
# トポロジーの異なるネットワークを探索して、代数的連結性の分布関数を策定し、ネットワークトポロジーの配列辞書を決める。
def top_search(nw_cond, cond_top, base_top_list):
	n_strand = nw_cond[2]
	n_cell = nw_cond[4]
	#
	init_dic = base_top_list[0]
	n_jp = base_top_list[3]
	#
	n_sampling = cond_top[1]
	n_try = cond_top[2]
	repeat = cond_top[3]
	#
	random_dir = str(n_strand) +"_chains_" + str(n_cell) + "_cells_"
	random_dir += str(n_sampling) + "_sampling_" + str(n_try) + "_trials_" + str(repeat) + "_times"
	if os.path.isdir(random_dir):
		print("#####\nRandom Calculation target dir exists!!\n")
		while True:
			choice = input("overwrite? [y/N]:").lower()
			if choice in ['y', 'ye', 'yes']:
				os.makedirs(random_dir, exist_ok = True)
				break
			else:
				sys.exit('Bye now !')
	else:
		os.makedirs(random_dir, exist_ok = True)
	#
	candidate_list = strand_exchange(init_dic, n_jp, cond_top, n_strand, n_cell)
	#
	with open(os.path.join(random_dir, 'init.pickle'), mode = 'wb') as f:
		pickle.dump(candidate_list, f)
	#
	print("##################################################")
	print("Trial, Sampling, Repeat = ", n_try, ", ", n_sampling, ", ", repeat)
	print("Total Sampling = ", n_try*n_sampling*repeat)
	print("Initial Candidates = ", len(candidate_list))
	print("##################################################")

	return candidate_list, random_dir

#####################################################
# 任意のストランドを選択し、ストランドの繋ぎ変えを行う
def strand_exchange(init_dic, n_jp, cond_top, n_strand, n_cell):
	pre_sampling = cond_top[0] 
	n_sampling = cond_top[1]
	n_try = cond_top[2]
	repeat = cond_top[3]
	f_pool = cond_top[4]
	
	if mp.cpu_count() < f_pool:
		f_pool = int(mp.cpu_count()/2) - 1
	p = mp.Pool(f_pool)

	pre_list = []
	args = []
	# "self.random_reduce()" により任意のストランドを除去したものをサンプリング
	print("Searching for Initial Random Structure of ", pre_sampling)
	print("Please Wait a little bit!")
	if f_pool > 1:
		for i in range(pre_sampling):
			args.append([random_reduce(init_dic, n_strand, n_cell, n_jp).copy(), i, "Pre", init_dic, n_jp, n_try])
		result = p.map(search, args)
		for data in result:
			pre_list.extend(data)
	else:
		for i in range(pre_sampling):
			args = [random_reduce(init_dic, n_strand, n_cell, n_jp).copy(), i, "Pre", init_dic, n_jp, n_try]
			tmp = search(args)
			pre_list.extend(tmp)
	print("##################################################")
	print("Pre_Search Result:", len(pre_list))
	print("##################################################")

	# pre_list を alg_const の昇順に並べ替えて、sortbyalgを作成
	tmp_array = np.array(pre_list)
	sortbyalg = tmp_array[np.argsort(tmp_array[:,1])]
	#
	candidate_list = []
	step = int(len(sortbyalg)/n_sampling)
	for count in range(repeat):
		args = []
		for i in range(n_sampling):
			# サンプリング数が十分なときに、sortbyalg を代数的連結性の低い方から分割してサンプリング
			if len(sortbyalg) > n_sampling*repeat:
				part = sortbyalg[step*i:step*(i+1)]
				target =list(random.choice(part))
			else:# サンプリング数が少ないときは、sortbyalg全体からランダムサンプリング
				target =list(random.choice(sortbyalg))
			if f_pool > 1:
				args.append([target[0].copy(), i, str(count+1)+"/"+str(repeat), init_dic, n_jp, n_try])
			else:
				args = [target[0].copy(), i, str(count+1)+"/"+str(repeat), init_dic, n_jp, n_try]
				candidate_list.extend(search(args))
		
		if f_pool > 1:
			result = p.map(search, args)
			for data in result:
				candidate_list.extend(data)
		print("##################################################")
		print(str(count)+"/"+str(repeat), "Search Result:", len(candidate_list))
		print("##################################################")

	return candidate_list

########################################################
# ストランドの繋ぎ変えを行う
def search(args):
	target_dic, x, id, init_dic, n_jp, n_try = args
	alg_const = calc_lap_mat(target_dic, n_jp)
	print(id, "Sampling ID =", x,  "Initial Algebratic Connectivity =", alg_const)
	#
	result = []
	count = 0
	failed = 0
	show = 5
	while count < n_try:
		# 現状のストランドのリストの中からランダムに一つ選択し、"selected_strand"とする。
		selected_strand = random.choice(list(target_dic.keys()))
		# 繋ぎ変え得るストランドのリストを現状のネットワークと比較して、交換可能なセットを見つける。
		tmp_dic, alg_const = find_pair(selected_strand, target_dic, init_dic, n_jp)
		if alg_const != 0:
			count += 1
			result.append([tmp_dic.copy(), alg_const])
			failed = 0
			if count != 0 and round(show*count/n_try) == show*count//n_try and round(show*count/n_try) == -(-show*count//n_try):
				print(id, "Sampling ID =", x, "count = ", count)
		else:
			failed +=1

		if failed >= n_try:
			print("##########################################")
			print(id, "Sampling ID =", x,  " FAILED!! with ", failed, "th trials.")
			print("##########################################")
			count = n_try
			failed = 0

	return result

#########################################################
# 任意のストランドを選択し、所望の分岐数にストランドを消去
def random_reduce(init_dic, n_strand, n_cell, n_jp):
	# alg_const_init = calc_lap_mat(init_dic, n_jp)
	# flag = 1
	# while flag == 1:
	tmp_dic = copy.deepcopy(init_dic)
	del_bond = []
	# 消去対象のストランドをリストアップ
	# ユニットセル中のストランドから必要な数だけ抽出
	tmp_del = random.sample(range(8), (8 - n_strand))
	# 全セルに渡って、消去するストランドのリストを作成
	for i in range(n_cell**3):
		del_bond.extend(list(8*i + np.array(tmp_del)))
	# ストランドを消去した辞書とその代数的連結性を得る。
	for target in del_bond:
		tmp_dic.pop(target)
		# alg_const = calc_lap_mat(tmp_dic, n_jp)

		# if alg_const_init <= alg_const:
			# break
	return tmp_dic

################################################################
# 交換可能なストランドのペアを見つける
def find_pair(selected_strand, target_dic, init_dic, n_jp):
	deleted_dic = dict(init_dic.items() - target_dic.items())
	# 選択したストランドの両端のjp（"jp_pairs"）を見つけ、
	# それと繋がり得る可能性のあるjpをリストアップ（"connected_list"）
	# 繋ぎ変え得るストランドを見つける（"possible_strand"）
	possible_strand = get_connected_jp(init_dic, selected_strand, n_jp)
	# 繋ぎ変え得るストランドのリストを現状のネットワークと比較して、交換可能なセットを見つける。
	found_dic = {}
	result_dic = {}
	for p_str in possible_strand:
		if p_str in target_dic:
			# 現行のストランドリストの中にいないものを選択
			count = 0
			tmp_dic = {}
			for tt in target_dic[selected_strand]:
				for uu in target_dic[p_str]:
					tmp_str = []
					if tt != uu:
						if (tt,uu) in deleted_dic.values():
							tmp_str = [k for k, v in init_dic.items() if v == (tt, uu)]
						elif (uu,tt) in deleted_dic.values():
							tmp_str = [k for k, v in init_dic.items() if v == (uu, tt)]
						if tmp_str != []:
							tmp_dic[tmp_str[0]] = tuple(init_dic[tmp_str[0]])
							count +=1
			if count == 2 and tmp_dic != {}:
				found_dic.update(tmp_dic)
				target_dic.pop(p_str)
				target_dic.pop(selected_strand)
				target_dic.update(tmp_dic)
				result_dic = target_dic
				alg_const = calc_lap_mat(result_dic, n_jp)
				return result_dic, alg_const
			else:
				alg_const = 0
		else:
				alg_const = 0
	return result_dic, alg_const

#################################################################
# 任意のjpと連結したjpを見つけ、繋がり得る可能性のあるjpをリストアップ
def get_connected_jp(init_dic, selected_strand, n_jp):
	# 選択したストランドの両端のjp（"jp_pairs"）を見つけ、
	jp_pairs = tuple(init_dic[selected_strand])
	# ラプラシアン行列を利用して、それと繋がり得る可能性のあるjpをリストアップ（"connected_list"）
	lap_mat = make_lap_mat(init_dic, n_jp)
	connected_list = []
	for t_jp in jp_pairs:
		t_list = list(lap_mat[t_jp])
		connected_list.append(find_index(t_list))
	# 繋ぎ変え得るストランドを見つける。
	possible_strand = []
	for i in connected_list[0]:
		for j in connected_list[1]:
			link = get_link(init_dic, i, j)
			if link:
				possible_strand.append(link[0])
	return possible_strand

############################################################################
# ラプラシアン行列を利用して、任意のjpとストランドによりつながる８個のjpを見つける。
def find_index(t_list):
	con_list = [i for i, x in enumerate(t_list) if x == -1.0]
	return con_list

#######################################################
# 任意の二つのjpがストランドによりつながっている可能性を調査
# 繋がり得るものであれば、そのストランドのID番号を返す。
def get_link(init_dic, jp0, jp1):
	val_tpl = tuple([jp0, jp1])
	link = [k for k, v in init_dic.items() if v == val_tpl]
	if link:
		return link
	elif link == []:
		val_tpl = tuple([jp1, jp0])
		link = [k for k, v in init_dic.items() if v == val_tpl]
		if link:
			return link
	else:
		return []

####################################################
# 任意のネットワークトポロジーから、ラプラシアン行列を作成
def calc_lap_mat(topl_dic, n_jp):
	lap_mat = make_lap_mat(topl_dic, n_jp)
	# 固有値を計算
	la, v = np.linalg.eig(lap_mat)
	# 固有値の二番目の値を見つける。
	alg_const = round(sorted(np.real(la))[1], 3)
	return alg_const

####################################################
# 任意のネットワークトポロジーから、ラプラシアン行列を作成
def make_lap_mat(topl_dic, n_jp):
	lap_mat = np.zeros((n_jp, n_jp))
	for i in topl_dic:
		lap_mat[topl_dic[i][0], topl_dic[i][1]] = -1
		lap_mat[topl_dic[i][1], topl_dic[i][0]] = -1
		lap_mat[topl_dic[i][0], topl_dic[i][0]] += 1
		lap_mat[topl_dic[i][1], topl_dic[i][1]] += 1
	return lap_mat

#########################################################
# 過去の探索データを使って、代数的連結性の分布関数を選択
def top_select(restart):
	with open(os.path.join(restart, 'init.pickle'), mode = 'rb') as f:
		candidate_list = pickle.load(f)
	print("##################################################")
	print("Reloaded Candidates = ", len(candidate_list))
	print("##################################################")

	return candidate_list, restart


#####################################################################
# ヒストグラム中の最大頻度を与えるネットワークトポロジーの配列辞書を決める。
def nw_search(candidate_list, random_dir, hist_bin, multi, target_dir):
	tmp_list = []
	val_list = []
	data_list = []

	# ヒストグラムを作成
	histdata = list(np.array(candidate_list)[:, 1])
	cond = ["", histdata, hist_bin, "True", ['Arg. Con.', 'Freq.']]
	x, val = make_hist_all(cond, random_dir)

	# 最頻値のレンジを決める
	val_range = find_range(x, val, multi)
	# 上記のレンジに入る配列をピックアップ
	for i in candidate_list:
		if i[1] >= val_range[0] and i[1] <= val_range[1]:
			tmp_list.append(i)
	random.shuffle(tmp_list)
	# 
	count = 0
	for i, selected_list in enumerate(tmp_list):
		if selected_list[0] not in val_list:
			val_list.append(selected_list[1])
			data_list.append(selected_list[0])
			count += 1
		if len(val_list) == multi:
			u = UDFManager(os.path.join(target_dir, 'target_condition.udf'))
			u.put(random_dir, 'TargetCond.Model.RandomData')
			u.put(val_list, 'TargetCond.Model.SelectedValue[]')
			u.write()
			with open(os.path.join(random_dir, 'selected_val.dat'), 'w') as f:
				f.write("Selected arg. con.\n\n")
				for i in val_list:
					print("Selected arg. con.", round(i, 4))
					f.write("arg. con. = " + str(round(i, 4)) + '\n')
			return data_list
	#
	print("No effective list was found for multi numbers of", multi, "!  Try again!!")
	sys.exit()


#######################
# 最大頻度の範囲を決める。
def find_range(x, val, multi):
	index = np.where(val == max(val))[0][0]
	f_range = 0
	value = 0
	if index < len(val) -1:
		while value < multi:
			value = 0
			f_range += 1
			for i in range(2*f_range + 1):
				if val[index - f_range + i] != 0:
					value += 1
		val_range = [ x[index - f_range], x[index + f_range] ]
		print("##################################################")
		print("Most frequent range = ", round(val_range[0], 5), " to ", round(val_range[1], 5))
		key = input("OK? input y or n: ")
		if key == 'y' or key == 'Y':
			pass
		else:
			print("##################################################")
			print("Input new range:")
			low = float(input("low=: "))
			high = float(input("high=: "))
			val_range = list([low, high])
	else:
		print("##################################################")
		print("Input new range:")
		low = float(input("low=: "))
		high = float(input("high=: "))
		val_range = list([low, high])
	return val_range

###########################
# ヒストグラムのグラフの作成
def make_hist_all(cond_list, random_dir):
	# base = cond_list[0]
	list = cond_list[1]
	bin = cond_list[2]
	norm = cond_list[3]
	leg = cond_list[4]

	f_dat = "nw_hist.dat"
	f_plt = "make_hist.plt"
	f_png = "histgram.png"
	# ヒストグラムのデータ作成
	bin_width, hist_data, val, x = make_hist_data(list, norm, bin)
	# ヒストグラムのデータを書き出し 
	write_data(hist_data, random_dir, f_dat)
	# グラフを作成
	make_graph(bin_width, random_dir, f_plt, leg, f_png, f_dat)
	return x, val

# ヒストグラムのデータ作成
def make_hist_data(list, norm, bin):
	# ヒストグラムを作成
	weight = np.ones(len(list))/float(len(list))
	if norm:
		val, x = np.histogram(list, bins=bin, weights=weight)
	else:
		val, x = np.histogram(list, bins=bin)
	# グラフ用にデータを変更
	bin_width = (x[1]-x[0])
	mod_x = (x + bin_width/2)[:-1]
	hist_data = np.stack([mod_x, val], axis = 1)
	return bin_width, hist_data, val, x

# ヒストグラムのデータを書き出し 
def write_data(hist_data, random_dir, f_dat):
	os.makedirs(random_dir, exist_ok=True)
	with open(os.path.join(random_dir, f_dat), 'w') as f:
		f.write("# Histgram data:\n\n")
		for line in hist_data:
			f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
	return

# グラフを作成
def make_graph(bin_width, random_dir, f_plt, leg, f_png, f_dat):
	make_script(bin_width, random_dir, leg, f_png, f_dat, f_plt)
	cwd = os.getcwd()
	os.chdir(random_dir)
	if platform.system() == "Windows":
		subprocess.call(f_plt, shell=True)
	elif platform.system() == "Linux":
		subprocess.call('gnuplot ' + f_plt, shell=True)
	os.chdir(cwd)
	return
	
# 必要なスクリプトを作成
def make_script(bin_width, random_dir, leg, f_png, f_dat, f_plt):
	script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
	script += '# \ndata = "' + f_dat + '" \nset output "' + f_png + ' "\n'
	script += '#\nset size square\n# set xrange [0:1.0]\n#set yrange [0:100]\n'
	script += '#\nset xlabel "' + leg[0] + '"\nset ylabel "' + leg[1] + '"\n\n'
	script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
	script += '#\nplot data w boxes noti'
	print("OK")
	with open(os.path.join(random_dir, f_plt), 'w') as f:
		f.write(script)
	return
	
################################################################################
# ターゲットとなるネットワーク全体の辞書を決める。
################################################################################
def make_data_dic(nw_cond, base_top_list, top_dic_list):
	n_segments = nw_cond[3]
	n_sc = nw_cond[5]
	jp_xyz_dic = base_top_list[1]
	atom_jp = base_top_list[2]
	vector_dic = base_top_list[4]

	calcd_data_dic_list = []
	for str_top_dic in top_dic_list:
		atom_all = []
		pos_all = {}
		#
		strand_xyz, bond_all, atom_strand, angle_all = make_strands(str_top_dic, jp_xyz_dic, vector_dic, n_segments, n_sc)
		atom_all.extend(atom_jp)
		atom_all.extend(atom_strand)
		pos_all.update(jp_xyz_dic)
		pos_all.update(strand_xyz)
		#
		calcd_data_dic = {"atom_all":atom_all, "bond_all":bond_all, "pos_all":pos_all, "angle_all":angle_all}
		calcd_data_dic_list.append(calcd_data_dic)
	return calcd_data_dic_list

# 一本のストランド中の各アトムのxyzリストとボンドリストを作成
def make_strands(str_top_dic, jp_xyz_dic, vector_dic, n_segments, n_sc):
	strand_xyz = {}
	bond_all = {}
	atom_strand = []
	angle_all = []
	seq_atom_id = len(jp_xyz_dic)
	bond_id = 0
	# 
	for strand_id in str_top_dic:
		# 各ストランドの始点と終点のjpごとのXYZで処理する
		start_id = str_top_dic[strand_id][0]
		end_id = str_top_dic[strand_id][1]
		vector = vector_dic[strand_id]
		start_xyz = np.array(jp_xyz_dic[start_id])
		end_xyz = np.array(jp_xyz_dic[end_id])
		tmp_xyz, tmp_bond, tmp_atom_st, tmp_angle, seq_atom_id, bond_id = calc_single_strand_rnd(start_id, end_id, vector, start_xyz, end_xyz, seq_atom_id, bond_id, n_segments, n_sc)
		#
		strand_xyz.update(tmp_xyz)
		bond_all.update(tmp_bond)
		atom_strand.extend(tmp_atom_st)
		angle_all.append(tmp_angle)
	return strand_xyz, bond_all, atom_strand, angle_all

# 一本のストランド中の各アトムのxyzリストとボンドリストを作成
def calc_single_strand_rnd(start_id, end_id, vector, start_xyz, end_xyz, seq_atom_id, bond_id, n_segments, n_sc):
	tmp_xyz = {}
	tmp_bond = {}
	tmp_angle = []
	tmp_atom_st = []
	# 始点のアトムのIDを設定
	s_id = start_id
	tmp_angle.append(s_id)
	# ストランドの鎖長分のループ処理
	unit_len = 1./(n_segments + 1)
	ortho_vec = find_ortho_vec(vector)
	mod_o_vec = np.linalg.norm(vector)*ortho_vec
	for seg in range(n_segments):	
		#
		pos = tuple(start_xyz + vector*(seg + 1)*unit_len)
		tmp_xyz[seq_atom_id] = pos
		if seg == 0 or seg == n_segments - 1:
			tmp_atom_st.append([seq_atom_id, 1, 1])
		else:
			tmp_atom_st.append([seq_atom_id, 2, 2])
		e_id = seq_atom_id
		#
		if seg == 0:
			bond = 0
		else:
			bond = 1
		tmp_bond[bond_id] = tuple([bond, [s_id, e_id]])
		bond_id += 1
		tmp_angle.append(e_id)
		s_id = e_id
		seq_atom_id += 1
		#
		if n_sc != 0:
			sc_s_id = s_id
			for i in range(n_sc):
				tmp_xyz[seq_atom_id] = tuple(np.array(pos)  +  (i + 1)*mod_o_vec*unit_len)
				tmp_atom_st.append([seq_atom_id, 2, 1])
				sc_e_id = seq_atom_id
				#
				bond = 2
				tmp_bond[bond_id] = tuple([bond, [sc_s_id, sc_e_id]])
				sc_s_id = sc_e_id
				seq_atom_id += 1
				bond_id += 1
	e_id = end_id
	bond = 0
	tmp_bond[bond_id] = tuple([bond, [s_id, e_id]])
	tmp_angle.append(e_id)
	bond_id += 1
	return tmp_xyz, tmp_bond, tmp_atom_st, tmp_angle, seq_atom_id, bond_id