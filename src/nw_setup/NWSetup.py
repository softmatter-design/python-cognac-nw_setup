#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################################
from ast import arg
import numpy as np
import copy
import random
import platform
import subprocess
import sys
import os
import pickle
from multiprocessing import Pool
from UDFManager import UDFManager
################################################################################
class SelectSet:
	def __init__(self, nw_cond, sim_cond, rnd_cond, target_dir):
		self.nw_model = nw_cond[0]
		self.strand = nw_cond[1]
		self.n_strand = nw_cond[2]
		self.n_segments = nw_cond[3]
		self.n_cell = nw_cond[4]
		self.n_sc = nw_cond[5]

		self.multi = sim_cond[1]

		self.restart = rnd_cond[0] 
		self.cond_top = rnd_cond[1]
		self.n_hist = rnd_cond[2]

		self.target_dir = target_dir
	############################
	def select_set(self):
		# ネットワークを設定
		if self.nw_model == "Regular":
			calcd_data_dic = self.regnw_setup()
		elif self.nw_model == "Random":
			calcd_data_dic = self.rndnw_setup()

		return calcd_data_dic
	
	#######################################
	def regnw_setup(self):
		calcd_data_dic = self.calc_all()
		return calcd_data_dic

	def rndnw_setup(self):
		self.base_top_list = self.make_8chain_dic()

		if self.restart == '':
			# トポロジーの異なるネットワークを探索して、任意の多重度のネットワークポリマーの代数的連結性の分布関数を策定
			candidate_list, random_dir = self.top_search()
		else:
			# 設定したリスタートファイルを読み込んで、リストを作成
			candidate_list, random_dir = self.top_select()

		top_dic_list = self.nw_search(candidate_list, random_dir)

		###########################################
		# ターゲットとなるネットワーク全体の辞書を設定。
		calcd_data_dic = self.make_data_dic(top_dic_list)
		return calcd_data_dic


	################################################
	## REGULAR NW SETUP
	################################################
	def calc_all(self):
		# 架橋点 JP を設定
		jp_xyz, strand_se_xyz = self.calc_jp_strands()
		#
		calcd_data_dic = self.set_atom(jp_xyz, strand_se_xyz)
		# print(len(calcd_data_dic))
		return calcd_data_dic
	##########################################
	# JPおよびサブチェインの始点と終点のXYZを設定
	def calc_jp_strands(self):
		# jp_xyz は、JPの座標のリスト
		# strand_se_xyz は、サブチェインの出発点と終点のリスト
		if self.strand == "3_Chain_S":
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

		elif self.strand == "3_Chain_D":
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

		elif self.strand == "4_Chain":
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

		elif self.strand == "6_Chain":
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

		elif self.strand == "8_Chain":
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
	def set_atom(self, jp_xyz, strand_se_xyz):
		calcd_data_dic={}
		for mul in (range(self.multi)):
			# atom_all = []
			# pos_all = {}
			for mol, jp in enumerate(jp_xyz):
				atom_all = []
				pos_all = {}
				# システム全体にわたるジャンクションポイントのxyzとIDの辞書を作成
				jp_id_dic, jp_xyz_dic, atom_jp = self.set_jp_id_reg(jp)
				atom_all.extend(atom_jp)
				pos_all.update(jp_xyz_dic)
				# print(jp_xyz_dic)
				# サブチェイン中の各アトムのxyzリストとボンドリストを作成
				strand_xyz, bond_all, atom_sc, angle_all = self.set_strands(jp_id_dic, strand_se_xyz[mol], mol)
				#
				atom_all.extend(atom_sc)
				pos_all.update(strand_xyz)
				#
				calcd_data_dic[mul+mol] = {"atom_all":atom_all, "bond_all":bond_all, "pos_all":pos_all, "angle_all":angle_all}
		return calcd_data_dic

	###################################################
	# システム全体にわたるJPのxyzとIDの辞書を作成
	def set_jp_id_reg(self, jp_xyz):
		jp_id_dic = {}
		jp_xyz_dic = {}
		atom_jp = []
		jp_id = 0
		for z in range(self.n_cell):
			for y in range(self.n_cell):
				for x in range(self.n_cell):
					base_xyz = np.array([x,y,z])
					for jp in jp_xyz:
						jp_id_dic[tuple(np.array(jp) + base_xyz)] = (jp_id)
						jp_xyz_dic[(jp_id)] = tuple(np.array(jp) + base_xyz)
						atom_jp.append([jp_id, 0, 0])
						jp_id += 1
		return jp_id_dic, jp_xyz_dic, atom_jp

	#########################################################
	# サブチェイン中の各アトムのxyzリストとボンドリストを作成
	def set_strands(self, jp_id_dic, strand_se_xyz, mol):
		strand_xyz = {}
		bond_all = {}
		atom_sc = []
		angle_all = []
		sub_id = len(jp_id_dic)
		bond_id = 0
		for z in range(self.n_cell):
			for y in range(self.n_cell):
				for x in range(self.n_cell):
					b_xyz = (x,y,z)
					for se_xyz in strand_se_xyz:
						tmp_xyz, tmp_bond, new_sub_id, new_bond_id, tmp_atom_sc, tmp_angle = self.calc_single_strand_reg(jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, mol)
						strand_xyz.update(tmp_xyz)
						bond_all.update(tmp_bond)
						atom_sc.extend(tmp_atom_sc)
						angle_all.append(tmp_angle)
						sub_id = new_sub_id
						bond_id = new_bond_id
		return strand_xyz, bond_all, atom_sc, angle_all

	###############################################################
	# 一本のサブチェイン中の各アトムのxyzリストとボンドリストを作成
	def calc_single_strand_reg(self, jp_id_dic, sub_id, bond_id, b_xyz, se_xyz, mol):
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
		unit_len = 1./(self.n_segments + 1)
		ortho_vec = self.find_ortho_vec(vec)
		mod_o_vec = np.linalg.norm(vec)*ortho_vec
		# 始点のアトムのIDを設定
		mod_xyz = list(start_xyz)[:]
		for dim in range(3):
			if mod_xyz[dim] == self.n_cell:
				mod_xyz[dim] = 0
		s_id = jp_id_dic[tuple(mod_xyz)]
		tmp_angle.append(s_id)
		# 終点のアトムのIDを周期境界条件で変更
		mod_xyz = list(end_xyz)[:]
		for dim in range(3):
			if mod_xyz[dim] == self.n_cell:
				mod_xyz[dim] = 0
		E_id = jp_id_dic[tuple(mod_xyz)]
		# サブチェインの鎖長分のループ処理
		for seg in range(self.n_segments):
			pos = tuple(start_xyz + vec*(seg + 1)/(self.n_segments + 1.))
			tmp_xyz[sub_id] = pos

			if seg == 0 or seg == self.n_segments - 1:
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
			
			if self.n_sc != 0:
				sc_s_id = s_id
				for i in range(self.n_sc):
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
	def find_ortho_vec(self, list):
		vec = np.array(list).reshape(-1,1)
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
	def make_8chain_dic(self):
		# ユニットセルでの、jp およびサブチェインの始点と終点のXYZを設定
		self.jp_xyz = [
				[0.,0.,0.], 
				[0.5, 0.5, 0.5]
				]
		self.strand_se_xyz = [
					[[0.5, 0.5, 0.5], [0, 0, 0]],
					[[0.5, 0.5, 0.5], [1, 0, 0]],
					[[0.5, 0.5, 0.5], [0, 1, 0]],
					[[0.5, 0.5, 0.5], [1, 1, 0]],
					[[0.5, 0.5, 0.5], [0, 0, 1]],
					[[0.5, 0.5, 0.5], [1, 0, 1]],
					[[0.5, 0.5, 0.5], [0, 1, 1]],
					[[0.5, 0.5, 0.5], [1, 1, 1]]
					]
		self.start_jp = [0.5, 0.5, 0.5]
		# システム全体にわたるピボットのxyzとIDの辞書を作成
		jp_id_dic, jp_xyz_dic, atom_jp, n_jp = self.set_jp_id_rnd()
		# ストランドの結合状態を記述
		init_8ch_dic, vector_dic = self.set_strands_8(jp_id_dic)
		#
		base_top_list = [init_8ch_dic, jp_xyz_dic, atom_jp, n_jp, vector_dic]

		return base_top_list 

	##########################################
	# システム全体にわたるjpのxyzとIDの辞書を作成
	def set_jp_id_rnd(self):
		jp_id = 0
		jp_id_dic = {}
		jp_xyz_dic = {}
		atom_jp = []
		for z in range(self.n_cell):
			for y in range(self.n_cell):
				for x in range(self.n_cell):
					base_xyz = np.array([x,y,z])
					for jp in self.jp_xyz:
						jp_id_dic[tuple(np.array(jp) + base_xyz)] = (jp_id)
						jp_xyz_dic[jp_id] = tuple(np.array(jp) + base_xyz)
						atom_jp.append([jp_id, 0, 0])
						jp_id += 1
		n_jp = jp_id
		return jp_id_dic, jp_xyz_dic, atom_jp, n_jp

	#####################################################
	# ストランドの結合状態を記述
	def set_strands_8(self, jp_id_dic):
		init_8ch_dic = {}
		vector_dic = {}
		str_id = 0
		for z in range(self.n_cell):
			for y in range(self.n_cell):
				for x in range(self.n_cell):
					base_xyz = np.array([x,y,z])
					for xyz in self.strand_se_xyz:
						# サブチェインの末端間のベクトルを設定
						start_xyz = np.array(self.start_jp) + base_xyz
						end_xyz = np.array(xyz[1]) + base_xyz
						vector = end_xyz - start_xyz
						# 始点のアトムのIDを設定
						start_id = jp_id_dic[tuple(start_xyz)]
						# 終点のアトムのIDを周期境界条件で変更
						mod_end_xyz = list(end_xyz)[:]
						for i in range(3):
							if mod_end_xyz[i] == self.n_cell:
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
	def top_search(self):
		self.init_dic = self.base_top_list[0]
		self.n_jp = self.base_top_list[3]
		#
		self.pre_sampling = self.cond_top[0] 
		self.n_sampling = self.cond_top[1]
		self.n_try = self.cond_top[2]
		self.repeat = self.cond_top[3]
		self.f_pool = self.cond_top[4]
		#
		random_dir = str(self.n_strand) +"_chains_" + str(self.n_cell) + "_cells_"
		random_dir += str(self.n_sampling) + "_sampling_" + str(self.n_try) + "_trials_" + str(self.repeat) + "_times"
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
		candidate_list = self.strand_exchange()
		#
		with open(os.path.join(random_dir, 'init.pickle'), mode = 'wb') as f:
			pickle.dump(candidate_list, f)
		#
		print("##################################################")
		print("Trial, Sampling = ", self.n_try, ", ", self.n_sampling)
		print("Total Sampling = ", self.n_try*self.n_sampling)
		print("Initial Candidates = ", len(candidate_list))
		print("##################################################")

		return candidate_list, random_dir

	#####################################################
	# 任意のストランドを選択し、ストランドの繋ぎ変えを行う
	def strand_exchange(self):
		tmp_list = []
		p = Pool(self.f_pool)
		# "self.random_reduce()" により任意のストランドを除去したものをサンプリング
		print("Searching for Initial Random Structure of ", self.pre_sampling)
		print("Please Wait a little bit!")
		args = [[self.random_reduce(), i, "Pre"] for i in range(self.pre_sampling)] 
		result = p.map(self.search, args)
		for data in result:
			tmp_list.extend(data)
		print("##################################################")
		print("Pre_Search Result:", len(tmp_list))
		print("##################################################")
		#
		for repeat in range(self.repeat):
			# tmp_list を alg_const の昇順に並べ替えて、その分布に応じてサンプリング
			tmp_array = np.array(tmp_list)
			sortbyalg = tmp_array[np.argsort(tmp_array[:,1])]
			args = []
			step = int(len(sortbyalg)/self.n_sampling)
			for i in range(self.n_sampling):
				args.append([list(sortbyalg[i*step])[0], i, str(repeat + 1)+"/"+str(self.repeat)])
			# 
			result = p.map(self.search, args)
			tmp_list = []
			for data in result:
				tmp_list.extend(data)
			print("##################################################")
			print(str(repeat)+"/"+str(self.repeat), "Search Result:", len(tmp_list))
			print("##################################################")

		return tmp_list

	########################################################
	# ストランドの繋ぎ変えを行う
	def search(self, args):
		dic, x, id = args
		alg_const = self.calc_lap_mat(dic)
		print(id, "Sampling ID =", x,  "Initial Algebratic Connectivity =", alg_const)
		#
		result = []
		count = 0
		failed = 0
		show = 5
		while count < self.n_try:
			# 現状のストランドのリストの中からランダムに一つ選択し、"selected_strand"とする。
			selected_strand = random.choice(list(dic.keys()))
			# 繋ぎ変え得るストランドのリストを現状のネットワークと比較して、交換可能なセットを見つける。
			tmp_dic, alg_const = self.find_pair(selected_strand, dic)
			if alg_const != 0:
				count += 1
				result.append([tmp_dic, alg_const])
				failed = 0
				if count != 0 and round(show*count/self.n_try) == show*count//self.n_try and round(show*count/self.n_try) == -(-show*count//self.n_try):
					print(id, "Sampling ID =", x, "count = ", count)
			else:
				failed +=1
			if failed >= self.n_try:
				print("##########################################")
				print(id, "Sampling ID =", x,  " FAILED!! with ", failed, "th trials.")
				print("##########################################")
				count = self.n_try
				failed = 0
		# 
		return result

	#########################################################
	# 任意のストランドを選択し、所望の分岐数にストランドを消去
	def random_reduce(self):
		alg_const = 0
		while alg_const <= 0:
			tmp_dic = copy.deepcopy(self.init_dic)
			del_bond = []
			# 消去対象のストランドをリストアップ
			# ユニットセル中のストランドから必要な数だけ抽出
			tmp_del = random.sample(range(8), (8 - self.n_strand))
			# 全セルに渡って、消去するストランドのリストを作成
			for i in range(self.n_cell**3):
				del_bond.extend(list(8*i + np.array(tmp_del)))
			# ストランドを消去した辞書とその代数的連結性を得る。
			for target in del_bond:
				tmp_dic.pop(target)
			alg_const = self.calc_lap_mat(tmp_dic)
		return tmp_dic

	################################################################
	# 交換可能なストランドのペアを見つける
	def find_pair(self, selected_strand, dic):
		deleted_dic = dict(self.init_dic.items() - dic.items())
		# 選択したストランドの両端のjp（"jp_pairs"）を見つけ、
		# それと繋がり得る可能性のあるjpをリストアップ（"connected_list"）
		# 繋ぎ変え得るストランドを見つける（"possible_strand"）
		possible_strand = self.get_connected_jp(self.init_dic, selected_strand)
		# 繋ぎ変え得るストランドのリストを現状のネットワークと比較して、交換可能なセットを見つける。
		found_dic = {}
		for p_str in possible_strand:
			if p_str in dic:
				# 現行のストランドリストの中にいないものを選択
				count = 0
				tmp_dic = {}
				for tt in dic[selected_strand]:
					for uu in dic[p_str]:
						tmp_str = []
						if tt != uu:
							if (tt,uu) in deleted_dic.values():
								tmp_str = [k for k, v in self.init_dic.items() if v == (tt, uu)]
							elif (uu,tt) in deleted_dic.values():
								tmp_str = [k for k, v in self.init_dic.items() if v == (uu, tt)]
							if tmp_str != []:
								tmp_dic[tmp_str[0]] = tuple(self.init_dic[tmp_str[0]])
								count +=1
				if count == 2 and tmp_dic != {}:
					found_dic.update(tmp_dic)
					dic.pop(p_str)
					dic.pop(selected_strand)
					dic.update(tmp_dic)
					alg_const = self.calc_lap_mat(dic)
					return dic, alg_const
				else:
					alg_const = 0
			else:
					alg_const = 0
		return dic, alg_const

	#################################################################
	# 任意のjpと連結したjpを見つけ、繋がり得る可能性のあるjpをリストアップ
	def get_connected_jp(self, dic, selected_strand):
		# 選択したストランドの両端のjp（"jp_pairs"）を見つけ、
		jp_pairs = tuple(dic[selected_strand])
		# ラプラシアン行列を利用して、それと繋がり得る可能性のあるjpをリストアップ（"connected_list"）
		lap_mat = self.make_lap_mat(dic)
		connected_list = []
		for t_jp in jp_pairs:
			t_list = list(lap_mat[t_jp])
			connected_list.append(self.find_index(t_list))
		# 繋ぎ変え得るストランドを見つける。
		possible_strand = []
		for i in connected_list[0]:
			for j in connected_list[1]:
				link = self.get_link(dic, i, j)
				if link:
					# print(link)
					possible_strand.append(link[0])
		return possible_strand

	############################################################################
	# ラプラシアン行列を利用して、任意のjpとストランドによりつながる８個のjpを見つける。
	def find_index(self, t_list):
		con_list = [i for i, x in enumerate(t_list) if x == -1.0]
		return con_list

	#######################################################
	# 任意の二つのjpがストランドによりつながっている可能性を調査
	# 繋がり得るものであれば、そのストランドのID番号を返す。
	def get_link(self, dic, jp0, jp1):
		val_tpl = tuple([jp0, jp1])
		# print("val", val_tpl)
		link = [k for k, v in dic.items() if v == val_tpl]
		# print(link)
		if link:
			return link
		elif link == []:
			val_tpl = tuple([jp1, jp0])
			link = [k for k, v in dic.items() if v == val_tpl]
			if link:
				return link
		else:
			return []

	####################################################
	# 任意のネットワークトポロジーから、ラプラシアン行列を作成
	def calc_lap_mat(self, topl_dic):
		lap_mat = self.make_lap_mat(topl_dic)
		# 固有値を計算
		la, v = np.linalg.eig(lap_mat)
		# 固有値の二番目の値を見つける。
		alg_const = round(sorted(np.real(la))[1], 3)
		return alg_const

	####################################################
	# 任意のネットワークトポロジーから、ラプラシアン行列を作成
	def make_lap_mat(self, topl_dic):
		lap_mat = np.zeros((self.n_jp, self.n_jp))
		for i in topl_dic:
			lap_mat[topl_dic[i][0], topl_dic[i][1]] = -1
			lap_mat[topl_dic[i][1], topl_dic[i][0]] = -1
			lap_mat[topl_dic[i][0], topl_dic[i][0]] += 1
			lap_mat[topl_dic[i][1], topl_dic[i][1]] += 1
		return lap_mat

	#########################################################
	# 過去の探索データを使って、代数的連結性の分布関数を選択
	def top_select(self):
		with open(os.path.join(self.restart, 'init.pickle'), mode = 'rb') as f:
			candidate_list = pickle.load(f)
		print("##################################################")
		print("Reloaded Candidates = ", len(candidate_list))
		print("##################################################")

		return candidate_list, self.restart


	#####################################################################
	# ヒストグラム中の最大頻度を与えるネットワークトポロジーの配列辞書を決める。
	def nw_search(self, candidate_list, random_dir):
		tmp_list = []
		val_list = []
		data_list = []

		# ヒストグラムを作成
		histdata = list(np.array(candidate_list)[:, 1])
		cond = ["", histdata, self.n_hist, "True", ['Arg. Con.', 'Freq.']]
		x, val = self.make_hist_all(cond, random_dir)

		# 最頻値のレンジを決める
		val_range = self.find_range(x, val)
		# 上記のレンジに入る配列をピックアップ
		for i in candidate_list:
			if i[1] >= val_range[0] and i[1] <= val_range[1]:
				tmp_list.append(i)
		random.shuffle(tmp_list)
		# 
		count = 0
		for i, selected_list in enumerate(tmp_list):
			# print(i, count, selected_list[0], val_list)
			if selected_list[0] not in val_list:
				val_list.append(selected_list[1])
				data_list.append(selected_list[0])
				count += 1
			if len(val_list) == self.multi:
				u = UDFManager(os.path.join(self.target_dir, 'target_condition.udf'))
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
		print("No effective list was found for multi numbers of", self.multi, "!  Try again!!")
		sys.exit()


	#######################
	# 最大頻度の範囲を決める。
	def find_range(self, x, val):
		index = np.where(val == max(val))[0][0]
		f_range = 0
		value = 0
		if index < len(val) -1:
			while value < self.multi:
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
	def make_hist_all(self, cond_list, random_dir):
		self.base = cond_list[0]
		self.list = cond_list[1]
		self.bins = cond_list[2]
		self.norm = cond_list[3]
		self.leg = cond_list[4]

		self.f_dat = "nw_hist.dat"
		self.f_plt = "make_hist.plt"
		self.f_png = "histgram.png"
		# ヒストグラムのデータ作成
		bin_width, hist_data, val, x = self.make_hist_data()
		# ヒストグラムのデータを書き出し 
		self.write_data(hist_data, random_dir)
		# グラフを作成
		self.make_graph(bin_width, random_dir)
		return x, val

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
		return bin_width, hist_data, val, x

	# ヒストグラムのデータを書き出し 
	def write_data(self, hist_data, random_dir):
		os.makedirs(random_dir, exist_ok=True)
		with open(os.path.join(random_dir, self.f_dat), 'w') as f:
			f.write("# Histgram data:\n\n")
			for line in hist_data:
				f.write(str(line[0]) + '\t' + str(line[1])  + '\n')
		return

	# グラフを作成
	def make_graph(self, bin_width, random_dir):
		self.make_script(bin_width, random_dir)
		cwd = os.getcwd()
		os.chdir(random_dir)
		if platform.system() == "Windows":
			subprocess.call(self.f_plt, shell=True)
		elif platform.system() == "Linux":
			subprocess.call('gnuplot ' + self.f_plt, shell=True)
		os.chdir(cwd)
		return
		
	# 必要なスクリプトを作成
	def make_script(self, bin_width, random_dir):
		script = 'set term pngcairo font "Arial,14" \nset colorsequence classic \n'
		script += '# \ndata = "' + self.f_dat + '" \nset output "' + self.f_png + ' "\n'
		script += '#\nset size square\n# set xrange [0:1.0]\n#set yrange [0:100]\n'
		script += '#\nset xlabel "' + self.leg[0] + '"\nset ylabel "' + self.leg[1] + '"\n\n'
		script += 'set style fill solid 0.5\nset boxwidth ' + str(bin_width) + '\n'
		script += '#\nplot data w boxes noti'
		print("OK")
		with open(os.path.join(random_dir, self.f_plt), 'w') as f:
			# script = self.script_content(bin_width)
			f.write(script)
		return
		
	################################################################################
	# ターゲットとなるネットワーク全体の辞書を決める。
	################################################################################
	def make_data_dic(self, top_dic_list):
		self.jp_xyz_dic = self.base_top_list[1]
		self.atom_jp = self.base_top_list[2]
		self.vector_dic = self.base_top_list[4]

		calcd_data_dic_list = []
		for str_top_dic in top_dic_list:
			atom_all = []
			pos_all = {}
			#
			strand_xyz, bond_all, atom_strand, angle_all = self.make_strands(str_top_dic)
			atom_all.extend(self.atom_jp)
			atom_all.extend(atom_strand)
			pos_all.update(self.jp_xyz_dic)
			pos_all.update(strand_xyz)
			#
			calcd_data_dic = {"atom_all":atom_all, "bond_all":bond_all, "pos_all":pos_all, "angle_all":angle_all}
			calcd_data_dic_list.append(calcd_data_dic)
		return calcd_data_dic_list

	# 一本のストランド中の各アトムのxyzリストとボンドリストを作成
	def make_strands(self, str_top_dic):
		strand_xyz = {}
		bond_all = {}
		atom_strand = []
		angle_all = []
		seq_atom_id = len(self.jp_xyz_dic)
		bond_id = 0
		# 
		for strand_id in str_top_dic:
			# 各ストランドの始点と終点のjpごとのXYZで処理する
			start_id = str_top_dic[strand_id][0]
			end_id = str_top_dic[strand_id][1]
			vector = self.vector_dic[strand_id]
			start_xyz = np.array(self.jp_xyz_dic[start_id])
			end_xyz = np.array(self.jp_xyz_dic[end_id])
			tmp_xyz, tmp_bond, tmp_atom_st, tmp_angle, seq_atom_id, bond_id = self.calc_single_strand_rnd(start_id, end_id, vector, start_xyz, end_xyz, seq_atom_id, bond_id)
			#
			strand_xyz.update(tmp_xyz)
			bond_all.update(tmp_bond)
			atom_strand.extend(tmp_atom_st)
			angle_all.append(tmp_angle)
		return strand_xyz, bond_all, atom_strand, angle_all

	# 一本のストランド中の各アトムのxyzリストとボンドリストを作成
	def calc_single_strand_rnd(self, start_id, end_id, vector, start_xyz, end_xyz, seq_atom_id, bond_id):
		tmp_xyz = {}
		tmp_bond = {}
		tmp_angle = []
		tmp_atom_st = []
		# 始点のアトムのIDを設定
		s_id = start_id
		tmp_angle.append(s_id)
		# ストランドの鎖長分のループ処理
		unit_len = 1./(self.n_segments + 1)
		ortho_vec = self.find_ortho_vec(vector)
		mod_o_vec = np.linalg.norm(vector)*ortho_vec
		for seg in range(self.n_segments):	
			#
			pos = tuple(start_xyz + vector*(seg + 1)*unit_len)
			tmp_xyz[seq_atom_id] = pos
			if seg == 0 or seg == self.n_segments - 1:
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
			if self.n_sc != 0:
				sc_s_id = s_id
				for i in range(self.n_sc):
					tmp_xyz[seq_atom_id] = tuple(np.array(pos)  +  (i + 1)*mod_o_vec*unit_len)
					# print(tuple(np.array(pos)  +  (i + 1)*ortho_vec*unit_len))
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
	#
	def find_ortho_vec(self, list):
		vec = np.array(list).reshape(-1,1)
		# 線形独立である新たな三次元ベクトルを見つける。
		rank = 0
		while rank != 2:
			a = np.array(np.random.rand(3)).reshape(-1,1)
			target = np.hstack((vec, a))
			rank = np.linalg.matrix_rank(target)
		# QR分解により
		q, r = np.linalg.qr( target )
		# print(q[:,1])
		ortho_vec = q[:,1]
		return ortho_vec


