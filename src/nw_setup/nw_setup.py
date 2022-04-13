#!/usr/bin/env python
# -*- coding: utf-8 -*-
################################################################
import nw_setup
################################################################
################################################################
def main():
    ################################################################
    # 設定条件を読み込み、ネットワークポリマーの諸量を計算
    basic_cond, nw_cond, sim_cond, sim_cond2, rnd_cond, target_cond, target_dir = nw_setup.setupcondition()

    ###################
    # ネットワークを設定
    calcd_data_dic = nw_setup.select_set(nw_cond, sim_cond2, rnd_cond, target_dir)

    ##################
    # baseUDF の作成
    baseudf = nw_setup.MakeInitUDF(basic_cond, sim_cond, sim_cond2, target_cond, calcd_data_dic)
    baseudf.setup_baseudf(target_dir)
    
    ###############
    # シミュレーションを設定
    setup = nw_setup.SetUpUDF(basic_cond, sim_cond, sim_cond2, target_dir)
    setup.setup_udf()

################################################################################
#      Main     #
################################################################################
if __name__=='__main__':
	main()