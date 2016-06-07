import method
import numpy as np
def var_init():

    method.N = 46  # # sites
    method.M = 28  # # v sites
    method.Ntot = 511  # total length

    method.n_Ag = 141  # # panel seq
    method.n_B = 10000  # max # BCs

    method.N_pr = 2  # # repl before sel
    method.N_re = 80  # # cycles
    method.N_pr_BB = 9  # # prol

    method.pA = 0.2  # probab of aa mut
    method.pL = 0.3  # probab of lethal mut

    method.frac_rc = 0.9  # % recycled BCs

    method.Bin_E_sel = 40

    method.sht = 0.  # shift for dX range

    method.R = 0.3  # scaling hk range

    method.n_v_site = 28

    method.n_v1_site = 11
    method.n_v2_site = 11

    method.Ec = 10.0  # threshold to avoid big weights

    method.all_MC_flag = 1

    method.f = 1.1  # C[j]*=f

    # Variables
    method.s_np = np.zeros((method.n_Ag + 1, method.Ntot + 1), dtype=int)  # full panel seq
    method.s0 = np.zeros((method.n_Ag + 1, method.N + 1), dtype=int)  # binding region

    method.s = np.zeros((6, method.N + 1), dtype=int)  # FDC-Ag seq ,s[Ag type,site]
    method.n_FDC_Ag = 0

    method.s_im_0 = np.zeros((method.N + 1,), dtype=int)  # immunogens
    method.s_im_1 = np.zeros((method.N + 1,), dtype=int)
    method.s_im_2 = np.zeros((method.N + 1,), dtype=int)

    # epitopes
    method.label_np = np.array([97, 276, 278, 279, 280, 282, 283, 365, 368, 371, 372,
                         425, 430, 455, 457, 461, 463, 465, 467, 469, 474, 134, 137, 138, 139, 140, 141, 143, 144,
                         145,
                         146, 147, 148, 149, 150, 151, 152, 460, 462, 464, 197,
                         234, 276, 386, 392, 463], dtype=int)

    method.h_GC_0 = np.zeros((5, method.N + 1), dtype=float)  # h val. for sd BCs

    method.h_GC = np.zeros((method.n_B + 1, method.N + 1), dtype=float)  # GC BCs
    method.h_sel = np.zeros((method.n_B + 1, method.N + 1), dtype=float)  # sel BCs
    method.h_m = np.zeros((method.n_B * 10 + 1, method.N + 1), dtype=float)  # Memory BCs

    method.h_s = np.zeros((method.n_B + 1, method.N + 1), dtype=float)  # surv BCs: those survived from lethal mutations

    method.avg_hk = np.zeros((method.N + 1,), dtype=float)  # # avg h val(sel.BCs)n
    method.var_hk = np.zeros((method.N + 1,), dtype=float)

    method.n_mt = np.zeros((method.n_B + 1,), dtype=int)
    # # mut s: non lethaly mutated BCs ., sel:selected, m:memory
    method.n_mt_s = np.zeros((method.n_B + 1,), dtype=int)
    method.n_mt_sel = np.zeros((method.n_B + 1,), dtype=int)
    method.n_mt_m = np.zeros((method.n_B * 10 + 1,), dtype=int)

    method.n_mt_b = np.zeros((method.n_B + 1,), dtype=int)
    method.n_mt_d = np.zeros((method.n_B + 1,), dtype=int)  # b:beneficial, d:deleterious
    method.n_mt_b_s = np.zeros((method.n_B + 1,), dtype=int)
    method.n_mt_d_s = np.zeros((method.n_B + 1,), dtype=int)
    method.n_mt_b_sel = np.zeros((method.n_B + 1,), dtype=int)
    method.n_mt_d_sel = np.zeros((method.n_B + 1,), dtype=int)
    method.n_mt_b_m = np.zeros((method.n_B * 10 + 1,), dtype=int)
    method.n_mt_d_m = np.zeros((method.n_B * 10 + 1,), dtype=int)
    method.En = np.array([1.2, 1.1, 1.05], dtype=float)  # neutralization threshold???

    method.N_A_aa = np.zeros((method.n_B + 1,), dtype=int)  # # & label of FDC Ag seen by aa-mut BCs
    method.Ag_label_aa = np.zeros((method.n_B + 1, 6), dtype=int)

    method.N_A = np.zeros((method.n_B + 1,), dtype=int)  # # & label of FDC Ag seen by survived BC( aa-mut + unmut BCs)
    method.Ag_label = np.zeros((method.n_B + 1, 6), dtype=int)

    method.Ag_label_sel = np.zeros((method.n_B + 1, 3), dtype=int)

    method.flag_aa = np.zeros((method.n_B + 1,), dtype=int)  # indicator of aa-mut BC 1: mutated
    method.v_site_label = np.array(
        [4, 11, 16, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45,
         46],
        dtype=int)

    method.v1_label = np.array([4, 11, 16, 38, 39, 40, 22, 23, 24, 25, 26], dtype=int)
    method.v2_label = np.array([27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37], dtype=int)
    method.c_site_label = np.zeros((method.N,), dtype=int)

    method.n_c_site = 0
    method.n_GC_0 = 0  # # init BCs
    method.n_GC = 0  # # curr BCs
    method.n_sel = 0
    method.n_m = 0
    method.n_s = 0

    method.hl_0 = -1. * 0.6 * 0.3  # -0.18
    method.hu_0 = 5. * 0.6 * 0.3  # 0.9

    method.hl = -1. * 0.6 * method.R  # h-range
    method.hu = 5. * 0.6 * method.R
    method.dxl = -6.5  # dX-range
    method.dxu = 1.5

    method.E0 = 60. * 0.6 * method.R * 1.0  # 10.8 kcal/mol, act.thresh.

    method.flag = np.zeros((method.n_B + 1,), dtype=int)
    method.C = np.zeros((6,), dtype=float)

    method.act_Ag = np.zeros((method.n_B + 1, method.N + 1), dtype=int)

    method.hv = np.zeros((method.n_B + 1,), dtype=float)  # interaction of sel. BCs
    method.hv1 = np.zeros((method.n_B + 1,), dtype=float)
    method.hv2 = np.zeros((method.n_B + 1,), dtype=float)
    method.hc = np.zeros((method.n_B + 1,), dtype=float)
    method.hv_s = np.zeros((method.n_B + 1,), dtype=float)  # interaction of surv. BCs
    method.hc_s = np.zeros((method.n_B + 1,), dtype=float)

    method.br = np.zeros((method.n_B + 1,), dtype=float)  # breadth
    method.n_ntAg = np.zeros((method.n_B + 1,), dtype=int)  # # neut Ags out of panels

    method.n = np.zeros((method.Bin_E_sel + 1,), dtype=int)  # # bins for Energy histogram
    method.P = np.zeros((method.Bin_E_sel + 1,), dtype=float)  # E_sel hist
    method.Ej = np.zeros((method.n_B + 1,), dtype=float)  # E sel BC. val

    method.max_mut = np.zeros((method.n_B + 1,), dtype=int)  # distinct mut
    method.max_mut_s = np.zeros((method.n_B + 1,), dtype=int)
    method.max_mut_sel = np.zeros((method.n_B + 1,), dtype=int)

    method.mut_site = np.zeros((method.n_B + 1, method.N + 1), dtype=int)  # mutated site of B cell?
    method.mut_site_s = np.zeros((method.n_B + 1, method.N + 1), dtype=int)  # mutated site of survived B cell?
    method.mut_site_sel = np.zeros((method.n_B + 1, method.N + 1), dtype=int)  # mutated site of selected B cell?

    method.W = np.zeros((method.n_B * 10 + 1,), dtype=float)  # weight of all competing cells
    method.actAg_MC = np.zeros((method.n_B * 10 + 1,), dtype=int)  # Ag seen by memorgy B cells, randomly assigned

    method.E_sel_min = 100.
    method.E_sel_max = -10.
    method.E_sel = np.zeros((method.n_B + 1,), dtype=float)

    method.h = np.zeros((method.N + 1,), dtype=float)
    method.E_ij = np.zeros((6,), dtype=float)  # never used
    # # GCR cycle (<80)
    method.avg_Ps_s = 0.  # avg sel probab she doesn't even use this

    method.seed = 100  # used for rd # generation
    method.n_round = 0  # # of rounds of complete GC cycle on 1 instantiation of GC lp label
    method.t_tot = 0  # tot_# of GCRS
    method.flag_succ = 0  # t>1&&n_GC>=n_GC_0*pow(2,N_pr_BB))||t==N_re
    method.flag_end = 0  # no modification
    method.flag_dX = 0  # mut type 1:affinity +  -1:affinity -   2:shielding motif
    method.t = 0

