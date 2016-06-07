import math
import sys

import GCR_initiation
import differentiation_recycle
import method
import output_reset
import proliferation_SHM
import removeLethal
import selection_mixture_encounter
import var_init


#include <stdio.method.h>              # This program simulates dynamics of affinity maturation in germinal centers
#include <stdlib.method.h>             # implemented in its entirety by Shenshen Wang on June 15, 2015
#include "time.method.h"	 	        # for seeding the rd # generator

# # Functions
#            # for SHM             # weights of competing cells
#
# #1. 3 Activated B blasts that have exceeded threshold, and their energies(Ei),# blasts in initiated GC, for each blasts
# # method.h for v sites and c sites, # B blasts after 9 cycle of proliferation before SHM kicks in
# #2. GC cycle and population just after initiation,for each proliferation print method.h values for 1st cell& n_Gc+1 cell (Should have the same hs)
# ##GC after each prol. and mut/BC after SHM
# #method.n_s, method.n_sel, method.n_m: # surv BC, # sel BC, # memorgy BC @ the end of simulationn_sel, avg_hv, avg_hc: # sel BCs, method.hv,method.hc averaged over sel BCs & corresponding sites
# # nth GC cycle, Ag type, nth bin in which sel BC belongs to, Eng of that sel BC, #sel BCs in the bin, method.n_mt_sel,nth bin in which belongs to,
# # middle pt Eng in the bin,#sel BCs in bin, % out of total sel BCs
# #nth GC cycle,# sel BCs,total_b_mut in a clone,total_d_mut_in_a_clone,b_mtperclone,d_mtperclone,nt thresh.,ratio,avg_br[l], stdev_br[l], avg_#_ntAg[l], stdev_#_ntAg[l], max_#_ntAg[l], min_#_ntAg[l]         # # cells in diff cpt
# #method.t_tot,method.t,n_CB,avg_mt,avg_b_mt,avg_d_mt,max_mt, min_mt, method.n_m
# #method.t_tot,method.t,# BCS after 9 rds of prof,(n_mut_tot/method.n_sel),(n_mut_b_tot/method.n_sel),(n_mut_d_tot/method.n_sel),avg_#_distinct_mut/sel.BC, max_n_mt in a pool of BC, min_n_mt in a pool of BC, # memory BCs/
# ##
# #method.t,avg_hc_s, avg_hv_s, avg_hc, avg_hv, avg_hv1, avg_hv2, sqrt(var_hc), sqrt(var_hv), sqrt(var_hv1), sqrt(var_hv2)
# #*method.t, Csites method.h val. avg. over surv.BCs& # c sites,Vsites method.h val. avg. over surv.BCs& # v sites,method.C sites method.h val. avg. over sel.BCs& # c sites,V sites method.h val. avg. over sel.BCs& # v sites,
# #,v1 sites method.h val. avg. over sel.BCs& # v1 sites,v2 sites method.h val. avg. over sel.BCs& # v2 sites,sqrt(var_hc),sqrt(var_hv), sqrt(var_hv1), sqrt(var_hv2)
# #only here, no file assigned
# #each of Ag panel sequence (46sites))
# #store labels of c&v sites, for each affinity affecting mutated BCs store #Ag variants in FDC, BC tag, #Ag seen by BC
# # after selection process store method.N_A_aa of survived BC, method.n_s,method.N_A_aa and Ag_lable for corresponding s_BC,
# #method.t, nt thresh, ratio, avg_br, sqrt(var_br), avg_# neutralized Ag panel per sel. BCs, sqrt(var_n_ntAg), Max # ntd Ag per among BCs/cycle, Max # ntd Ag per among BCs/cycle
# #Ag panel full sequence, only to be read            # Em[j] & (method.E_ij[j]-Em[j])
# #only created but output has assinged
# #method.t_tot,method.t, # avg_br, sqrt(var_br), method.flag_succ
# #method.t_tot,method.t,# & sites avg.h_v ,stdev_h,method.flag_succ
# #method.t_tot,method.t,# & sites avg.h_v ,stdev_h,method.flag_succ           # success indicator & duration
# #method.t_tot, method.t, method.flag_succ, only store if it has to exit GC( with method.flag_succ =1)        # # rc BCs that can neut WT, v1, v2
# #method.t_tot, method.t, method.n_sel, n_nt_WT, n_nt_v1, n_nt_v2, avg_E_WT, avg_E_v1, avg_E_v2
# # method.t_tot,method.t # sel BCs,# sel. BCs neutralizing Ag WT,# sel. BCs neutralizing Ag V1 ,# sel. BCs neutralizing Ag V2,avg_E_WT, avg_E_v1, avg_E_v2
# ##method.t_tot, method.t, E_dyn:dynamic threshold
# #method.t_tot, method.t, number avg. of method.h stored in the order of tp 1 c -> v1 -> v2 -> tp 2 method.C for every 5 rds or method.t==1
# #method.t_tot, method.t, number stdev of method.h stored in the order of tp 1 c -> v1 -> v2 -> tp 2 method.C for every 5 rds or method.t==1
# ##method.t_tot, method.t, site avg of c1,c2,v1,v2          # sel probab
# ##method.t_tot, method.t, ni, method.n_s,n_a, P_a, P_Th, Psel
# method.t_tot, method.t, # sel BCs , surv.BCs, Ag_sel.BCs, Avg P_a, avg_P_Th, avg_P_a*avg_P_Th
# #method.t_tot, method.t, ith survived cells which got selected, Eij binding E, Eij-Eo from threshold,Prob. sel by Ag,Prob.sel by Th
def main():
    num_nan_error = 0
    for i in range(100):
        sim_ID = i
        var_init.var_init()
        try:
            run_sim(sim_ID)
        except ZeroDivisionError:
            num_nan_error +=1
            print 'ZeroDivisionError'
            sys.exc_clear()
    Sim_summary = open("simulation_summary.txt","w")
    Sim_summary.write("Num_Run=%d,Num_ZeroDivisionError= %d" %(100, num_nan_error))
    Sim_summary.close()

def run_sim(sim_ID):

    panel_full=open("gp120_panel_full.txt", "rb")
    init_BC=open("init_BC_N46_M28.txt", "w")
    E_sel_Hist=open("E_sel_Hist.txt", "w")
    sel_BC=open("sel_BC.txt", "w")
    n_mut=open("n_mut_breadth.txt", "w")
    output_ns=open("ns_vs_t.txt", "w")
    avg_hc_hv=open("avg_hc_hv.txt", "w")
    nt_seq=open("nt_seq.txt", "w")
    c_v_site=open("c_v_site.txt", "w")
    nt_br=open("nt_br.txt", "w")

    mean_E=open("mean_E.txt", "w")
    br_t=open("br_t.txt", "w")
    hc_t=open("hc_t.txt", "w")
    hv_t=open("hv_t.txt", "w")
    succ_ts=open("succ_ts.txt", "w")
    nt_WT_v1_v2=open("nt_WT_v1_v2.txt", "w")
    meanE_vs_t=open("meanE_vs_t.txt", "w")
    avg_hk_t=open("avg_hk_t.txt", "w")
    std_hk_t=open("std_hk_t.txt", "w")
    avg_h_c_v_vs_t=open("h_traj.txt", "w")
    avg_Psel=open("avg_Psel.txt", "w")
    sel_scatter=open("sel_scatter.txt", "w")


    output_ns.write( "# method.t_tot\tmethod.t\tn_CB\tmethod.n_s\tmethod.n_sel\tn_CC\tnr\tmethod.C[Ag]\tmut_t(mut/selB)\tmut_b(benmut/selB)\tmut_d(delmut/selB)\tavg_mt\t max_mt\tmin_mt\tmethod.n_m(memB)\n")
    avg_hc_hv.write("# method.t\t method.h[4]\t method.h[11]\t method.h[16]\t method.h[38]\t method.h[39]\t method.h[40]\t method.hc_s(avg_hc_s)\t method.hv_s(avg_hv_s)\t hc_sel(avg_hc)\t hv_sel(avg_hv)\t method.hv1(avg_hv1)\t method.hv2(avg_hv2)\t sd_hc(stdev_hc)\t sd_hv(stdev_hv)\t sd_hv1(stdev_hv1)\t sd_hv2(stdev_hv2)\n")
    nt_br.write("# method.t\t method.En(Neut.aff.thr)\t ratio\t avg_br\t stdev_br\t avg_n_ntAg(Avg.#.Panel_seq_neut)\t stdev_n_ntAg\t max_n_ntAg(Max#Panel_seq.neut)\t min_n_ntAg(Min#Panel_seq.neut)\n")
    br_t.write("#method.t_tot\t method.t\t avg_br\t std_br\t method.flag_succ\n")
    hc_t.write("#method.t_tot\t method.t\t avg_hc\t std_hc\t method.flag_succ\n")
    hv_t.write("#method.t_tot\t method.t\t avg_hv\t std_hv\t method.flag_succ\n")
    succ_ts.write("#ts\t method.t\t method.flag_succ\n")
    nt_WT_v1_v2.write( "#method.t_tot\t method.t\t method.n_sel\t n_nt_WT(#selB_neut.WT))\t n_nt_v1(#selB_neut.v1)\t n_nt_v2(#selB_neut.v2)\t avg_E_WT(Avg.E.selB_neut.WT)\t avg_E_v1(Avg.E.selB_neut.v1)\t avg_E_v2(Avg.E.selB_neut.v2)\n")
    meanE_vs_t.write("#method.t_tot\t method.t\t E_dyn\n")
    avg_hk_t.write("#method.t_tot\t method.t\t method.avg_hk[k]\n")
    std_hk_t.write("#method.t_tot\t method.t\t std_hk[k]\n")
    avg_h_c_v_vs_t.write("#method.t_tot\t method.t\t h_c1\t h_c2\t h_v1\t h_v2\n")
    avg_Psel.write("#method.t_tot\t method.t\t ni\t method.n_s\t n_a\t P_a\t P_Th\t Psel\n")
    sel_scatter.write("#method.t_tot\t method.t\t i\t method.E_ij\t Eij-method.E0\t P_a\t P_Th\n")

	# unsigned rng_seed = gen_random_seed() # initial value for random number generator
#	 printf("Initializing the state for random # generator with method.seed %lu...\n" % (rng_seed))
#	 init_genrand( rng_seed )
    method.n_round += 1  # initialized as 0
    for method.t in range(method.N_re+1):	# method.N_re: # cycles (~80)
        if method.t>0 and method.n_GC>=method.n_GC_0*math.pow(2,method.N_pr_BB):
            break	# method.n_GC=current BCs, method.N_pr_BB: #prol(~9), breaks when GC population recovers to initial level
        if method.t ==1:
            continue
        if method.t==0: # Activation
            GCR_initiation.GCR_initiation(init_BC, c_v_site, panel_full, nt_seq)
            method.t += 1
        # GCR cycle
        method.t_tot += 1    # initialized as 0
        print "\n\nt=%d, method.t_tot=%d\n\n" % (method.t, method.t_tot)
        c_v_site.write("\n# method.t=%d\n" % (method.t))

        proliferation_SHM.proliferation_SHM(init_BC,output_ns,c_v_site)
        removeLethal.removeLethal(c_v_site)
        selection_mixture_encounter.selection_mixture_encounter(sel_scatter, avg_Psel, meanE_vs_t)
        differentiation_recycle.differentiation_recycle(succ_ts, nt_WT_v1_v2, avg_hc_hv, meanE_vs_t, output_ns)

        if method.flag_end==1:
            E_sel_Hist.close()
            sel_BC.close()
            n_mut.close()
            init_BC.close()
            output_ns.close()
            c_v_site.close()
            nt_br.close()
            mean_E.close()
            br_t.close()
            hc_t.close()
            hv_t.close()
            succ_ts.close()
            nt_WT_v1_v2.close()
            meanE_vs_t.close()
            avg_hk_t.close()
            std_hk_t.close()
            avg_h_c_v_vs_t.close()
            sel_scatter.close()
            avg_Psel.close()

            return 0 # if no GCBs ??? But no manipulation of flag_end always 0

        output_reset.output_reset(n_mut, output_ns, nt_WT_v1_v2, nt_br, br_t, avg_hc_hv, hc_t, hv_t, avg_hk_t, std_hk_t,
                         avg_h_c_v_vs_t, E_sel_Hist, sel_BC,sim_ID)


    method.n_round += 1

    for k in range(1,method.N+1): # method.N # ep sites 18+6+22
        if (k==4 or k==11 or k==16 or k==38 or k==39 or k==40 or k==22 or k==23 or k==24 or k==25 or k==26):
            method.s[1,k]=-1
        else:
            method.s[1,k]=1
    		 # create 11 mutated sites out of 22 var sites(v1)
    for j in range(1,method.n_FDC_Ag+1):
        method.C[j]=0.3*method.f
    	 #different periods (WT, v1, v2), n_FDC for different types of variants

    for method.t in range(1,method.N_re+1): #injecting v1
        method.t_tot += 1
        print ("\n\nt=%d, method.t_tot=%d\n\n" % (method.t, method.t_tot))
        c_v_site.write("\n# method.t=%d\n" % (method.t))
        if (method.t>1 and method.n_GC>=method.n_GC_0*math.pow(2,method.N_pr_BB)): # also when it recovers the initial population
            break

        proliferation_SHM.proliferation_SHM(init_BC,output_ns,c_v_site)
        removeLethal.removeLethal(c_v_site)
        selection_mixture_encounter.selection_mixture_encounter(sel_scatter, avg_Psel, meanE_vs_t)
        differentiation_recycle.differentiation_recycle(succ_ts, nt_WT_v1_v2, avg_hc_hv, meanE_vs_t, output_ns)

        if method.flag_end==1:
        #  goto End # terminates
            E_sel_Hist.close()
            sel_BC.close()
            n_mut.close()
            init_BC.close()
            output_ns.close()
            c_v_site.close()
            nt_br.close()
            mean_E.close()
            br_t.close()
            hc_t.close()
            hv_t.close()
            succ_ts.close()
            nt_WT_v1_v2.close()
            meanE_vs_t.close()
            avg_hk_t.close()
            std_hk_t.close()
            avg_h_c_v_vs_t.close()
            sel_scatter.close()
            avg_Psel.close()
            return 0

        output_reset.output_reset(n_mut, output_ns, nt_WT_v1_v2, nt_br, br_t, avg_hc_hv, hc_t, hv_t, avg_hk_t, std_hk_t,
                     avg_h_c_v_vs_t, E_sel_Hist, sel_BC,sim_ID)

    method.n_round += 1
    for k in range(1,method.N+1):	# V2 sequential injection
        if (k==27 or k==28 or k==29 or k==30 or k==31 or k==32 or k==33 or k==34 or k==35 or k==36 or k==37):
            method.s[1,k]=-1
        else:
            method.s[1,k]=1

    for j in range(1,method.n_FDC_Ag+1):
        method.C[j]=0.3*method.f

    for method.t in range(1,method.N_re+1):
        method.t_tot += 1
        print "\n\nt=%d, method.t_tot=%d\n\n" % (method.t, method.t_tot)
        c_v_site.write("\n# method.t=%d\n" % (method.t))

        if ((method.t>1 and method.n_GC>=method.n_GC_0*math.pow(2,method.N_pr_BB)) or method.t== method.N_re) : #what is about method.N_re?
            method.flag_succ=1
            succ_ts.write("%d\t %d\t %d\n" % (method.t_tot, method.t, method.flag_succ))
            E_sel_Hist.close()
            sel_BC.close()
            n_mut.close()
            init_BC.close()
            output_ns.close()
            c_v_site.close()
            nt_br.close()
            mean_E.close()
            br_t.close()
            hc_t.close()
            hv_t.close()
            succ_ts.close()
            nt_WT_v1_v2.close()
            meanE_vs_t.close()
            avg_hk_t.close()
            std_hk_t.close()
            avg_h_c_v_vs_t.close()
            sel_scatter.close()
            avg_Psel.close()
            return 0

        proliferation_SHM.proliferation_SHM(init_BC,output_ns,c_v_site)
        removeLethal.removeLethal(c_v_site)
        selection_mixture_encounter.selection_mixture_encounter(sel_scatter, avg_Psel, meanE_vs_t)
        differentiation_recycle.differentiation_recycle(succ_ts, nt_WT_v1_v2, avg_hc_hv, meanE_vs_t, output_ns)

        if method.flag_end==1:

            E_sel_Hist.close()
            sel_BC.close()
            n_mut.close()
            init_BC.close()
            output_ns.close()
            c_v_site.close()
            nt_br.close()
            mean_E.close()
            br_t.close()
            hc_t.close()
            hv_t.close()
            succ_ts.close()
            nt_WT_v1_v2.close()
            meanE_vs_t.close()
            avg_hk_t.close()
            std_hk_t.close()
            avg_h_c_v_vs_t.close()
            sel_scatter.close()
            avg_Psel.close()
            return 0

        output_reset.output_reset(n_mut,output_ns,nt_WT_v1_v2,nt_br,br_t,avg_hc_hv,hc_t,hv_t,avg_hk_t,std_hk_t,avg_h_c_v_vs_t,E_sel_Hist,sel_BC,sim_ID)

    E_sel_Hist.close()
    sel_BC.close()
    n_mut.close()
    init_BC.close()
    output_ns.close()
    c_v_site.close()
    nt_br.close()
    mean_E.close()
    br_t.close()
    hc_t.close()
    hv_t.close()
    succ_ts.close()
    nt_WT_v1_v2.close()
    meanE_vs_t.close()
    avg_hk_t.close()
    std_hk_t.close()
    avg_h_c_v_vs_t.close()
    sel_scatter.close()
    avg_Psel.close()
    return 0

if __name__ == "__main__":
    import cProfile
    cProfile.run('main()')
