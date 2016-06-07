import method
import proliferation_SHM
import calc_all_Wi
import differentiation_recycle
import GCR_initiation
import output_reset
import pmut
import removeLethal
import selection_mixture_encounter

def E_sel_hist(j,E_sel_Hist): # sort BCs into 40 bins of Energy histogram, store it in E_sel_Hist: 1. Calc E for Ag type j,2.sort it into bins,
# 3. store it in E_sel_Hist


    E_max=-10.
    E_min=1000.
    mut_ctr=0
    Ent=1.05*method.E0
    avg_n_nt=0.
    avg_br=0.
    var_n_nt=0.
    var_br=0.
    Eij=0.

    for i in range(1, method.n_sel+1): # for selected BCs initialize
        method.Ej[i]=0.

    for l in range(method.Bin_E_sel): #for #of bins:40 initialize method.n, which counts how many falls into such bin
        method.n[l]=0

    for i in range(1,method.n_sel+1):
        for k in range(1,method.N+1):
            method.Ej[i]+=method.h_sel[i,k]*method.s[j,k]

    E_min=0.
    E_max=45.

    delta_E=(E_max-E_min)/float(method.Bin_E_sel)

    for i in range(1, method.n_sel+1): #bin each BC in one of 40 bins, if method.Ej falls into the largest bin,then fprint such BC,method.Ej, #of such
    # BCS in the bin, # of muts in the BC

        for l in range(method.Bin_E_sel):
            E_l=E_min+l*delta_E

            if l == (method.Bin_E_sel-1):
                if (method.Ej[i] >= E_l) and (method.Ej[i] <=E_max):
                    E_sel_Hist.write("l=%d\n" % (l))
                    method.n[l]+=1

                    mut_ctr += method.n_mt_sel[i]  ## mut in selected BC
                    E_sel_Hist.write("method.Ej[%d]=%.3f, method.n[%d]=%d, method.n_mt_sel[%d]=%d\n" % (i, method.Ej[i], l, method.n[l], i, mut_ctr))
                    mut_ctr =0

            else:
                if (method.Ej[i] >= E_l) and (method.Ej[i] < (E_l+delta_E)):
                    method.n[l]+=1

    E_sel_Hist.write( "# l\t E_l\t method.n[l]\t method.P[l]\n")

    for l in range(method.Bin_E_sel):  #fprint histogram of all BCs vs energy
        method.P[l] = float(method.n[l])/float(method.n_sel)
        E_sel_Hist.write("%d\t %.3f\t %d\t %.3f\n" % (l, E_min+(l+0.5)*delta_E, method.n[l], method.P[l]))

    E_sel_Hist.write(" method.E_sel_min=%.3f, method.E_sel_max=%.3f, del_E=%.3f\n" % (E_min, E_max, delta_E))
    E_sel_Hist.write("Total # of sel BCs is %d\n\n" % (method.n_sel))
