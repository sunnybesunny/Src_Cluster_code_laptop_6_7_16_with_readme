import math
import method
import proliferation_SHM
import differentiation_recycle
import E_sel_hist
import GCR_initiation
import output_reset
import pmut
import removeLethal
import selection_mixture_encounter

def calc_all_Wi(): #calculate weight based on method.Ec, returns # of sel.BCs

    n_i=0
    Eij=0.

    for i in range(1,method.n_s+1):  # for surv,BCs
        n_i += 1
        method.W[n_i]=0.         # method.W:weight of competing cells

        if method.N_A[1]==1: #if BCs see 1, (all BCs see the same number of Ags)
            Eij=0.
            for k in range(1, method.N+1):
                Eij+=method.h_s[i,k]*method.s[method.Ag_label[i,1],k]
            method.W[n_i]=method.C[method.Ag_label[i,1]]*math.exp((Eij-method.Ec)/0.6) #Calc weight

        if method.n_FDC_Ag>1 and method.N_A[1]==method.n_FDC_Ag: # if BCs see more than 1 Ags
            for j in range(1,method.n_FDC_Ag+1):
                Eij=0.
                for k in range(1,method.N+1):
                    Eij+=method.h_s[i,k]*method.s[j,k]
                method.W[n_i]+=method.C[j]/float(method.n_FDC_Ag)*math.exp((Eij-method.Ec)/0.6) #since Ctot remains the same

    if method.all_MC_flag==1: #memory cell flag, if 1, include MCs into competition

        for i in range(1,method.n_m+1):	# method.n_m is 0
            n_i += 1
            method.W[n_i]=0.

            if method.N_A[1]==1:
                method.actAg_MC[i]=1+int(method.rand_gen()*method.n_FDC_Ag)
                Eij=0.
                for k in range(1,method.N+1):
                    Eij+=method.h_m[i,k]*method.s[method.actAg_MC[i],k]
                method.W[n_i]=method.C[method.actAg_MC[i]]*math.exp((Eij-method.Ec)/0.6)

            if method.n_FDC_Ag>1 and method.N_A[1]==method.n_FDC_Ag:
                for j in range(1,method.n_FDC_Ag+1):
                    Eij=0.
                    for k in range(1,method.N+1):
                        Eij+=method.h_m[i,k]*method.s[j,k]
                    method.W[n_i]+=method.C[j]/float(method.n_FDC_Ag*math.exp((Eij-method.Ec)/0.6))
    return n_i
