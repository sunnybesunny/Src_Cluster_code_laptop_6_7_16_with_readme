import math
import method
import proliferation_SHM
import calc_all_Wi
import differentiation_recycle
import E_sel_hist
import GCR_initiation
import output_reset
import pmut
import removeLethal
import pdb

def selection_mixture_encounter(sel_scatter,avg_Psel,meanE_vs_t):
#for surv. BCs calc W based on act.thresh(10.8) 2.Prob.Ag internal.
#generate rd -> if selec. by Ag-> if comp>0, then compete for Tfh (calcualte pth) 3.calculate dynamic threshol. by sum all
#survived BCs assuming they only see Ag typ 1 + weight of memory cells ...

    sum1=0.
    avg_P_a=0.
    avg_P_Th=0.
    sum_W=0. #total weight
    n_a=0
    Eij=0.
    Ei=0.
    sum_pr=0.
    n_i=0
    ni=0
    method.n_sel=0

    ni=calc_all_Wi.calc_all_Wi()

    for i in range(1,method.n_s+1): #for surv. BCs
        P_a=0.
        P_Th=0.
        # activation
        sum1=0.

        if method.N_A[1]==1: #calc method.W based on act.thresh(10.8)
            Eij=0.
            for k in range(1,method.N+1):
                Eij+=method.h_s[i,k]*method.s[method.Ag_label[i,1],k]
                 #potential improvement
            sum1=method.C[method.Ag_label[i,1]]*math.exp((Eij-method.E0)/0.6)

        if method.n_FDC_Ag>1 and method.N_A[1]==method.n_FDC_Ag:
            for j in range(1,method.n_FDC_Ag+1):
                Eij=0.
                for k in range(1,method.N+1):
                    Eij+=method.h_s[i,k]*method.s[j,k]
                sum1+=method.C[j]/float(method.n_FDC_Ag)*math.exp((Eij-method.E0)/0.6)

        P_a=sum1/(1.+sum1) #Prob.Ag internal.
        avg_P_a+=P_a

        # competition
        # pdb.set_trace()
        if method.rand_gen()<=P_a: #if selected by Ag
            n_a += 1
#            pdb.set_trace()
            if ni==1:
                P_Th=1.
                avg_P_Th+=P_Th
                if method.rand_gen()<=P_Th:
                    method.n_sel += 1
                    for k in range(1,method.N+1): #
                        method.h_sel[method.n_sel,k]=method.h_s[i,k]
                    for j in range(1,method.N_A[i]+1):
                        method.Ag_label_sel[method.n_sel,j]=method.Ag_label[i,j]
                    method.n_mt_sel[method.n_sel]=method.n_mt_s[i]
                    method.n_mt_b_sel[method.n_sel]=method.n_mt_b_s[i]
                    method.n_mt_d_sel[method.n_sel]=method.n_mt_d_s[i]
                    method.max_mut_sel[method.n_sel]=method.max_mut_s[i]
                    for l in range(1,method.N+1):
                        if method.mut_site_s[i,l]==0:
                            break
                        method.mut_site_sel[method.n_sel,l]=method.mut_site_s[i,l]
                 #only 1 BC

            elif ni>1: #if competitors>0
                sum_W=0.
                for l in range(1,ni+1):
                    if l==i:
                        continue
                    sum_W+=method.W[l] #sum all weights except own
                ratio2=method.C[1]*method.W[i]/(sum_W/float(ni-1))
                P_Th=ratio2/(1.+ratio2) #prob of getting Tfh
                avg_P_Th+=P_Th


        else:
            continue # lose competition

        if method.rand_gen() <= P_Th:
            method.n_sel += 1
            for k in range(1, method.N + 1):  #
                method.h_sel[method.n_sel, k] = method.h_s[i, k]
            for j in range(1, method.N_A[i] + 1):
                method.Ag_label_sel[method.n_sel, j] = method.Ag_label[i, j]
            method.n_mt_sel[method.n_sel] = method.n_mt_s[i]
            method.n_mt_b_sel[method.n_sel] = method.n_mt_b_s[i]
            method.n_mt_d_sel[method.n_sel] = method.n_mt_d_s[i]
            method.max_mut_sel[method.n_sel] = method.max_mut_s[i]
            for l in range(1, method.N + 1):
                if method.mut_site_s[i, l] == 0:
                    break
                method.mut_site_sel[method.n_sel, l] = method.mut_site_s[i, l]
                # only 1 BC

        if method.N_A[1]==1: # if see only 1 ag
            sel_scatter.write("%d\t %d\t %d\t %.2f\t %.3f\t %.3f\t %.3f\n" % (method.t_tot, method.t, i, Eij, Eij-method.E0, P_a, P_Th))
        #i

    sel_scatter.write("\n")

    if method.n_s>0: # surv BC >0
        avg_P_a/=float(method.n_s)
    else:
        avg_P_a=0.

    if n_a>0: # sel BC >0
        avg_P_Th/=float(n_a)
    else:
        avg_P_Th=0.

    avg_Psel.write( "%d\t %d\t %d\t %d\t %d\t %.2f\t %.2f\t %.3f\n" % (method.t_tot, method.t, ni, method.n_s, n_a, avg_P_a, avg_P_Th, avg_P_a*avg_P_Th))

    for i in range(1,method.n_s+1): #for all surv. BCs

        n_i += 1 #counts until n_i = method.n_s
        for k in range(1,method.N+1): #sum E for 1st Ag seen by BC
            Ei+=method.h_s[i,k]*method.s[method.Ag_label[i,1],k]
        sum_pr+=math.exp(Ei/0.6) # sum of surv BCs affin. to Ag type #1
        Ei=0.

    if method.all_MC_flag==1:
        for i in range(1,method.n_m+1):
            n_i += 1
            for k in range(1,method.N+1): #for all memory cells
                Ei+=method.h_m[i,k]*method.s[method.actAg_MC[i],k] #calc E values and sum up all the weights
            if i<=method.n_m and i>(method.n_m-10): # method.n_m-10<i<method.n_m
                print "method.actAg_MC[%d]=%d\n" % (i, method.actAg_MC[i])
            sum_pr+=math.exp(Ei/0.6)
            Ei=0.

    E_dyn=0.6*math.log(sum_pr/float(n_i)) #dynamic threshold
    meanE_vs_t.write("\n%d\t %d\t %.2f\t " % (method.t_tot, method.t, E_dyn))