import method
import proliferation_SHM
import calc_all_Wi
import differentiation_recycle
import E_sel_hist
import GCR_initiation
import output_reset
import pmut
import selection_mixture_encounter

def removeLethal(c_v_site): #for mutation survived. BCs , take method.h_GC,n #mt,#mt_ben,#mt_del,#distinct mut 2.pass method.N_A_aa & method.Ag_label_aa for aa-mutated BCs
 # generate method.N_A & method.Ag_label for un-mutated BCs /*randomly choose label of Ag variant, if 1>1 then repeat assign rd, if it hasn'.t already been assigned,repeat rd assignment if the same rd is chosen */

    #
    n_L=0
    flag_Ag=0
    method.n_s=0

    for i in range(1,method.n_GC+1): # cnt # of lethal mutations
        if method.flag[i] ==1:
            n_L += 1
            method.flag[i]=0
            continue

        else: #for survived. BCs(from mutation) , take method.h_GC,n #mt,#mt_ben,#mt_del,#distinct mut

            method.n_s += 1

            for k in range(1,method.N+1): #method.h_s: surv BCs
                method.h_s[method.n_s,k]=method.h_GC[i,k]

            method.n_mt_s[method.n_s]=method.n_mt[i]
            method.n_mt_b_s[method.n_s]=method.n_mt_b[i] #
            method.n_mt_d_s[method.n_s]=method.n_mt_d[i] #
            method.max_mut_s[method.n_s]=method.max_mut[i]

            for l in range(1,method.N+1):
                if method.mut_site[i,l]==0:
                    break #if no mut
                method.mut_site_s[method.n_s,l]=method.mut_site[i,l]
                #transfer array

         # pass method.N_A_aa & method.Ag_label_aa for aa-mutated BCs
            if method.flag_aa[i]==1: # for aa-mutated BCs
                method.N_A[method.n_s]=method.N_A_aa[i] #how many Ag BCs see

                for l in range(1, method.N_A_aa[i]+1):
                     method.Ag_label[method.n_s,l]=method.Ag_label_aa[i,l] #transfer label of Ag seen
                     c_v_site.write("l=%d, method.Ag_label[%d,%d]=%d\n" % (l, method.n_s, l, method.Ag_label[method.n_s,l]))

         # generate method.N_A & method.Ag_label for un-mutated BCs

            if method.flag_aa[i]==0:
                method.N_A[method.n_s]=1
                l=1
                condition= True
                while condition:
                #understand break and continue clearly
                    rd=1+int(method.rand_gen()*method.n_FDC_Ag) #randomly choose label of Ag variant, if 1>1 then repeat assign rd, if it hasn'method.t already been assigned,repeat rd assignment if the same rd is chosen
                    if l>1:
                        for n in range(1,l):
                            if rd==method.Ag_label[method.n_s,n]:
                                flag_Ag=1
                                break
                        if flag_Ag==1:
                            flag_Ag=0
                            condition = l <= method.N_A[method.n_s]
                            continue

                    method.Ag_label[method.n_s,l]=rd
                    c_v_site.write("l=%d, method.Ag_label[%d,%d]=%d\n" % (l, method.n_s, l, method.Ag_label[method.n_s,l]))

                    l += 1
                    condition = l<=method.N_A[method.n_s]

    print "frac of BCs bearing lethal mut is %d/%d=%.3f\n" % (n_L, method.n_GC, float(n_L)/float(method.n_GC))
    print "frac of surv BCs (from mutation) is %d/%d=%.3f\n" % (method.n_s, method.n_GC, float(method.n_s)/float(method.n_GC))
