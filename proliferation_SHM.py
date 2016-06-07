import math
import method
import calc_all_Wi
import differentiation_recycle
import E_sel_hist
import GCR_initiation
import output_reset
import pmut
import removeLethal
import selection_mixture_encounter
#import pdb

def proliferation_SHM(init_BC,output_ns,c_v_site):  #0. start with BC= 3*2^9 1. replicate twice 2. append state matrices 3.for every GC mutate seq i with prob 0.14 per seq per div
#4. select variants of FDC Ags to be seen for BCR (all or none )5. Repeat sampling delX until temp method.h is bounded draw delX from dist
#,Determine Ej_max, if mutation of BCR to any Ag (E_l+dX)> E_jmax then,l is the type for the most strongly binding Ag,if not, randomly choose
#a Ag which to stronlgy bind(Act_Ag_aa) 6. update new h accordingly corresponingly to the chosen epitopes 7.count # mut ben, # mut del,
##novel mut and corresonding sites: ben_mut: affinity improving mut ,del_mut: affinity degrading mut


    n_am=0
    aa_ctr=0
    n_b=0
    pm=0.5*3.*float(method.N/500.)
    ct_h=0
    flag_Ag=0
    E_l=0.
    E_j=0.
    Ej_max=-100.
    diff=0.
    alpha=2.
    flag_k=0 #v site indicator
    flag_mut=0 #mut indicator
    ctr_p1=0
    ctr_n1=0
    ctr_2=0

    print "# of B cells before prol is method.n_GC=%d\n" % (method.n_GC)

    output_ns.write("%d\t %d\t %d\t " % (method.t_tot, method.t, method.n_GC))
     # # cells in diff cpt
    init_BC.write("\n\n# method.t=%d, method.n_GC=%d\n" % (method.t, method.n_GC))

    # replicate BCs
    for m in range(method.N_pr):	    # method.N_pr 2 # repl before sel  # 1. replicate twice 2. append state matrices

        init_BC.write("\nAt the %d prol round\n\n" % (m))

        for i in range(1,method.n_GC+1):
            for k in range(1,method.N+1):
                method.h_GC[i+method.n_GC,k]=method.h_GC[i,k]   # append method.h_GC values for 46 sites of new B

            method.n_mt[i+method.n_GC]=method.n_mt[i]     # # mut
            method.n_mt_b[i+method.n_GC]=method.n_mt_b[i]
            method.n_mt_d[i+method.n_GC]=method.n_mt_d[i]

            method.max_mut[i+method.n_GC]=method.max_mut[i]

            for l in range(1,method.N+1):
                if method.mut_site[i,l]==0:
                    break
                method.mut_site[i+method.n_GC,l]=method.mut_site[i,l]

        if method.t==0:
            init_BC.write("for i=1\n")
            for k in range(1,method.N+1):
                init_BC.write("%.3f\t " % (method.h_GC[1,k]))
            init_BC.write("\n")
            init_BC.write("for i=method.n_GC+1=%d\n" % (method.n_GC+1))

            for k in range(1,method.N+1):
                init_BC.write("%.3f\t " % (method.h_GC[method.n_GC+1,k]))
            init_BC.write("\n")

        method.n_GC*=2

        # set method.N_A_aa[i]    # # of Ags that BCR can see in a round of selection& method.Ag_label_aa: label of them 0,1,2
        for i in range(1,method.n_GC+1): # initialize the array : # mut =1 , type mut =0
            method.N_A_aa[i]=1 # in this case, only sees 1
            for l in range(1,method.n_FDC_Ag+1):# for # of Ag variants in GC
                method.Ag_label_aa[i,l]=0      # in this case only recognizes WT?? what does this mean?
            method.flag_aa[i]=0

        for i in range(1,method.n_GC+1): # for every GC
            rd1= method.rand_gen()
            if rd1<=pm: #mutate seq i with prob 0.14 per seq per div
                rd2= method.rand_gen() # probabilit of mutation
                if rd2>=(1.-method.pL):
                    method.flag[i]=1	#	method.pL 0.3  probab of lethal mut If a lethal mut
                    continue

                if rd2 <=method.pA: #	aa-mut method.pA 0.2  probab of aa mut

                    c_v_site.write("GC B cells with Affinity affecting mutation: method.n_FDC_Ag=%d, # of FDC Ag seen(method.N_A_aa[%d])=%d\n" % (method.n_FDC_Ag, i, method.N_A_aa[i]))
                    c_v_site.write( "Ag labels are\n")

                    # select  FDC Ags to be seen

                    l=1
                    condition = True
                    while condition:
                        #Randomly assign the label of Ag variants that B cell sees/ round of GC

                        rd=1+int(method.rand_gen()*method.n_FDC_Ag)
                         # Returns a uniform random // (int)()to the biggest integer * select type of Ag to be seen
                        if l>1:
                            for n in range(1,l):  # if finds rd==method.Ag_label_aa[i,n] has assigned previously, then generate new rd ( choose different aa to mutate every single time)
                                if rd==method.Ag_label_aa[i,n]:
                                    flag_Ag=1
                                    break
                            if flag_Ag==1:
                                flag_Ag=0
                                condition = l <= method.N_A_aa[i]
                                continue
                                 # randomly choose label of Ag variant, if 1>1 then repeat assign rd, if it hasn'method.t already been assigned,repeat rd assignment if the same rd is chosen

                        method.Ag_label_aa[i,l]=rd# determine the label of Ag
                        # int method.N_A_aa[method.n_B+1]=0 # & label of FDC Ag seen by aa-mut BCs
                        # int method.Ag_label_aa[method.n_B+1,6]=0 stores Ag type selection for every GC B

                        c_v_site.write("# of Ag seen:l=%d, Label of Ag seen: method.Ag_label_aa[%d,%d]=%d\n" % (l, i, l, method.Ag_label_aa[i,l]))
                        l += 1
                        condition = l<=method.N_A_aa[i]

                    condition = True
                    while condition:
                    # Repeat sampling delX until temp method.h is bounded 1. draw delX from dist 2. Determine Ej_max
    #                3. if mutation of BCR to any Ag (E_l+dX)> E_jmax then  l is the type for the most strongly binding Ag(in reality
    #                this is true as long as there is one dominant Ag type)
    #                4. if not, randomly choose a Ag which to stronlgy bind 5. update method.h= method.h+ delh
                        flag_h=1
                        method.flag_dX=0
                        p= method.rand_gen()
                        dX=pmut.Pmut(p) #get delX drawn from distribution

                        #find ref seq for mut for method.N_A_aa =1, Act_ag= method.Ag_label_aa
#                        pdb.set_trace()
                        for l in range(1,int(method.N_A_aa[i])+1):#What's the point of method.Ej max?
                            E_l=0.
                            for k in range(1,method.N+1):
                                E_l+=method.h_GC[i,k]*method.s[method.Ag_label_aa[i,l]][k]
                                #GC B cell method.h values*FDC_Ag seq identity
                            for n in range(1,method.N_A_aa[i]+1): #find max(Eij'):1.Calculate E for all Ags that BCR can see 2.Set the largest as the Emax)=ini.-100
                                if n==l:
                                    condition = flag_h == 0
                                    continue   # goes to n != l

                                for k in range(1,method.N+1):
                                    E_j+=method.h_GC[i,k]*method.s[method.Ag_label_aa[i,n]][k]

                                if E_j>Ej_max:
                                    Ej_max=E_j

                                E_j=0.

                            if (E_l+dX>Ej_max) and (E_l+dX-Ej_max>diff):
                                #found j,if mutated El+delX > Ej_max, set it as the most strongly binding Ag&& the one that has the most difference,init.diff=0*/

                                for k in range(1,method.N+1):
                                    method.act_Ag[i,k]=method.s[method.Ag_label_aa[i,l]][k]

                                method.flag_aa[i]=1
                                diff=E_l+dX-Ej_max

                            Ej_max=-100.
                        diff=0.

                        if method.flag_aa[i]==0: # if all interaction strenght are same or has negative delX ,then randomly choose one of them
                            rd=1+int(method.rand_gen()*method.N_A_aa[i])

                            for k in range(1,method.N+1):
                                 method.act_Ag[i,k]=method.s[method.Ag_label_aa[i,rd]][k]

                            method.flag_aa[i]=1

                        #1) make the mutation 2) calcualte method.h= method.h + del method.h = method.h + dX/method.s(I don'method.t understand this part) */
                        k=int(method.N*method.rand_gen())+1
                        if k==(method.N+1):
                            k=method.N
                        temp_h=method.h_GC[i,k]+dX/float(method.act_Ag[i,k])
                        #bounding
                        if temp_h<-1.5 or temp_h>1.5:
                            flag_h=0
                        condition = flag_h==0


                    #v site or not
#                    pdb.set_trace()
                    flag_k=0

                    for l in range(method.n_v_site): #k is the chosen site of mutation, if it's in Var set then method.flag==1 and break
                        k_v=method.v_site_label[l]
                        if k==k_v and k<=40:
                            flag_k=1
                            break

                    if method.act_Ag[i,k]==1: # if strain has WT on that residue flag_k=0, mutated var site flag_k= 1
                        flag_k=0

                    #MUTATED v site w |method.h'_k|<|h_k|

                    if flag_k==1 and math.fabs(temp_h)<math.fabs(method.h_GC[i,k]):
                         #bounded between [-1.5,inf), aff enhancement
                        k_c=method.c_site_label[int(method.rand_gen())*method.n_c_site]
                        dh_c=alpha*(math.fabs(method.h_GC[i,k])-math.fabs(temp_h))
                        method.h_GC[i,k_c]+=dh_c
                        method.flag_dX=1
                        ctr_p1 += 1

                    #MUTATED v site w |method.h'_k|>|h_k|

                    if flag_k==1 and math.fabs(temp_h)>math.fabs(method.h_GC[i,k]): #bounded between [-1.5,inf),aff deduction
                        k_c=method.c_site_label[int(method.rand_gen())*method.n_c_site]

                        dh_c=alpha*(math.fabs(method.h_GC[i,k])-math.fabs(temp_h))

                        if dh_c<(-1.5-method.h_GC[i,k_c]):
                            dh_c=-1.5-method.h_GC[i,k_c]

                        method.h_GC[i,k_c]+=dh_c
                        method.flag_dX=-1
                        ctr_n1 += 1


                    #type II c sites

                    if k>40: #mutated or unmutated V sites(structure motifs)- shielding mutation, bounded by [-1.5,1.5]

                        k_c=method.c_site_label[int(method.rand_gen())*method.n_c_site]
                        dh_c=alpha*(math.fabs(method.h_GC[i,k])-math.fabs(temp_h))

                        if dh_c>(1.5-method.h_GC[i,k_c]):
                            dh_c=1.5-method.h_GC[i,k_c]
                        if dh_c<(-1.5-method.h_GC[i,k_c]):
                            dh_c=-1.5-method.h_GC[i,k_c]

                        method.h_GC[i,k_c]+=dh_c
                        method.flag_dX=2
                        ctr_2 += 1


                    #paratope interacting method.W/ conserved sites or method.W/ unmutated V sites
                    method.h_GC[i,k]=temp_h
#                    pdb.set_trace()
                    n_am += 1
                    method.n_mt[i]+= 1

                    #ben(affinity improving) & del(affinity degrading) mut
                    if dX>1e-6:
                        method.n_mt_b[i]+= 1
                    if dX<-1e-6:
                        method.n_mt_d[i]+= 1

                    if method.max_mut[i]>0: #method.max_mut stores for all GC Bs,method.mut_site: mutated Paratope locus
                        for n in range(1,method.max_mut[i]+1): #check whether it has already counted mut,break if it is not novel
                            if k==method.mut_site[i,n]:
                                flag_mut=1
                                break
                    if flag_mut==0: #novel mutated site store # mut,mutated paratop */
                        method.max_mut[i]+= 1
                        method.mut_site[i,method.max_mut[i]]=k

                    else:
                        flag_mut=0

                    aa_ctr += 1
                    if dX>1e-6:
                        n_b += 1

                    if (i==method.n_GC) and (k==method.N) and (m==method.N_pr):
                        aa_ctr=0
                        n_b=0

    print "At method.t=%d, # BCs after %d rounds of prol is method.n_GC=%d\n" % (method.t, method.N_pr, method.n_GC)
    print "# aa mut per BC is n_am/method.n_GC=%d/%d=%.3f\n" % (n_am, method.n_GC, float(n_am)/float(method.n_GC))

    init_BC.write("# BC after %d rounds of prol is method.n_GC=%d\n" % (method.N_pr, method.n_GC))

    init_BC.write("# aa mut per BC is n_am/method.n_GC=%d/%d=%.3f\n" % (n_am, method.n_GC, float(n_am)/float(method.n_GC)))
