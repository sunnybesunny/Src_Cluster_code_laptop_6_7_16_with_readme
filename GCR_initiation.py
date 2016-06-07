import numpy as np
import method
import proliferation_SHM
import calc_all_Wi
import differentiation_recycle
import E_sel_hist
import output_reset
import pmut
import removeLethal
import selection_mixture_encounter
#import pdb

def GCR_initiation(init_BC, c_v_site,panel_full, nt_seq):
  # 1. creates 3 activated blasts with random site energies(Ei>=(method.E0-0.5*method.R)&&Ei<(method.E0+0.5*method.R))
  # 2. proliferate 9 times
# 3. creates panel sequence array for 46 sites/ total sites
  # 4. initialize Act_ag

    Ei=0.
    flag_in=0
    n_v=0

     # c & v labels Create c site label complementary to v site label and store the site number in c_v_site
    for k in range(1,method.N+1): 	# method.N # sites 46
        for l in range(0,method.n_v_site):  # method.n_v_site 28 : not conserved region
            if k==method.v_site_label[l]:
                flag_in=1
                break      # variable site label
        if flag_in==1:
            flag_in=0
            continue
        c_v_site.write("k=%d, l=%d\n" % (k, (l+1)))
        method.c_site_label[method.n_c_site]=k
        method.n_c_site += 1 # method.n_c_site 0
    for l in range(method.n_c_site):
        c_v_site.write("%d\t " % (method.c_site_label[l]))
    c_v_site.write("\n\n")
    for l in range(method.n_v_site):
        c_v_site.write("%d\t " % (method.v_site_label[l]))
    c_v_site.write("\n\n")

    # FDC Ag
    for k in range(1,method.N+1):

        if k==1 or k==2 or k==3 or k==5 or k==6 or k==7 or k==8 or k==9 or\
            k==10 or k==12 or k==13 or k==14 or k==15 or k==17 or k==18 or k==19\
            or k==20 or k==21:
            method.s[1,k]=1
        else:
            method.s[1,k]=1

    method.n_FDC_Ag += 1

    for j in range(1,method.n_FDC_Ag+1):
        method.C[j]=0.3*method.f

    # Create WT, V1, V2 Immunogens

    for k in range(1,method.N+1):

        if k==1 or k==2 or k==3 or k==5 or k==6 or k==7 or k==8 or k==9 or\
            k==10 or k==12 or k==13 or k==14 or k==15 or k==17 or k==18 or k==19\
            or k==20 or k==21:
            method.s_im_0[k]=1
        else:
            method.s_im_0[k]=1

    for k in range(1,method.N+1):

        if k==4 or k==11 or k==16 or k==38 or k==39 or k==40 or k==22 or k==23\
            or k==24 or k==25 or k==26:
            method.s_im_1[k]=-1
        else:
            method.s_im_1[k]=1

    for k in range(1,method.N+1):

        if k==27 or k==28 or k==29 or k==30 or k==31 or k==32 or k==33 or k==34\
            or k==35 or k==36 or k==37:
            method.s_im_2[k]=-1
        else:
            method.s_im_2[k]=1

    # B blasts for 10000 B cells, generate site energy randomly and only those exceed threshol survive
    for i in range(1,method.n_B+1):

        if method.n_GC==3:
            break
        for k in range(1,method.N+1): # Generate site energy randomly
            method.h[k]=method.hl+(method.hu-method.hl)*method.rand_gen() #method.hl,method.hu: range of method.h for intiail BCR sequence
            Ei+=method.h[k]

        if Ei>=(method.E0-0.5*method.R) and Ei<(method.E0+0.5*method.R): # Check activation threshold, generate GC for each surv B blast until generating 3

            method.n_GC += 1
            for k in range(1,method.N+1):
                method.h_GC[method.n_GC,k]=method.h[k]
                method.h_GC_0[method.n_GC,k]=method.h[k]

            init_BC.write("Activated B blast #%d has E[%d]=%.3f\n" % (method.n_GC, method.n_GC, Ei))

        Ei=0.
#    pdb.set_trace()
    print "Generated %d B blasts\n" % (method.n_GC) # why do you need this method.n_GC=3 always?

    method.n_GC_0=method.n_GC
    init_BC.write( "\n# of B blasts is %d\n" % (method.n_GC_0))

    # output seq of B blasts: store method.h values for 3 activated B blasts
    for i in range(1,method.n_GC_0+1):

        init_BC.write("\nB blast #%d\n" % (i))

        init_BC.write("v sites: ")

        for l in range(method.n_v_site):
            k_v=method.v_site_label[l]
            init_BC.write("%.3f\t " % (method.h_GC_0[i,k_v]))
        init_BC.write("\n")

        init_BC.write("c sites: ")
        for l in range(method.n_c_site):

            k_c=method.c_site_label[l]
            init_BC.write("%.3f\t " % (method.h_GC_0[i,k_c]))
        init_BC.write( "\n")

    for m in range(method.N_pr_BB): # B blasts prolif 9 cycles + append the status matri.
        for i in range(1,method.n_GC+1):

            for k in range(1,method.N+1):
                method.h_GC[i+method.n_GC,k]=method.h_GC[i,k]
        method.n_GC*=2

    init_BC.write( "\n# of B blasts that initiate GCR is %d\n" % (method.n_GC))

    for i in range(1,method.n_GC+1):	#initializing method.act_Ag n_GCbyN matrix method.act_Ag: Ag that B cells most strongly binds to
        for k in range(1,method.N+1):
            method.act_Ag[i,k]=1

    # panel seq & binding segment open up full panel sequence and store method.s value that corresponds to 46 sites/ total sequence for each panel sequence
    method.s_np = np.loadtxt(panel_full)
    off_set = np.zeros((1,method.s_np.shape[1]))
    method.s_np=np.concatenate((off_set,method.s_np),axis=0)
    off_set = np.zeros((method.s_np.shape[0], 1))
    method.s_np = np.concatenate((off_set, method.s_np), axis=1)

    for j in range(1,method.n_Ag+1): # method.n_Ag #panel sequence method.s_np full panel sequence
        # for k in range(1,method.Ntot+1): # #total lengths
        #
        #     fscanf(panel_full, "%d" % (&method.s_np[j,k]))
        for l in range(method.N):
            k_np=method.label_np[l] #method.label_np: actual location of epitopes
            n_v += 1
            method.s0[j, n_v] = method.s_np[j, k_np]
        # print method.s_np
        # print method.s0
#            try: method.s0[j,n_v]=method.s_np[j,k_np]
#            except  IndexError:
#                print method.s0.shape, method.s_np.shape,j ,n_v,k_np
#                raise
        for k in range(1,n_v+1):
            nt_seq.write( "%d " % (method.s0[j,k]))
        nt_seq.write("\n")
        n_v=0
#    pdb.set_trace()
    panel_full.close()
    nt_seq.close()
