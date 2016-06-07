import math
import method
import proliferation_SHM
import calc_all_Wi
import differentiation_recycle
import E_sel_hist
import GCR_initiation
import pmut
import removeLethal
import selection_mixture_encounter

def output_reset(n_mut,output_ns,nt_WT_v1_v2,nt_br,br_t,avg_hc_hv,hc_t,hv_t,avg_hk_t,std_hk_t,avg_h_c_v_vs_t,E_sel_Hist,sel_BC,sim_ID): #1. calc clone total b_mut,d_mut,#_mut 2.Average novel mutations in a clone 3. affinity: calc affinity for each Ag
#->Comp method.W/ threshold ->calc Avg of affinity beyond threshold 4.breadth calc: for all selected BCs->calc affin b/w each of pan.seq.->
# if Eij> act.thresh * nt.factor # nt_Ag ++ 5.for each residue in conserved region,cal num avg h_k & var_h_k over sel BCs
# 6.# avg & stdev of method.h stored in the order of tp 1 c -> v1 -> v2 -> tp 2 method.C 7.site avg over the same category of epitopes
# 7.E hist of sel BCs 8. h_c, h_v averaged over the same type of sites + averaged over survived BCs 9.  h_c, h_v1,h_v2 averaged over
# the same type of sites + averaged over selected BCs 10.variance of method.hv,method.hv1,method.hv2

    f_name =  "stat_num_sim_%d.txt" %(sim_ID)
    for_stat = open(f_name, "w")
    for_stat.write("#t_tot\t t\t avg_br\t avg_hv\t avg_hc\n")

    avg_hv=0.
    avg_hc=0.
    var_hv=0.
    var_hc=0.   # clone-to-clone var
    avg_hv1=0.
    avg_hv2=0.
    var_hv1=0.
    var_hv2=0.
    avg_hv_s=0.
    avg_hc_s=0.
    Et=0.  # total affinity of all sel BCs to an 'avg' FDC-Ag
    ct=0
    Eij=0.

    avg_br=0.
    var_br=0. # mean & var of breadth
    avg_n_ntAg=0.
    var_n_ntAg=0.
    max_n_ntAg=0
    min_n_ntAg=1000 # max & min # neut Ags for sel BCs
    n_mut_tot=0
    ct_b=0
    ct_d=0
    n_mut_b_tot=0
    n_mut_d_tot=0 # total # ben & del mut
    i_max_br=0
    i_min_br=0
    avg_h_k=0.
    var_h_k=0.
    E_nt=11.34 #neut. threshold
    n_nt_WT=0
    n_nt_v1=0
    n_nt_v2=0 ## rc BCs that can neut WT, v1, v2
    E_WT=0.
    E_v1=0.
    E_v2=0.    #affinity of rc BCs
    avg_E_WT=0.
    avg_E_v1=0.
    avg_E_v2=0.
    max_n_mt=1
    min_n_mt=method.N
    avg_max_mut=0. #avg # distinct mut
    avg_h_c1=0.
    avg_h_c2=0.
    avg_h_v1=0.
    avg_h_v2=0.

    print "At method.t=%d, method.n_s=%d, method.n_sel=%d, method.n_m=%d\n" % (method.t, method.n_s, method.n_sel, method.n_m)

    n_mut.write("\n# method.t=%d\n" % (method.t))
    n_mut.write("\n# method.n_sel(after Tfh selection)=%d\n" % (method.n_sel))

    sel_BC.write("\n# method.n_s=%d, method.n_sel=%d, method.n_m=%d\n" %(method.n_s, method.n_sel, method.n_m))

    # total # mut
    for i in range(1,method.n_sel+1):  #calc clone total b_mut,d_mut,#_mut
        ct_b += method.n_mt_b_sel[i]	#clone total beneficial mutation
        ct_d += method.n_mt_d_sel[i]	#clone total deleterioious mutation
        ct += method.n_mt_sel[i]
        n_mut_tot += ct
        n_mut_b_tot += ct_b
        n_mut_d_tot += ct_d
        ct=0
        ct_b=0
        ct_d=0

    # distinct mut

    for i in range(1,method.n_sel+1): #Average novel mutations in a clone,max_n_mt= max # of distinct mut, min_n_mt= min # of distinct. mut
        avg_max_mut+=float(method.max_mut_sel[i])
        if method.max_mut_sel[i]>max_n_mt:
            max_n_mt=method.max_mut_sel[i]
        if method.max_mut_sel[i]<min_n_mt:
            min_n_mt=method.max_mut_sel[i]
#i
    avg_max_mut/=float(method.n_sel)

    n_mut.write("Total # of mut on %d sel BCs is %d, # mut per sel BC is %.3f\n" % (method.n_sel, n_mut_tot, float(n_mut_tot)/float(method.n_sel)))
    n_mut.write("\nTotal # of ben mut is n_b=%d, total # of del mut is n_d=%d\n" % (n_mut_b_tot, n_mut_d_tot))

    n_mut.write("# of ben mut per sel B cell is m_b=%.3f, # of del mut per sel B cell is m_d=%.3f\n" % (float(n_mut_b_tot)/float(method.n_sel), float(n_mut_d_tot)/float(method.n_sel)))

    output_ns.write("%.3f\t %.3f\t %.3f\t %.2f\t %d\t %d\t %d\n" % (float(n_mut_tot)/float(method.n_sel),float(n_mut_b_tot)/float(method.n_sel), float(n_mut_d_tot)/float(method.n_sel),
                      avg_max_mut, max_n_mt, min_n_mt, method.n_m))

    method.n_mut_tot=0
    method.n_mut_b_tot=0
    method.n_mut_d_tot=0

    # affinity

    if method.n_sel>0:
        for i in range(1,method.n_sel+1):
            for k in range(1,method.N+1): #1.calc affinity for each Ag 2. Comp method.W/ threshold 3.calc Avg of affinity beyond threshold
                E_WT+=method.h_sel[i,k]*method.s_im_0[k]
                E_v1+=method.h_sel[i,k]*method.s_im_1[k]
                E_v2+=method.h_sel[i,k]*method.s_im_2[k]

            if E_WT >= E_nt:
                n_nt_WT += 1
                avg_E_WT += E_WT
            if E_v1 >= E_nt:
                n_nt_v1 += 1
                avg_E_v1 += E_v1
            if E_v2 >= E_nt:
                n_nt_v2 += 1
                avg_E_v2 += E_v2

            E_WT=0.
            E_v1=0.
            E_v2=0.

        if (n_nt_WT==0):
            avg_E_WT=0.
        else:
            avg_E_WT/=float(n_nt_WT)
        if (n_nt_v1==0):
            avg_E_v1=0.
        else:
            avg_E_v1/=float(n_nt_v1)
        if (n_nt_v2==0):
            avg_E_v2=0.
        else:
            avg_E_v2/=float(n_nt_v2)

        nt_WT_v1_v2.write("\n%d\t %d\t %d\t %d\t %d\t %d\t %.2f\t %.2f\t %.2f\t "  % (method.t_tot, method.t, method.n_sel, n_nt_WT, n_nt_v1, n_nt_v2, avg_E_WT, avg_E_v1, avg_E_v2))

    # breadth

    nt_br.write("method.t is %d " % (method.t))
    br_t.write("\n %d\t%d\t " % (method.t_tot, method.t))
    for_stat.write("\n%d\t %d\t "%(method.t_tot, method.t))
    n_mut.write("\n# method.En(Neut.threshold)\t ratio\t avg_br[l]\t stdev_br[l]\t avg_n_ntAg[l](Avg.#.panel_seq.neut)\t stdev_n_ntAg[l]\t max_n_ntAg[l](Max.#.panel_seq.neut)\t min_n_ntAg[l]\n")

    for l in range(3):

        if method.En[l]==0.:
            break #neut. factor

        for i in range(1,method.n_sel+1):
            method.n_ntAg[i]=0 #initialize

        for i in range(1,method.n_sel+1): #1. for all selected BCs, calc affin b/w each of pan.seq. 3. if Eij> act.thresh * nt.factor # nt_Ag ++
            for j in range(1,method.n_Ag+1):     #method.n_Ag:# panel seq.
                for k in range(1,method.N+1):
                    Eij+=method.h_sel[i,k]*method.s0[j,k]

                if Eij>=(method.En[l]*method.E0):
                    method.n_ntAg[i]+=1  #
                Eij=0.

            if method.n_ntAg[i]>max_n_ntAg:
                max_n_ntAg=method.n_ntAg[i]
                i_max_br=i

            if method.n_ntAg[i]<min_n_ntAg:
                min_n_ntAg=method.n_ntAg[i]
                i_min_br=i

            avg_n_ntAg+=float(method.n_ntAg[i])

            method.br[i]=float(method.n_ntAg[i])/float(method.n_Ag)
            avg_br+=method.br[i]

        avg_n_ntAg/=float(method.n_sel)
        avg_br/=float(method.n_sel)

        for i in range(1, method.n_sel+1):
            var_n_ntAg+=(float(method.n_ntAg[i])-avg_n_ntAg)**2
            var_br+=(method.br[i]-avg_br)**2
            method.br[i]=0.
            method.n_ntAg[i]=0
        var_n_ntAg/=float(method.n_sel)
        var_br/=float(method.n_sel)

        n_mut.write("%.2f\t %.2f\t %.5f\t %.5f\t %.1f\t %.1f\t %d\t %d\n" % (method.En[l]*method.E0, method.En[l], avg_br, math.sqrt(var_br), avg_n_ntAg, math.sqrt(var_n_ntAg), max_n_ntAg, min_n_ntAg))

        if l==2: #method.En=1.05*method.E0=11.3,at last iteraction, write to file
            nt_br.write("%.2f\t %.2f\t %.5f\t %.5f\t %.1f\t %.1f\t %d\t %d\n" % (method.En[l]*method.E0, method.En[l], avg_br, math.sqrt(var_br), avg_n_ntAg, math.sqrt(var_n_ntAg), max_n_ntAg, min_n_ntAg))
            br_t.write("%.3f\t %.3f\t%d\t" % (avg_br, math.sqrt(var_br), method.flag_succ))
            for_stat.write("%.3f\t" % (avg_br))
        avg_br=0.
        var_br=0.
        avg_n_ntAg=0.
        var_n_ntAg=0.
        max_n_ntAg=0
        min_n_ntAg=1000

    # avg method.h for c & v sites
    # only store avg method.h of structural motif type2 method.C
    avg_hc_hv.write("\n %d\t " % (method.t)) # structural motif
    hc_t.write("\n  %d\t %d\t " % (method.t_tot, method.t))
    hv_t.write("\n %d\t %d\t " % (method.t_tot, method.t))

    for l in range(method.n_c_site): #for each residue in conserved region,cal avg h_k & var_h_k over all BCs
        k_c=method.c_site_label[l]

        for i in range(1,method.n_sel+1):
            avg_h_k+=method.h_sel[i,k_c]
        avg_h_k/=float(method.n_sel)

        for i in range(1,method.n_sel+1):
            var_h_k+= (method.h_sel[i,k_c]-avg_h_k)**2
        var_h_k/=float(method.n_sel)
        avg_h_k=0.
        var_h_k=0.

    for l in range(method.n_v_site):
        k_v=method.v_site_label[l]
        for i in range(1,method.n_sel+1):
            avg_h_k+=method.h_sel[i,k_v]
        avg_h_k/=float(method.n_sel)

        for i in range(1,method.n_sel+1):
            var_h_k+=(method.h_sel[i,k_v]-avg_h_k)**2
        var_h_k/=float(method.n_sel)

        if (k_v==4 or k_v==11 or k_v==16 or k_v==38 or k_v==39 or k_v==40) : #if structural motif store it in hc_hv
            avg_hc_hv.write("%.3f\t" % (avg_h_k))
        avg_h_k=0.
        var_h_k=0.

    # number avg & stdev of method.h stored in the order of tp 1 c -> v1 -> v2 -> tp 2 method.C for every 5 rds or method.t==1
    if (method.t_tot % 5==0 or method.t==1): #every 5 rds or on 1st rd
        k_o=0
        for k in range(1,method.N+1):
            method.avg_hk[k]=0.
            method.var_hk[k]=0.  #init.
        #type I c
        for l in range(method.n_c_site):
            k_c=method.c_site_label[l]
            k_o += 1

            for i in range(1,method.n_sel+1):
                method.avg_hk[k_o]+=method.h_sel[i,k_c]
            method.avg_hk[k_o]/= float(method.n_sel)

            for i in range(1, method.n_sel+1):
                method.var_hk[k_o]+=(method.h_sel[i,k_c]-method.avg_hk[k_o])**2
            method.var_hk[k_o]/=float(method.n_sel)

        #v1
        for l in range(method.n_v1_site):
            k_v1=method.v1_label[l]
            k_o += 1

            for i in range(1,method.n_sel+1):
                method.avg_hk[k_o]+=method.h_sel[i,k_v1]
            method.avg_hk[k_o]/=float(method.n_sel)

            for i in range(1,method.n_sel+1):
                method.var_hk[k_o]+=(method.h_sel[i,k_v1]-method.avg_hk[k_o])**2
            method.var_hk[k_o]/=float(method.n_sel)
        #v2
        for l in range(method.n_v2_site):
            k_v2=method.v2_label[l]
            k_o += 1

            for i in range(1,method.n_sel+1):
                method.avg_hk[k_o]+=method.h_sel[i,k_v2]
            method.avg_hk[k_o]/=float(method.n_sel)

            for i in range(1,method.n_sel+1):
                method.var_hk[k_o]+=(method.h_sel[i,k_v2]-method.avg_hk[k_o])**2
            method.var_hk[k_o]/=float(method.n_sel)
        #type II c
        for k in range(41,47):
            k_o += 1
            for i in range(1,method.n_sel+1):
                method.avg_hk[k_o]+=method.h_sel[i,k]
            method.avg_hk[k_o]/=float(method.n_sel)

            for i in range(1,method.n_sel+1):
                method.var_hk[k_o]+=(method.h_sel[i,k]-method.avg_hk[k_o])**2
            method.var_hk[k_o]/=float(method.n_sel)

        avg_hk_t.write("\n %d\t %d\t " % (method.t_tot, method.t))
        std_hk_t.write("\n %d\t %d\t " % (method.t_tot, method.t))

        for k in range(1,method.N+1):
            avg_hk_t.write("%.3f\t " % (method.avg_hk[k]))
            std_hk_t.write("%.3f\t " % (math.sqrt(method.var_hk[k])))
    # site and clone (sel cells) avg method.h
    # for method.h traj
    #avg over the same category of epitopes
    k_o=0

    for l in range(method.n_c_site):
        k_o += 1
        avg_h_c1+=method.avg_hk[k_o]
    avg_h_c1/=float(method.n_c_site)

    for l in range(method.n_v1_site):
        k_o += 1
        avg_h_v1+=method.avg_hk[k_o]
    avg_h_v1/=float(method.n_v1_site)

    for l in range(method.n_v2_site):
        k_o += 1
        avg_h_v2+=method.avg_hk[k_o]
    avg_h_v2/=float(method.n_v2_site)

    for k in range(1,7):
        k_o += 1
        avg_h_c2+=method.avg_hk[k_o]
    avg_h_c2/=6.
    avg_h_c_v_vs_t.write("%d\t %d\t %.3f\t %.3f\t %.3f\t %.3f\n" % (method.t_tot, method.t, avg_h_c1, avg_h_c2, avg_h_v1, avg_h_v2))

    # E hist of sel BCs
    E_sel_Hist.write("\n# method.t=%d\n\n" % (method.t))

    if method.n_round==1:
        for j in range(1,method.n_FDC_Ag+1):
            E_sel_Hist.write("# For Ag j=%d\n" % (j))
            E_sel_hist.E_sel_hist(j,E_sel_Hist)

    if method.n_round>1:
        for j in range(method.n_FDC_Ag+1):
            if j==0: #for WT Ag, make sure method.s value is all 1
                for k in range(1,method.N+1):
                    method.s[j,k]=1
            E_sel_Hist.write("# For Ag j=%d\n" % (j))
            E_sel_hist.E_sel_hist(j,E_sel_Hist)

    # method.h-char of surv BCs
    for i in range(1,method.n_s+1): # h_c, h_v averaged over the same type of sites + averaged over survived BCs
        for l in range(method.n_v_site):
            k_v = method.v_site_label[l]
            method.hv_s[i] += method.h_s[i, k_v]

        for l in range(method.n_c_site):
            k_c= method.c_site_label[l]
            method.hc_s[i] += method.h_s[i, k_c]

        for k in range(1,method.N+1): #initialize method.h_s
            method.h_s[i, k] =0.

        avg_hv_s += method.hv_s[i]/float(method.n_v_site)
        avg_hc_s += method.hc_s[i]/float(method.n_c_site)
        method.hv_s[i]=0.
        method.hc_s[i]=0. # initialize method.hv_s,method.hc_s

    avg_hv_s/=float(method.n_s)
    avg_hc_s/=float(method.n_s)
    method.n_s=0 #initialize method.n_s

    for i in range(1,method.n_sel+1):# h_c, h_v averaged over the same type of sites over: What is this block for?
        for l in range(method.n_v_site):
            k_v=method.v_site_label[l]
            method.hv[i]+=method.h_sel[i,k_v]
        for l in range(method.n_c_site):
            k_c=method.c_site_label[l]
            method.hc[i]+=method.h_sel[i,k_c]
        method.hv[i]/=float(method.n_v_site)
        method.hc[i]/=float(method.n_c_site)
#i
    for i in range(1,method.n_sel+1): #initialize method.hv, method.hc for what???
        method.hv[i]=0.
        method.hc[i]=0.
    # method.h-char of sel BCs
    for i in range(1,method.n_sel+1): # h_c, h_v1,h_v2 averaged over the same type of sites + averaged over selected BCs:
        for l in range(method.n_v_site):
            k_v=method.v_site_label[l]
            method.hv[i]+=method.h_sel[i,k_v]
        for l in range((method.n_v_site-6)/2):
            k_v1=method.v1_label[l]
            k_v2=method.v2_label[l]
            method.hv1[i]+=method.h_sel[i,k_v1]
            method.hv2[i]+=method.h_sel[i,k_v2]
        for l in range(method.n_c_site):
            k_c=method.c_site_label[l]
            method.hc[i]+=method.h_sel[i,k_c]
        for k in range(1,method.N+1):#initialize method.h_sel
            method.h_sel[i,k]=0.

        avg_hv+=method.hv[i]/float(method.n_v_site)
        avg_hv1+=method.hv1[i]/float(((method.n_v_site-6)/2))
        avg_hv2+=method.hv2[i]/float(((method.n_v_site-6)/2))
        avg_hc+=method.hc[i]/float(method.n_c_site)

    avg_hv/=float(method.n_sel)
    avg_hv1/=float(method.n_sel)
    avg_hv2/=float(method.n_sel)
    avg_hc/=float(method.n_sel)

    for i in range(1,method.n_sel+1): #variance of method.hv,method.hv1,method.hv2
        var_hv+=(method.hv[i]/float(method.n_v_site)-avg_hv)**2
        var_hv1+=(method.hv1[i]/float((method.n_v_site-6)/2)-avg_hv1)**2
        var_hv2+=(method.hv2[i]/float((method.n_v_site-6)/2)-avg_hv2)**2
        var_hc+=(method.hc[i]/float(method.n_c_site)-avg_hc)**2

    var_hv/=float(method.n_sel)
    var_hv1/=float(method.n_sel)
    var_hv2/=float(method.n_sel)
    var_hc/=float(method.n_sel)

    sel_BC.write("\n# For %d sel BCs, avg_hv=%.3f, avg_hc=%.3f\n" % (method.n_sel, avg_hv, avg_hc))
    for_stat.write("%.3f\t %.3f\t" %(avg_hv, avg_hc))
    method.n_sel=0 #initialize method.n_sel

    avg_hc_hv.write("%.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t %.3f\t "  % (avg_hc_s, avg_hv_s, avg_hc, avg_hv, avg_hv1, avg_hv2, math.sqrt(var_hc), math.sqrt(var_hv), math.sqrt(var_hv1), math.sqrt(var_hv2)))

    hc_t.write("%.3f\t %.3f\t %d\t " % (avg_hc, math.sqrt(var_hc), method.flag_succ))
    hv_t.write("%.3f\t %.3f\t %d\t " % (avg_hv, math.sqrt(var_hv), method.flag_succ))
    for_stat.close()
