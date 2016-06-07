import method
import proliferation_SHM
import calc_all_Wi
import E_sel_hist
import GCR_initiation
import output_reset
import pmut
import removeLethal
import selection_mixture_encounter
#import pdb
from method import rand_gen

def differentiation_recycle(succ_ts, nt_WT_v1_v2, avg_hc_hv, meanE_vs_t, output_ns): # 1.initialize method.h_GC,method.n_GC=0,2.if recycle, method.n_GC++ then
#transfer method.h,method.n_mt_b,method.n_mt_d,method.mut_site 3.if not method.n_m++ then transfer method.h,method.n_mt_b,method.n_mt_d,method.mut_site
#4. end method.flag if method.n_GC=0 & method.t> 1GCR

    Ei=0.
    Ar=0.
    for i in range(1, method.n_GC+1):
        for k in range(1,method.N+1): #method.h_s has method.h_GC
            method.h_GC[i,k]=0.

    method.n_GC=0  #method.n_GC is now initialized

    for i in range(1, method.n_sel+1):#method.n_sel: # sel after Tfh
        rd=method.rand_gen()
        if rd <= method.frac_rc: #method.frac_rc= 0.9, if recycles,append method.h_GC,method.n_mt,method.n_mt_b,method.n_mt_d,method.max_mut
            method.n_GC += 1
            for k in range(1, method.N+1):
                method.h_GC[method.n_GC,k]=method.h_sel[i,k]
            method.n_mt[method.n_GC]=method.n_mt_sel[i]
            method.n_mt_b[method.n_GC]=method.n_mt_b_sel[i]
            method.n_mt_d[method.n_GC]=method.n_mt_d_sel[i]
            method.max_mut[method.n_GC]=method.max_mut_sel[i]

            for l in range(1,method.N+1):  # transfer method.mut_site from method.mut_site_sel
                if method.mut_site_sel[i,l]==0: break
                method.mut_site[method.n_GC,l]=method.mut_site_sel[i,l]
        else:  # if not recycled becomes m BCs
            method.n_m += 1
            for k in range(1,method.N+1):  # store method.h_sel,method.n_mt_sel,method.n_mt_b_sel,method.n_mt_d_sel in equivalent m matrix
                method.h_m[method.n_m, k]=method.h_sel[i, k]
            method.n_mt_m[method.n_m] = method.n_mt_sel[i]
            method.n_mt_b_m[method.n_m] = method.n_mt_b_sel[i]
            method.n_mt_d_m[method.n_m] = method.n_mt_d_sel[i]

    if (method.t > 1 and method.n_GC ==0): # if >1 GCR && method.n_GC ==0, ends GC
        succ_ts.write("%d %d %d\n" % (method.t_tot, method.t, method.flag_succ))
        nt_WT_v1_v2.write("%d\n" % (method.flag_succ))
        avg_hc_hv.write("%d\n" % (method.flag_succ))
        meanE_vs_t.write("%d\n" % (method.flag_succ))
        method.flag_end=1

    output_ns.write("%d\t %d\t %d\t %.3f\t %.4f\t" % (method.n_s, method.n_sel, method.n_GC, float(method.n_sel)/float(method.n_s), method.C[1]))

    if method.n_sel>0:
        print "A frac of %d/%d=%.3f sel B cells are recycled\n" % (method.n_GC, method.n_sel, float(method.n_GC)/float(method.n_sel))
        print "A frac of %d/%d=%.3f become identical MCs and PCs\n" % (method.n_sel-method.n_GC, method.n_sel,(float(method.n_sel-method.n_GC)/float(method.n_sel)))
        print "# newly gen MCs and PCs is %d\n" % (method.n_m)
    else:
        print "Zero selected B cells"
        print "# newly gen MCs and PCs is %d\n" % (method.n_m)

