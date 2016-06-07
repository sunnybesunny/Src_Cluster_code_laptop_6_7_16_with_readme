import numpy as np
import random

# # declaration
    #
    # # Define Parameters as macros
    # # method.N= 46        # # sites
    # # method.M= 28        # # v sites
    # # method.Ntot= 511    # total length
    #
    # # method.n_Ag= 141       # # panel seq
    # # method.n_B	= 10000      # max # BCs
    # #
    # # method.N_pr= 2           # # repl before sel
    # # method.N_re= 80          # # cycles
    # # method.N_pr_BB	= 9        # # prol
    # #
    # # method.pA= 0.2       # probab of aa mut
    # # method.pL= 0.3       # probab of lethal mut
    # #
    # # method.frac_rc	= 0.9  # % recycled BCs
    # #
    # # method.Bin_E_sel= 40
    #
    # method.sht	     # shift for dX range
    # # method.R	= 0.3          # scaling hk range
    # #
    # # method.n_v_site	= 28
    # #
    # # method.n_v1_site	= 11
    # # method.n_v2_site	= 11
    # #
    # # method.Ec	= 10.0          # threshold to avoid big weights
    # #
    # # method.all_MC_flag	= 1
    # #
    # # method.f	= 1.1            # method.C[j]*=method.f
    #
    # # Variables
    # method.s_np  # full panel seq
    # method.s0    # binding region
    # method.s          # FDC-Ag seq ,method.s[Ag type,site]
    # method.n_FDC_Ag
    # method.s_im_0  # immunogens
    # method.s_im_1
    # method.s_im_2
    #
    # # epitopes
    # # method.label_np = np.array([97, 276, 278, 279, 280, 282, 283, 365, 368, 371, 372,
    # #  425, 430, 455, 457, 461, 463, 465, 467, 469, 474, 134, 137, 138, 139, 140, 141, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 460, 462, 464, 197,
    # #  234, 276, 386, 392, 463])
    #
    # method.h_GC_0       # method.h val. for sd BCs
    # method.h_GC       # GC BCs
    # method.h_sel    # sel BCs
    # method.h_m   # Memory BCs
    # method.h_s      # surv BCs: those survived from lethal mutations
    #
    # method.avg_hk # # avg method.h val(sel.BCs)method.n
    # method.var_hk
    # method.n_mt
    # # # mut method.s: non lethaly mutated BCs ., sel:selected, m:memory
    # method.n_mt_s
    # method.n_mt_sel
    # method.n_mt_m
    # method.n_mt_b
    # method.n_mt_d  # b:beneficial, d:deleterious
    # method.n_mt_b_s
    # method.n_mt_d_s
    # method.n_mt_b_sel
    # method.n_mt_d_sel
    # method.n_mt_b_m
    # method.n_mt_d_m
    # # method.En = np.array([1.2,1.1,1.05])     # neutralization threshold???
    #
    # method.N_A_aa            # # & label of FDC Ag seen by aa-mut BCs
    # method.Ag_label_aa
    # method.N_A             # # & label of FDC Ag seen by surv BC
    # method.Ag_label
    # method.Ag_label_sel
    # method.flag_aa          # indicator of aa-mut BC 1: mutated
    # # method.v_site_label = np.array([4, 11, 16, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46])
    #
    # # method.v1_label= np.array([4,11,16, 38, 39, 40, 22, 23, 24, 25, 26])
    # # method.v2_label=np.array([27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37])
    # method.c_site_label
    # method.n_c_site
    # method.n_GC_0                     # # init BCs
    # method.n_GC                       # # curr BCs
    # method.n_sel
    # method.n_m
    # method.n_s
    # # method.hl_0=-1.*0.6*0.3     # -0.18
    # # method.hu_0=5.*0.6*0.3     # 0.9
    # #
    # # method.hl=-1.*0.6*method.R         # method.h-range
    # # method.hu=5.*0.6*method.R
    # # method.dxl=-6.5             # dX-range
    # # method.dxu=1.5
    #
    # # method.E0=60.*0.6*method.R*1.0      # 10.8 kcal/mol, act.thresh.
    #
    # method.flag
    # method.C
    # method.act_Ag
    # method.hv    # interaction of sel. BCs
    # method.hv1
    # method.hv2
    # method.hc
    # method.hv_s	# interaction of surv. BCs
    # method.hc_s
    # method.br     # breadth
    # method.n_ntAg     # # neut Ags out of panels
    # method.n            # # bins for Energy histogram
    # method.P	 # method.E_sel hist
    # method.Ej        	#E sel BC. val
    # method.max_mut         # distinct mut
    # method.max_mut_s
    # method.max_mut_sel
    # method.mut_site
    # method.mut_site_s
    # method.mut_site_sel
    # method.W        # weight of all competing cells
    # method.actAg_MC    #  Ag seen by memorgy B cells, randomly assigned
    #
    # # method.E_sel_min=100.
    # # method.E_sel_max=-10.
    # method.E_sel
    # method.h
    # method.E_ij        #never used
    #                   			# # GCR cycle (<80)
    # method.avg_Ps_s           # avg sel probab she doesn'method.t even use this
    # method.seed                   #used for rd # generation
    # method.n_round                     # # of rounds of complete GC cycle on 1 instantiation of GC lp label
    # method.t_tot                   #tot_# of GCRS
    # method.flag_succ              #method.t>1&&method.n_GC>=method.n_GC_0*pow(2,method.N_pr_BB))||method.t==method.N_re
    # method.flag_end                #no modification
    # method.flag_dX                  # mut type 1:affinity +  -1:affinity -   2:shielding motif
    # method.t

# Define Parameters as macros
N= 46        # # sites
M= 28        # # v sites
Ntot= 511    # total length

n_Ag= 141       # # panel seq
n_B	= 10000      # max # BCs

N_pr= 2           # # repl before sel
N_re= 80          # # cycles
N_pr_BB	= 9        # # prol

pA= 0.2       # probab of aa mut
pL= 0.3       # probab of lethal mut

frac_rc	= 0.9  # % recycled BCs

Bin_E_sel= 40

sht	= 0.        # shift for dX range

R	= 0.3          # scaling hk range

n_v_site	= 28

n_v1_site	= 11
n_v2_site	= 11

Ec	= 10.0          # threshold to avoid big weights

all_MC_flag	= 1

f	= 1.1            # C[j]*=f

# Variables
s_np= np.zeros((n_Ag+1,Ntot+1), dtype=int)  # full panel seq
s0= np.zeros((n_Ag+1,N+1), dtype=int )     # binding region

s= np.zeros((6,N+1), dtype=int)        # FDC-Ag seq ,s[Ag type,site]
n_FDC_Ag=0

s_im_0= np.zeros((N+1,), dtype=int)  # immunogens
s_im_1= np.zeros((N+1,), dtype=int)
s_im_2= np.zeros((N+1,), dtype=int)

# epitopes
label_np = np.array([97, 276, 278, 279, 280, 282, 283, 365, 368, 371, 372,
 425, 430, 455, 457, 461, 463, 465, 467, 469, 474, 134, 137, 138, 139, 140, 141, 143, 144, 145, 146, 147, 148, 149, 150, 151, 152, 460, 462, 464, 197,
 234, 276, 386, 392, 463], dtype=int)

h_GC_0= np.zeros((5,N+1), dtype=float)        # h val. for sd BCs

h_GC = np.zeros((n_B+1,N+1), dtype=float)      # GC BCs
h_sel= np.zeros((n_B+1,N+1), dtype=float)     # sel BCs
h_m= np.zeros((n_B*10+1,N+1), dtype=float)   # Memory BCs

h_s= np.zeros((n_B+1,N+1), dtype=float )      # surv BCs: those survived from lethal mutations

avg_hk= np.zeros((N+1,), dtype=float) # # avg h val(sel.BCs)n
var_hk= np.zeros((N+1,), dtype=float)

n_mt= np.zeros((n_B+1,), dtype=int)
# # mut s: non lethaly mutated BCs ., sel:selected, m:memory
n_mt_s= np.zeros((n_B+1,), dtype=int)
n_mt_sel= np.zeros((n_B+1,), dtype=int)
n_mt_m= np.zeros((n_B*10+1,), dtype=int)

n_mt_b= np.zeros((n_B+1,), dtype=int)
n_mt_d= np.zeros((n_B+1,), dtype=int)  # b:beneficial, d:deleterious
n_mt_b_s= np.zeros((n_B+1,), dtype=int)
n_mt_d_s= np.zeros((n_B+1,), dtype=int)
n_mt_b_sel= np.zeros((n_B+1,), dtype=int)
n_mt_d_sel= np.zeros((n_B+1,), dtype=int)
n_mt_b_m= np.zeros((n_B*10+1,), dtype=int)
n_mt_d_m= np.zeros((n_B*10+1,), dtype=int)
En = np.array([1.2,1.1,1.05], dtype=float)     # neutralization threshold???

N_A_aa= np.zeros((n_B+1,), dtype=int)            # # & label of FDC Ag seen by aa-mut BCs
Ag_label_aa =np.zeros((n_B+1,6), dtype=int)

N_A= np.zeros((n_B+1,), dtype=int)               # # & label of FDC Ag seen by survived BC( aa-mut + unmut BCs)
Ag_label= np.zeros((n_B+1,6), dtype=int)

Ag_label_sel= np.zeros((n_B+1,3), dtype=int)

flag_aa= np.zeros((n_B+1,), dtype=int)           # indicator of aa-mut BC 1: mutated
v_site_label = np.array([4, 11, 16, 22, 23, 24, 25, 26, 27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40, 41, 42, 43, 44, 45, 46], dtype=int)

v1_label= np.array([4,11,16, 38, 39, 40, 22, 23, 24, 25, 26], dtype=int)
v2_label=np.array([27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37], dtype=int)
c_site_label = np.zeros((N,), dtype=int)

n_c_site=0
n_GC_0=0                      # # init BCs
n_GC=0                        # # curr BCs
n_sel=0
n_m=0
n_s=0

hl_0=-1.*0.6*0.3     # -0.18
hu_0=5.*0.6*0.3     # 0.9

hl=-1.*0.6*R         # h-range
hu=5.*0.6*R
dxl=-6.5             # dX-range
dxu=1.5

E0=60.*0.6*R*1.0      # 10.8 kcal/mol, act.thresh.

flag= np.zeros((n_B+1,), dtype=int)
C = np.zeros((6,), dtype=float)

act_Ag= np.zeros((n_B+1,N+1), dtype=int)

hv= np.zeros((n_B+1,), dtype=float)      # interaction of sel. BCs
hv1= np.zeros((n_B+1,), dtype=float)
hv2= np.zeros((n_B+1,), dtype=float)
hc= np.zeros((n_B+1,), dtype=float)
hv_s= np.zeros((n_B+1,), dtype=float)	# interaction of surv. BCs
hc_s= np.zeros((n_B+1,), dtype=float)

br= np.zeros((n_B+1,), dtype=float)       # breadth
n_ntAg= np.zeros((n_B+1,), dtype=int)       # # neut Ags out of panels

n=np.zeros((Bin_E_sel+1,), dtype=int)              # # bins for Energy histogram
P=np.zeros((Bin_E_sel+1,), dtype=float)	 # E_sel hist
Ej= np.zeros((n_B+1,), dtype=float )        	#E sel BC. val

max_mut= np.zeros((n_B+1,), dtype=int)          # distinct mut
max_mut_s= np.zeros((n_B+1,), dtype=int)
max_mut_sel= np.zeros((n_B+1,), dtype=int)

mut_site= np.zeros((n_B+1,N+1), dtype=int) # mutated site of B cell?
mut_site_s= np.zeros((n_B+1,N+1), dtype=int) # mutated site of survived B cell?
mut_site_sel= np.zeros((n_B+1,N+1), dtype=int) # mutated site of selected B cell?

W= np.zeros((n_B*10+1,), dtype=float)            # weight of all competing cells
actAg_MC= np.zeros((n_B*10+1,), dtype=int)    #  Ag seen by memorgy B cells, randomly assigned

E_sel_min=100.
E_sel_max=-10.
E_sel= np.zeros((n_B+1,), dtype=float)

h= np.zeros((N+1,), dtype=float)
E_ij = np.zeros((6,), dtype=float)            #never used
                  			# # GCR cycle (<80)
avg_Ps_s=0.                  # avg sel probab she doesn't even use this

seed=100                       #used for rd # generation
n_round=0                       # # of rounds of complete GC cycle on 1 instantiation of GC lp label
t_tot=0                       #tot_# of GCRS
flag_succ=0                  #t>1&&n_GC>=n_GC_0*pow(2,N_pr_BB))||t==N_re
flag_end=0                  #no modification
flag_dX=0                       # mut type 1:affinity +  -1:affinity -   2:shielding motif
t=0




def rand_gen():

    floatinfo = np.finfo(float)
    epsilon = floatinfo.eps
    a = np.random.uniform(0+epsilon, 1+epsilon)
    return a
