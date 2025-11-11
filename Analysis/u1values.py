import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
import FormFactor as FF
import csv
from scipy.integrate import quad
from scipy.integrate import simpson

plt.style.use('sty.mplstyle')
path_effects = [pe.Stroke(linewidth=3, foreground='k'), pe.Normal()]
path_effects1 = [pe.Stroke(linewidth=2, foreground='k'), pe.Normal()]

lw = 3
colors = {'Exp':'#DDDDDD', '3pF':'#2E2585', 'Helm':'#337538', 'KN':'#5DA899', 'adKN':'#94CBEC',
          'HF-SkE2':'#DCCD7D', 'Barbieri':'#C26A77', 'Payne':'#9F4A96'} #7E2954

#=================================================
#======   Generate u1 from data (** OK **)   =====
#=================================================
def Generateu1(Q2ref, Q, FF):
    x = Q**2	# Momentum transfer squared
    y = FF**2	# Form Factor squared
    u1 = np.ones_like(Q) 		# Make plots with this data
    u1Cpp = np.ones_like(Q2ref)	# Compute events using this data
    # Generate the u1 values using the simpson method
    for i in range(len(Q)):
        u1[i] = simpson(y[i::], x=x[i::])
    # Generate the u1 values from interpolation
    for i in range(len(Q2ref)):
        q2_ref = Q2ref[i]
        if ( q2_ref > x.max() ):
            u1Cpp[i] = 0
        elif ( x.min() <= q2_ref ):
            x1 = x[x <= q2_ref][-1]
            y1 = u1[x <= q2_ref][-1]
            x2 = x[x > q2_ref][0]
            y2 = u1[x > q2_ref][0]
            slope = ( y2 - y1 )/( x2 - x1 )
            u1Cpp[i] = slope*(q2_ref - x1) + y1
        else:
            x1, x2 = x[0], x[1]
            y1, y2 = u1[0], u1[1]
            slope = ( y2 - y1 )/ ( x2 - x1 )
            u1Cpp[i] = slope*(q2_ref - x1) + y1
    return u1Cpp, u1 

#=================================================
#=======   u1 values for Argon (** OK **)   ======
#=================================================
def Ar403pF_FF_Squared(Q2):
	return FF.Ar403pF_FF(Q2)**2

def Ar403pF_u1(Q2):
	main, _ = quad(Ar403pF_FF_Squared, Q2, np.inf)
	return main
Ar403pF_u1 = np.vectorize(Ar403pF_u1)
#=================================================
def Ar40HELM_FF_Squared(Q2):
	return FF.Ar40HELM_ff(Q2)**2

def Ar40HELM_u1(Q2):
	main, _ = quad(Ar40HELM_FF_Squared, Q2, np.inf)
	return main
Ar40HELM_u1 = np.vectorize(Ar40HELM_u1)
#=================================================
def Ar40KN_FF_Squared(Q2):
	return FF.Ar40KN_ff(Q2)**2

def Ar40KN_u1(Q2):
	main, _ = quad(Ar40KN_FF_Squared, Q2, np.inf)
	return main
Ar40KN_u1 = np.vectorize(Ar40KN_u1)
#=================================================
def Ar40adKN_FF_Squared(Q2):
	return FF.Ar40adKN_ff(Q2)**2

def Ar40adKN_u1(Q2):
	main, _ = quad(Ar40adKN_FF_Squared, Q2, np.inf)
	return main
Ar40adKN_u1 = np.vectorize(Ar40adKN_u1)

#=================================================
#=====   u1 values for Tungsten (** OK **)   =====
#=================================================
def W184WS_FF_Squared(Q2):
	return FF.W184WS_FF(Q2)**2

def W184WS_u1(Q2):
	main, _ = quad(W184WS_FF_Squared, Q2, np.inf)
	return main
W184WS_u1 = np.vectorize(W184WS_u1)
#=================================================
def W184HELMv1_FF_Squared(Q2):
	return FF.W184HELMv1_ff(Q2)**2
	
def W184HELMv1_u1(Q2):
	main, _ = quad(W184HELMv1_FF_Squared, Q2, np.inf)
	return main
W184HELMv1_u1 = np.vectorize(W184HELMv1_u1)
#=================================================
def W184HELMv2_FF_Squared(Q2):
	return FF.W184HELMv2_ff(Q2)**2
	
def W184HELMv2_u1(Q2):
	main, _ = quad(W184HELMv2_FF_Squared, Q2, np.inf)
	return main
W184HELMv2_u1 = np.vectorize(W184HELMv2_u1)
#=================================================
def W184KN_FF_Squared(Q2):
	return FF.W184KN_ff(Q2)**2

def W184KN_u1(Q2):
	main, _ = quad(W184KN_FF_Squared, Q2, np.inf)
	return main
W184KN_u1 = np.vectorize(W184KN_u1)
#=================================================
def W184adKN_FF_Squared(Q2):
	return FF.W184adKN_ff(Q2)**2

def W184adKN_u1(Q2):
	main, _ = quad(W184adKN_FF_Squared, Q2, np.inf)
	return main
W184adKN_u1 = np.vectorize(W184adKN_u1)

#=================================================
#=======   load (q2,u1) values (** OK **)   ======
#=================================================
Ar403pF_Q2Alt = np.loadtxt('../csv/InterpolationList/Argon/40Ar_InterpolationList3pF_Altmannshofer.txt', delimiter=',', usecols=0)
Ar403pF_u1Alt = np.loadtxt('../csv/InterpolationList/Argon/40Ar_InterpolationList3pF_Altmannshofer.txt', delimiter=',', usecols=1)

Ar40HFSKE2_QVis = FF.GeVem1_to_fm*np.loadtxt('../csv/FormFactors/40Ar_ChFF_HF-SkE2.dat', usecols=1) # Momentum transfer from fm^{-1} to GeV
Ar40HFSKE2_FFVis = np.abs( np.loadtxt('../csv/FormFactors/40Ar_ChFF_HF-SkE2.dat', usecols=10) )		# Charge Form Factor

Ar40NNLO_Barbieri_Q = np.loadtxt('../csv/FormFactors/40Ar_ChFF_NNLOsat_Barbieri.dat', usecols=0)
Ar40NNLO_Barbieri_FF = np.loadtxt('../csv/FormFactors/40Ar_ChFF_NNLOsat_Barbieri.dat', usecols=1)

Ar40NNLO_Payne_Q = 0.1973*np.loadtxt('../csv/FormFactors/Fch_40Ar_N2LO_Tminus1_Payne.dat', usecols=0)
Ar40NNLO_Payne_FF = np.loadtxt('../csv/FormFactors/Fch_40Ar_N2LO_Tminus1_Payne.dat', usecols=1)

W184WS_Q2Diego = np.loadtxt('../csv/InterpolationList/Tungsten/184W_InterpolationListWS_Diego.txt', delimiter=',', usecols=0)
W184WS_u1Diego = np.loadtxt('../csv/InterpolationList/Tungsten/184W_InterpolationListWS_Diego.txt', delimiter=',', usecols=1)

#=================================================
#========   Generate the data (** OK **)   =======
#=================================================
Ar403pF_u1Christian = Ar403pF_u1(Ar403pF_Q2Alt)
Ar40HELM_u1Christian = Ar40HELM_u1(Ar403pF_Q2Alt)
Ar40KN_u1Christian = Ar40KN_u1(Ar403pF_Q2Alt)
Ar40adKN_u1Christian = Ar40adKN_u1(Ar403pF_Q2Alt)
Ar40HFSKE2_u1Christian, Ar40HFSKE2_u1 = Generateu1(Ar403pF_Q2Alt, Ar40HFSKE2_QVis, Ar40HFSKE2_FFVis)
Ar40NNLO_Barbieri_u1Chris, Ar40NNLO_Barbieri_u1 = Generateu1(Ar403pF_Q2Alt, Ar40NNLO_Barbieri_Q, Ar40NNLO_Barbieri_FF)
Ar40NNLO_Payne_u1Chris, Ar40NNLO_Payne_u1 = Generateu1(Ar403pF_Q2Alt, Ar40NNLO_Payne_Q, Ar40NNLO_Payne_FF)

W184WS_u1Christian = W184WS_u1(W184WS_Q2Diego)
W184HELMv1_u1Christian = W184HELMv1_u1(W184WS_Q2Diego)
W184HELMv2_u1Christian = W184HELMv2_u1(W184WS_Q2Diego)
W184KN_u1Christian = W184KN_u1(W184WS_Q2Diego)
W184adKN_u1Christian = W184adKN_u1(W184WS_Q2Diego)

#=================================================
#=========   Writing to CSV (** OK **)   =========
#=================================================
with open("../csv/InterpolationList/Argon/40Ar_InterpolationList_3pF_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t Three-parameter Fermi Model (3pF) Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(Ar403pF_Q2Alt, Ar403pF_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Argon/40Ar_InterpolationList_Helm_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t Helm Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(Ar403pF_Q2Alt, Ar40HELM_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Argon/40Ar_InterpolationList_KN_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t Klein-Nystrand Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(Ar403pF_Q2Alt, Ar40KN_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Argon/40Ar_InterpolationList_adKN_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t (ad.) Klein-Nystrand Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(Ar403pF_Q2Alt, Ar40adKN_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Argon/40Ar_InterpolationList_HFSKE2_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t HF-SKE2 Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(Ar403pF_Q2Alt, Ar40HFSKE2_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Argon/40Ar_InterpolationList_NNLOBarbieri_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t NNLOsat-Barbieri Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(Ar403pF_Q2Alt, Ar40NNLO_Barbieri_u1Chris):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Argon/40Ar_InterpolationList_NNLOPayne_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t NNLOsat-Payne Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(Ar403pF_Q2Alt, Ar40NNLO_Payne_u1Chris):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Tungsten/184W_InterpolationListWS_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t Woods-Saxon Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(W184WS_Q2Diego, W184WS_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Tungsten/184W_InterpolationListHelmv1_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t Helm [v1] Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(W184WS_Q2Diego, W184HELMv1_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Tungsten/184W_InterpolationListHelmv2_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t Helm [v2] Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(W184WS_Q2Diego, W184HELMv2_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Tungsten/184W_InterpolationListKN_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t Klein-Nystrand Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(W184WS_Q2Diego, W184KN_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")

with open("../csv/InterpolationList/Tungsten/184W_InterpolationListadKN_Christian.txt", mode='w', newline="") as file:
	file.write("########################################################################################## \n")
	file.write(f"#####\t (ad.) Klein-Nystrand Interpolation List [Christian] \n")
	file.write("########################################################################################## \n")
	for q2val, u1val in zip(W184WS_Q2Diego, W184adKN_u1Christian):
		file.write(f"\t {q2val:.17f}, \t {u1val:.30f}, \n")


#=================================================
#=======   Plotting u1 values (** OK **)   =======
#=================================================
Ar40u1_fig = plt.figure()
Ar40FF_fig = plt.figure()
W184u1_fig = plt.figure()

Ar40u1_ax = Ar40u1_fig.add_subplot()
Ar40FF_ax = Ar40FF_fig.add_subplot()
W184u1_ax = W184u1_fig.add_subplot()

Qmax = np.sqrt(2.6) # GeV
Qval = np.linspace(0, Qmax, 500) # GeV

Ar40u1_ax.plot(Ar403pF_Q2Alt, Ar403pF_u1Alt, ls='-', lw=lw, color=colors['3pF'], label='3pF', path_effects=path_effects)
#Ar40u1_ax.plot(Ar403pF_Q2Alt, Ar403pF_u1Christian, lw=lw, color='#ff7f0e', path_effects=path_effects)
Ar40u1_ax.plot(Ar403pF_Q2Alt, Ar40HELM_u1Christian, ls='--', lw=lw, color=colors['Helm'], label='Helm', path_effects=path_effects)
Ar40u1_ax.plot(Ar403pF_Q2Alt, Ar40KN_u1Christian, ls=':', lw=lw, color=colors['KN'], label='KN', path_effects=path_effects)
Ar40u1_ax.plot(Ar403pF_Q2Alt, Ar40adKN_u1Christian, ls='-.', lw=lw, color=colors['adKN'], label='(ad.) KN', path_effects=path_effects)
Ar40u1_ax.plot(Ar403pF_Q2Alt, Ar40HFSKE2_u1Christian, ls=(5,(10,3)), lw=lw, color=colors['HF-SkE2'], label='HF-SkE2', path_effects=path_effects)
Ar40u1_ax.plot(Ar403pF_Q2Alt, Ar40NNLO_Barbieri_u1Chris, ls=(0,(2,2)), lw=lw, color=colors['Barbieri'], label='NNLOsat Barbieri', path_effects=path_effects)
Ar40u1_ax.plot(Ar403pF_Q2Alt, Ar40NNLO_Payne_u1Chris, ls='-', lw=lw, color=colors['Payne'], label='NNLOsat Payne', path_effects=path_effects)

Ar40u1_ax.set_ylim([1e-13, 1e-1])
Ar40u1_ax.set_xlim([0.0, 1.0])
Ar40u1_ax.set_yscale('log')
Ar40u1_ax.set_xlabel(r'$Q^2$ [GeV$^2$]')
Ar40u1_ax.set_ylabel(r'$u_1$ [GeV$^2$]')
Ar40u1_ax.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
Ar40u1_ax.yaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
Ar40u1_ax.legend(fontsize=20, frameon=False, framealpha=1, edgecolor='k', fancybox=True)
Ar40u1_ax.tick_params(which='both', right=False, top=False, tickdir='in')
Ar40u1_ax.set_title(r'Argon $^{40}$Ar $u_1$ values', pad = 20)
Ar40u1_fig.savefig('../Plots/Figures_PDF/rho_FF_u1/40Argon_u1vals.pdf', transparent=False, bbox_inches='tight')
Ar40u1_fig.savefig('../Plots/Figures_PNG/rho_FF_u1/40Argon_u1vals.png', transparent=False, bbox_inches='tight')

Ar40FF_ax.plot(Qval, np.abs(FF.Ar403pF_FF(Qval**2)), ls='-', lw=lw, color=colors['3pF'], label='3pF', path_effects=path_effects1)
Ar40FF_ax.plot(Qval, np.abs(FF.Ar40HELM_FF(Qval**2)), ls='--', lw=lw, color=colors['Helm'], label='Helm', path_effects=path_effects1)
Ar40FF_ax.plot(Qval, np.abs(FF.Ar40KN_FF(Qval**2)), ls=':', lw=lw, color=colors['KN'], label='KN', path_effects=path_effects1)
Ar40FF_ax.plot(Qval, np.abs(FF.Ar40adKN_FF(Qval**2)), ls='-.', lw=lw, color=colors['adKN'], label='(ad.) KN', path_effects=path_effects1)
Ar40FF_ax.plot(Ar40HFSKE2_QVis, Ar40HFSKE2_FFVis, ls=(5,(10,3)), lw=lw, color=colors['HF-SkE2'], label='HF-SkE2', path_effects=path_effects1)
Ar40FF_ax.plot(Ar40NNLO_Barbieri_Q, Ar40NNLO_Barbieri_FF, ls=(0,(2,2)), lw=lw, color=colors['Barbieri'], label='NNLOsat Barbieri', path_effects=path_effects1)
Ar40FF_ax.plot(Ar40NNLO_Payne_Q, Ar40NNLO_Payne_FF, ls='-', lw=lw, color=colors['Payne'], label='NNLOsat Payne', path_effects=path_effects1)
Ar40FF_ax.set_yscale('log')
Ar40FF_ax.set_ylim([1e-9, 1])
Ar40FF_ax.set_xlim([0, Qval.max()])
Ar40FF_ax.legend(fontsize=20, frameon=False, framealpha=1, edgecolor='k', fancybox=True)
Ar40FF_ax.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
Ar40FF_ax.yaxis.grid(True, ls='-', which='both', color='gray', alpha=0.35)
Ar40FF_ax.tick_params(which='both', right=False, top=False, tickdir='in')
Ar40FF_ax.set_xlabel(r'Momentum Transfer $Q$ [GeV]')
Ar40FF_ax.set_ylabel(r'Form Factor $|F(Q^2)|$')
Ar40FF_ax.set_title(r'Argon $^{40}$Ar Form Factor')
Ar40FF_fig.savefig('../Plots/Figures_PDF/rho_FF_u1/40ArgonFF_wholerange.pdf', transparent=False, bbox_inches='tight')
Ar40FF_fig.savefig('../Plots/Figures_PNG/rho_FF_u1/40ArgonFF_wholerange.png', transparent=False, bbox_inches='tight')


W184u1_ax.plot(W184WS_Q2Diego, W184WS_u1Diego, color='black', lw=5, label='WS [Diego C++]')
W184u1_ax.plot(W184WS_Q2Diego, W184WS_u1Christian, color='blue', lw=2, label='WS [Christian]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
W184u1_ax.plot(W184WS_Q2Diego, W184HELMv1_u1Christian, color='red', lw=2, label='Helm [v1]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
W184u1_ax.plot(W184WS_Q2Diego, W184HELMv2_u1Christian, color='green', lw=2, label='Helm [v2]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
W184u1_ax.plot(W184WS_Q2Diego, W184KN_u1Christian, color='purple', lw=2, label='KN', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
W184u1_ax.plot(W184WS_Q2Diego, W184adKN_u1Christian, color='orange', lw=2, label='(ad.) KN', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
W184u1_ax.set_yscale('log')
W184u1_ax.set_xlim([0, W184WS_Q2Diego.max()])
W184u1_ax.legend(fontsize=20, loc='lower left', frameon=False, framealpha=1, edgecolor='k', fancybox=True)
W184u1_ax.set_xlabel(r'$Q^2$ [GeV$^2$]')
W184u1_ax.set_ylabel(r'$u_1$ [GeV$^2$]')
W184u1_ax.set_title(r'Tungsten $^{184}$W $u_1$ values', pad=20)
W184u1_fig.savefig('../Plots/Figures_PDF/rho_FF_u1/184Tungsten_u1vals.pdf', transparent=False, bbox_inches='tight')
W184u1_fig.savefig('../Plots/Figures_PNG/rho_FF_u1/184Tungsten_u1vals.png', transparent=False, bbox_inches='tight')

plt.show()
