import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.integrate import quad
import FormFactor as FF

plt.style.use('sty.mplstyle')
path_effects = [pe.Stroke(linewidth=3.4, foreground='k'), pe.Normal()]

lw = 3
colors = {'Exp':'#DDDDDD', '3pF':'#2E2585', 'Helm':'#337538', 'KN':'#5DA899', 'adKN':'#94CBEC',
          'HF-SkE2':'#DCCD7D', 'Barbieri':'#C26A77', 'Payne':'#9F4A96'} #7E2954

#=================================================
#=   Extract Form Factor Ottermann  (** OK **)   =
#=================================================
Ar40_E1_Ottermann = 0.11575 # GeV
Ar40_E2_Ottermann = 0.24903 # GeV
Ar40_Theta1_Ottermann = np.array([55.017, 60.012, 65.008, 70.005, 75.021, 80.015, 85.020, 90.026, 95.040], dtype=np.float64)
Ar40_Theta2_Ottermann = np.array([30.019, 35.010, 40.010, 45.009, 50.012, 55.017, 60.012], dtype=np.float64)
Ar40_Exp_Xsection1_Ottermann = 1e-5*np.array([6302, 3417, 1867, 1027, 553.3, 295.1, 158.6, 82.81, 42.92], dtype=np.float64) # fm^{2}
Ar40_Exp_Xsection2_Ottermann = 1e-5*np.array([9129, 2310, 536.7, 100.0, 13.68, 2.963, 3.041], dtype=np.float64) # fm^{2}

#=================================================
#==   Extract numeric Form Factor  (** OK **)   ==
#=================================================
Ar40_HFSKE2_Qvals = FF.GeVem1_to_fm*np.loadtxt('../csv/FormFactors/40Ar_ChFF_HF-SkE2.dat', usecols=1) # Momentum transfer in GeV
Ar40_HFSKE2_FFvals = np.abs( np.loadtxt('../csv/FormFactors/40Ar_ChFF_HF-SkE2.dat', usecols=10) )  # Absolute value of the FF

Ar40_NNLO_Barbieri_Qvals = np.loadtxt('../csv/FormFactors/40Ar_ChFF_NNLOsat_Barbieri.dat', usecols=0)		# momentum transfer (GeV)
Ar40_NNLO_Barbieri_FFvals = np.loadtxt('../csv/FormFactors/40Ar_ChFF_NNLOsat_Barbieri.dat', usecols=1)	# form factor values

Ar40_NNLO_Payne_Qvals = 0.1973*np.loadtxt('../csv/FormFactors/Fch_40Ar_N2LO_Tminus1_Payne.dat', usecols=0)	# momentum transfer (GeV)
Ar40_NNLO_Payne_FFvals = np.loadtxt('../csv/FormFactors/Fch_40Ar_N2LO_Tminus1_Payne.dat', usecols=1)	# form factor values

#=================================================
#=======   Generating the Data (** OK **)   ======
#=================================================
Qmax = 0.7 # GeV
Qmax_incoh = 3.0 # GeV
Qvals = np.linspace(0, Qmax, 500)
Qvals_incoh = np.linspace(0, Qmax_incoh, 500)

Ar40_Q1_Ottermann, Ar40_FF1_Ottermann = FF.Compute_Exp_FF(Ar40_E1_Ottermann, Ar40_Theta1_Ottermann, Ar40_Exp_Xsection1_Ottermann)
Ar40_Q2_Ottermann, Ar40_FF2_Ottermann = FF.Compute_Exp_FF(Ar40_E2_Ottermann, Ar40_Theta2_Ottermann, Ar40_Exp_Xsection2_Ottermann)
Ar40_Q_Ottermann = np.concatenate( (Ar40_Q1_Ottermann, Ar40_Q2_Ottermann) )
Ar40_FF_Ottermann = np.concatenate( (Ar40_FF1_Ottermann, Ar40_FF2_Ottermann) )
Ar40FBSE_FF_vals = np.abs( FF.Ar40FBSE_FF(Qvals**2) ) 	# Absolute value of the [40Ar] Fourier-Bessel Form Factor
Ar403pF_FF_vals = np.abs( FF.Ar403pF_FF(Qvals**2) ) 	# Absolute value of the [40Ar] three-parameter Fermi Form Factor
Ar40HELM_FF_vals = np.abs( FF.Ar40HELM_FF(Qvals**2) ) 	# Absolute value of the [40Ar] Helm Form Factor
Ar40KN_FF_vals = np.abs( FF.Ar40KN_FF(Qvals**2) ) 		# Absolute value of the [40Ar] KN Form Factor
Ar40adKN_FF_vals = np.abs( FF.Ar40adKN_FF(Qvals**2) ) 	# Absolute value of the [40Ar] (ad.) KN Form Factor


W184WS_FF_values = np.abs( FF.W184WS_FF(Qvals**2) )			# Absolute value of the [184W] Woods-Saxon Form Factor
W184HELMv1_FF_values = np.abs( FF.W184HELMv1_FF(Qvals**2) )	# Absolute value of the [184W] Helmv1 Form Factor
W184HELMv2_FF_values = np.abs( FF.W184HELMv2_FF(Qvals**2) )	# Absolute value of the [184W] Helmv2 Form Factor
W184KN_FF_values = np.abs( FF.W184KN_FF(Qvals**2) )			# Absolute value of the [184W] KN Form Factor
W184adKN_FF_values = np.abs( FF.W184adKN_FF(Qvals**2) )		# Absolute value of the [184W] (ad.) KN Form Factor

PROTONGe_Alt_values = np.abs( FF.PROTONGe_Alt(Qvals_incoh**2) )	# Absolute value of the [proton] Alt. Ge Form Factor
PROTONGm_Alt_values = np.abs( FF.PROTONGm_Alt(Qvals_incoh**2) )	# Absolute value of the [proton] Alt. Gm Form Factor
NEUTRONGm_Alt_values = np.abs( FF.NEUTRONGm_Alt(Qvals_incoh**2) )	# Absolute value of the [neutron] Alt. Gm Form Factor
NEUTRONGe_Alt_values = np.abs( FF.NEUTRONGe_Alt(Qvals_incoh**2) )	# Absolute value of the [neutron] Alt. Ge Form Factor

PROTONGe_Kelly_values = np.abs( FF.PROTONGe_Kelly(Qvals_incoh**2) )	# Absolute value of the [proton] Kelly Ge Form Factor
PROTONGm_Kelly_values = np.abs( FF.PROTONGm_Kelly(Qvals_incoh**2) )	# Absolute value of the [proton] Kelly Gm Form Factor
NEUTRONGm_Kelly_values = np.abs( FF.NEUTRONGm_Kelly(Qvals_incoh**2) )	# Absolute value of the [neutron] Kelly Gm Form Factor
NEUTRONGe_Kelly_values = np.abs( FF.NEUTRONGe_Kelly(Qvals_incoh**2) )	# Absolute value of the [neutron] Kelly Ge Form Factor

#=================================================
#====   Plotting the Form Factors (** OK **)   ===
#=================================================
Ar40_fig = plt.figure()
W184_fig = plt.figure()
Proton_fig = plt.figure()
Neutron_fig = plt.figure()

Ar40_ax = Ar40_fig.add_subplot()
W184_ax = W184_fig.add_subplot()
Proton_ax = Proton_fig.add_subplot()
Neutron_ax = Neutron_fig.add_subplot()

Ar40_ax.plot(Ar40_Q_Ottermann, Ar40_FF_Ottermann, ls='', marker='s', markersize=2, color=colors['Exp'], path_effects=path_effects, label='Exp', zorder=100)
Ar40_ax.plot(Qvals, Ar403pF_FF_vals, color=colors['3pF'], ls='-', lw=lw, label="3pF", path_effects=path_effects)
Ar40_ax.plot(Qvals, Ar40HELM_FF_vals, color=colors['Helm'], ls="--", lw=lw, label="Helm", path_effects=path_effects)
Ar40_ax.plot(Qvals, Ar40KN_FF_vals, color=colors['KN'], ls=":", lw=lw, label="KN", path_effects=path_effects)
Ar40_ax.plot(Qvals, Ar40adKN_FF_vals, color=colors['adKN'], ls="-.", lw=lw, label="(ad.) KN", path_effects=path_effects)
Ar40_ax.plot(Ar40_HFSKE2_Qvals, Ar40_HFSKE2_FFvals, color=colors['HF-SkE2'], ls=(5,(10,3)), lw=lw, label="HF-SKE2", path_effects=path_effects)
Ar40_ax.plot(Ar40_NNLO_Barbieri_Qvals, Ar40_NNLO_Barbieri_FFvals, color=colors['Barbieri'], ls=(0,(2,2)), lw=lw, label="NNLOsat (Barbieri)", path_effects=path_effects)
Ar40_ax.plot(Ar40_NNLO_Payne_Qvals, Ar40_NNLO_Payne_FFvals, color=colors['Payne'], ls='-', lw=lw, label="NNLOsat (Payne)", path_effects=path_effects)
#Ar40_ax.plot(Qvals, Ar40FBSE_FF_vals, color='#000000', ls='-', lw=2.5, label="Fourier-Bessel FF")
Ar40_ax.tick_params(which='both', right=False, top=False)
Ar40_ax.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
Ar40_ax.yaxis.grid(True, ls='-', which='both', color='gray', alpha=0.35)
Ar40_ax.set_xlabel(r'Momentum Transfer $Q$ [GeV]')
Ar40_ax.set_ylabel(r'Form Factor $|F(Q^2)|$')
Ar40_ax.legend(fontsize=20, loc='upper right', frameon=False, framealpha=1, edgecolor='k', fancybox=True)
Ar40_ax.set_title(r'$^{40}$Ar Nuclear Form Factor')
Ar40_ax.set_yscale('log')
Ar40_ax.set_ylim([1e-4, 1])
Ar40_ax.set_xlim([0, Qmax])
Ar40_fig.savefig('../Plots/Figures_PDF/rho_FF_u1/40ArgonFF.pdf', transparent=False, bbox_inches='tight')
Ar40_fig.savefig('../Plots/Figures_PNG/rho_FF_u1/40ArgonFF.png', transparent=False, bbox_inches='tight')


W184_ax.plot(Qvals, W184WS_FF_values, color='blue', label='Woods-Saxon', lw=2, path_effects=path_effects)
W184_ax.plot(Qvals, W184HELMv1_FF_values, color='red', label='Helm [v1]', lw=2, path_effects=path_effects)
W184_ax.plot(Qvals, W184HELMv2_FF_values, color='green', label='Helm [v2]', lw=2, path_effects=path_effects)
W184_ax.plot(Qvals, W184KN_FF_values, color='purple', label='KN', lw=2, path_effects=path_effects)
W184_ax.plot(Qvals, W184adKN_FF_values, color='orange', label='(ad.) KN', lw=2, path_effects=path_effects)
W184_ax.set_yscale('log')
W184_ax.set_ylim([1e-4, 1])
W184_ax.set_xlim([0, Qmax])
W184_ax.xaxis.grid(True, which='major', ls='-', alpha=0.35)
W184_ax.yaxis.grid(True, which='both', ls='-', alpha=0.35)
W184_ax.legend(fontsize=20, frameon=False, framealpha=1, edgecolor='k', fancybox=True)
W184_ax.set_xlabel(r'Momentum Transfer $Q$ [Gev]')
W184_ax.set_ylabel(r'Form Factor $|F(Q^2)|$')
W184_ax.set_title(r'Tungsten $^{184}$W Nuclear Form Factor')
W184_fig.savefig('../Plots/Figures_PDF/rho_FF_u1/184TungstenFF.pdf', transparent=False, bbox_inches='tight')
W184_fig.savefig('../Plots/Figures_PNG/rho_FF_u1/184TungstenFF.png', transparent=False, bbox_inches='tight')


Proton_ax.plot(Qvals_incoh, PROTONGe_Alt_values, label=r'$G^p_E$ [Alt]')
Proton_ax.plot(Qvals_incoh, PROTONGm_Alt_values, label=r'$G^p_M$ [Alt]')
Proton_ax.plot(Qvals_incoh, PROTONGe_Kelly_values, label=r'$G^p_E$ [Kelly]')
Proton_ax.plot(Qvals_incoh, PROTONGm_Kelly_values, label=r'$G^p_M$ [Kelly]')
Proton_ax.set_yscale('log')
Proton_ax.set_xlim([0, Qmax_incoh])
Proton_ax.legend(fontsize=20, loc='upper right', frameon=False, framealpha=1, edgecolor='k', fancybox=True)
Proton_ax.tick_params(which='both', right=True, top=True)
Proton_ax.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
Proton_ax.yaxis.grid(True, ls='-', which='both', color='gray', alpha=0.35)
Proton_ax.set_xlabel(r'Momentum Transfer $Q$ [GeV]')
Proton_ax.set_ylabel(r'$|G^p(Q^2)|$')
Proton_ax.set_title(r'Proton Form Factors')
Proton_fig.savefig('../Plots/Figures_PDF/rho_FF_u1/protonFF.pdf', transparent=False, bbox_inches='tight')
Proton_fig.savefig('../Plots/Figures_PNG/rho_FF_u1/protonFF.png', transparent=False, bbox_inches='tight')


Neutron_ax.plot(Qvals_incoh, NEUTRONGe_Alt_values, label=r'$G^n_E$ [Alt]')
Neutron_ax.plot(Qvals_incoh, NEUTRONGm_Alt_values, label=r'$G^n_M$ [Alt]')
Neutron_ax.plot(Qvals_incoh, NEUTRONGe_Kelly_values, label=r'$G^n_E$ [Kelly]')
Neutron_ax.plot(Qvals_incoh, NEUTRONGm_Kelly_values, label=r'$G^n_M$ [Kelly]')
Neutron_ax.set_yscale('log')
Neutron_ax.set_xlim([0, Qmax_incoh])
Neutron_ax.legend(fontsize=20, loc='upper right', frameon=False, framealpha=1, edgecolor='k', fancybox=True)
Neutron_ax.tick_params(which='both', right=True, top=True)
Neutron_ax.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
Neutron_ax.yaxis.grid(True, ls='-', which='both', color='gray', alpha=0.35)
Neutron_ax.set_xlabel(r'Momentum Transfer $Q$ [GeV]')
Neutron_ax.set_ylabel(r'$|G^n(Q^2)|$')
Neutron_ax.set_title(r'Neutron Form Factors')
Neutron_fig.savefig('../Plots/Figures_PDF/rho_FF_u1/neutronFF.pdf', transparent=False, bbox_inches='tight')
Neutron_fig.savefig('../Plots/Figures_PNG/rho_FF_u1/neutronFF.png', transparent=False, bbox_inches='tight')

plt.show()
