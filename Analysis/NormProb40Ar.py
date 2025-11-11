import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patheffects as pe
from scipy.integrate import simpson
import csv

plt.style.use('sty.mplstyle')
font = 35

def StepFunction(xmin, xmax, y):
    x = np.concatenate((xmin[0:1], xmax))
    y = np.concatenate((y[0:1], y))
    return x, y

#=================================================
#======   BNB MiniBooNE NuMode (** OK **)   ======
#=================================================
Elow_MiniBooNE_nuMode = np.loadtxt("../csv/Fluxes/MiniBooNE_Fluxes_NuMode.dat.txt", usecols=0, comments="#") 				# Lower bin edge (GeV)
Ehigh_MiniBooNE_nuMode = np.loadtxt("../csv/Fluxes/MiniBooNE_Fluxes_NuMode.dat.txt", usecols=1, comments="#") 				# Upper bin edge (GeV)
numu_flux_MiniBooNE_nuMode = 2e5*np.loadtxt("../csv/Fluxes/MiniBooNE_Fluxes_NuMode.dat.txt", usecols=2, comments="#") 		# numu/m^2/POT/GeV
numubar_flux_MiniBooNE_nuMode = 2e5*np.loadtxt("../csv/Fluxes/MiniBooNE_Fluxes_NuMode.dat.txt", usecols=3, comments="#") 	# numubar/m^2/POT/GeV
#=================================================
Energy_MiniBooNE_nuMode = (Elow_MiniBooNE_nuMode + Ehigh_MiniBooNE_nuMode)/2
#=================================================
nu_total_flux_MiniBooNE_nuMode = simpson(numu_flux_MiniBooNE_nuMode, Energy_MiniBooNE_nuMode)
antinu_total_flux_MiniBooNE_nuMode = simpson(numubar_flux_MiniBooNE_nuMode, Energy_MiniBooNE_nuMode)
nu_normProb_MiniBooNE_nuMode = np.zeros_like(numu_flux_MiniBooNE_nuMode)
antinu_normProb_MiniBooNE_nuMode = np.zeros_like(numubar_flux_MiniBooNE_nuMode)
#=================================================
for i in range(len(numu_flux_MiniBooNE_nuMode)):
    bin_width = Ehigh_MiniBooNE_nuMode[i] - Elow_MiniBooNE_nuMode[i]
    nu_normProb_MiniBooNE_nuMode[i] = numu_flux_MiniBooNE_nuMode[i]*bin_width/nu_total_flux_MiniBooNE_nuMode
#=================================================
for i in range(len(numubar_flux_MiniBooNE_nuMode)):
    bin_width = Ehigh_MiniBooNE_nuMode[i] - Elow_MiniBooNE_nuMode[i]
    antinu_normProb_MiniBooNE_nuMode[i] = numubar_flux_MiniBooNE_nuMode[i]*bin_width/antinu_total_flux_MiniBooNE_nuMode
#=================================================
#=========   Writing to CSV (** OK **)   =========
#=================================================
with open("../csv/Fluxes/NormProb/nu_NormProb_BNB_MiniBooNE_nuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Neutrino Flux (NuMode): {nu_total_flux_MiniBooNE_nuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_MiniBooNE_nuMode, Ehigh_MiniBooNE_nuMode, nu_normProb_MiniBooNE_nuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

with open("../csv/Fluxes/NormProb/antinu_NormProb_BNB_MiniBooNE_nuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Antineutrino Flux (NuMode): {antinu_total_flux_MiniBooNE_nuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_MiniBooNE_nuMode, Ehigh_MiniBooNE_nuMode, antinu_normProb_MiniBooNE_nuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")


#=================================================
#====   BNB MiniBooNE AntiNuMode (** OK **)   ====
#=================================================
Elow_MiniBooNE_antinuMode = np.loadtxt("../csv/Fluxes/MiniBooNE_Fluxes_AnuMode.dat.txt", usecols=0, comments="#") 				# Lower bin edge (GeV)
Ehigh_MiniBooNE_antinuMode = np.loadtxt("../csv/Fluxes/MiniBooNE_Fluxes_AnuMode.dat.txt", usecols=1, comments="#") 				# Upper bin edge (GeV)
numu_flux_MiniBooNE_antinuMode = 2e5*np.loadtxt("../csv/Fluxes/MiniBooNE_Fluxes_AnuMode.dat.txt", usecols=2, comments="#") 		# numu/m^2/POT/GeV
numubar_flux_MiniBooNE_antinuMode = 2e5*np.loadtxt("../csv/Fluxes/MiniBooNE_Fluxes_AnuMode.dat.txt", usecols=3, comments="#") 	# numubar/m^2/POT/GeV
#=================================================
Energy_MiniBooNE_antinuMode = (Elow_MiniBooNE_antinuMode + Ehigh_MiniBooNE_antinuMode)/2
#=================================================
nu_total_flux_MiniBooNE_antinuMode = simpson(numu_flux_MiniBooNE_antinuMode, Energy_MiniBooNE_antinuMode)
antinu_total_flux_MiniBooNE_antinuMode = simpson(numubar_flux_MiniBooNE_antinuMode, Energy_MiniBooNE_antinuMode)
nu_normProb_MiniBooNE_antinuMode = np.zeros_like(numu_flux_MiniBooNE_antinuMode)
antinu_normProb_MiniBooNE_antinuMode = np.zeros_like(numubar_flux_MiniBooNE_antinuMode)
#=================================================
for i in range(len(numu_flux_MiniBooNE_antinuMode)):
    bin_width = Ehigh_MiniBooNE_antinuMode[i] - Elow_MiniBooNE_antinuMode[i]
    nu_normProb_MiniBooNE_antinuMode[i] = numu_flux_MiniBooNE_antinuMode[i]*bin_width/nu_total_flux_MiniBooNE_antinuMode
#=================================================
for i in range(len(numubar_flux_MiniBooNE_antinuMode)):
    bin_width = Ehigh_MiniBooNE_antinuMode[i] - Elow_MiniBooNE_antinuMode[i]
    antinu_normProb_MiniBooNE_antinuMode[i] = numubar_flux_MiniBooNE_antinuMode[i]*bin_width/antinu_total_flux_MiniBooNE_antinuMode
#=================================================
#=========   Writing to CSV (** OK **)   =========
#=================================================
with open("../csv/Fluxes/NormProb/nu_NormProb_BNB_MiniBooNE_antinuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Neutrino Flux (AntiNuMode): {nu_total_flux_MiniBooNE_antinuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_MiniBooNE_antinuMode, Ehigh_MiniBooNE_antinuMode, nu_normProb_MiniBooNE_antinuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

with open("../csv/Fluxes/NormProb/antinu_NormProb_BNB_MiniBooNE_antinuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Antineutrino Flux (AntiNuMode): {antinu_total_flux_MiniBooNE_antinuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_MiniBooNE_antinuMode, Ehigh_MiniBooNE_antinuMode, antinu_normProb_MiniBooNE_antinuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")


#=================================================
#=========   BNB SBND NuMode (** OK **)   ========
#=================================================
Elow_SBND_nuMode = np.loadtxt("../csv/Fluxes/SBND_Fluxes_NeutrinoMode.dat", usecols=0, comments="#") 				# Lower bin edge (GeV)
Ehigh_SBND_nuMode = np.loadtxt("../csv/Fluxes/SBND_Fluxes_NeutrinoMode.dat", usecols=1, comments="#") 				# Upper bin edge (GeV)
numu_flux_SBND_nuMode = 2e-5*np.loadtxt("../csv/Fluxes/SBND_Fluxes_NeutrinoMode.dat", usecols=2, comments="#") 		# numu/m^2/POT/GeV
numubar_flux_SBND_nuMode = 2e-5*np.loadtxt("../csv/Fluxes/SBND_Fluxes_NeutrinoMode.dat", usecols=3, comments="#") 	# numubar/m^2/POT/GeV
#=================================================
Energy_SBND_nuMode = (Elow_SBND_nuMode + Ehigh_SBND_nuMode)/2
#=================================================
nu_total_flux_SBND_nuMode = simpson(numu_flux_SBND_nuMode, Energy_SBND_nuMode)
antinu_total_flux_SBND_nuMode = simpson(numubar_flux_SBND_nuMode, Energy_SBND_nuMode)
nu_normProb_SBND_nuMode = np.zeros_like(numu_flux_SBND_nuMode)
antinu_normProb_SBND_nuMode = np.zeros_like(numubar_flux_SBND_nuMode)
#=================================================
for i in range(len(numu_flux_SBND_nuMode)):
    bin_width = Ehigh_SBND_nuMode[i] - Elow_SBND_nuMode[i]
    nu_normProb_SBND_nuMode[i] = numu_flux_SBND_nuMode[i]*bin_width/nu_total_flux_SBND_nuMode
#=================================================
for i in range(len(numubar_flux_SBND_nuMode)):
    bin_width = Ehigh_SBND_nuMode[i] - Elow_SBND_nuMode[i]
    antinu_normProb_SBND_nuMode[i] = numubar_flux_SBND_nuMode[i]*bin_width/antinu_total_flux_SBND_nuMode
#=================================================
#=========   Writing to CSV (** OK **)   =========
#=================================================
with open("../csv/Fluxes/NormProb/nu_NormProb_BNB_SBND_nuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Neutrino Flux (NuMode): {nu_total_flux_SBND_nuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_SBND_nuMode, Ehigh_SBND_nuMode, nu_normProb_SBND_nuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

with open("../csv/Fluxes/NormProb/antinu_NormProb_BNB_SBND_nuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Antineutrino Flux (NuMode): {antinu_total_flux_SBND_nuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_SBND_nuMode, Ehigh_SBND_nuMode, antinu_normProb_SBND_nuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

#=================================================
#=======   BNB SBND AntiNuMode (** OK **)   ======
#=================================================
Elow_SBND_antinuMode = np.loadtxt("../csv/Fluxes/SBND_Fluxes_AnuMode.dat", usecols=0, comments="#") 				# Lower bin edge (GeV)
Ehigh_SBND_antinuMode = np.loadtxt("../csv/Fluxes/SBND_Fluxes_AnuMode.dat", usecols=1, comments="#") 				# Upper bin edge (GeV)
numu_flux_SBND_antinuMode = 2e-5*np.loadtxt("../csv/Fluxes/SBND_Fluxes_AnuMode.dat", usecols=2, comments="#") 		# numu/m^2/POT/GeV
numubar_flux_SBND_antinuMode = 2e-5*np.loadtxt("../csv/Fluxes/SBND_Fluxes_AnuMode.dat", usecols=3, comments="#") 	# numubar/m^2/POT/GeV
#=================================================
Energy_SBND_antinuMode = (Elow_SBND_antinuMode + Ehigh_SBND_antinuMode)/2
#=================================================
nu_total_flux_SBND_antinuMode = simpson(numu_flux_SBND_antinuMode, Energy_SBND_antinuMode)
antinu_total_flux_SBND_antinuMode = simpson(numubar_flux_SBND_antinuMode, Energy_SBND_antinuMode)
nu_normProb_SBND_antinuMode = np.zeros_like(numu_flux_SBND_antinuMode)
antinu_normProb_SBND_antinuMode = np.zeros_like(numubar_flux_SBND_antinuMode)
#=================================================
for i in range(len(numu_flux_SBND_antinuMode)):
    bin_width = Ehigh_SBND_antinuMode[i] - Elow_SBND_antinuMode[i]
    nu_normProb_SBND_antinuMode[i] = numu_flux_SBND_antinuMode[i]*bin_width/nu_total_flux_SBND_antinuMode
#=================================================
for i in range(len(numubar_flux_SBND_antinuMode)):
    bin_width = Ehigh_SBND_antinuMode[i] - Elow_SBND_antinuMode[i]
    antinu_normProb_SBND_antinuMode[i] = numubar_flux_SBND_antinuMode[i]*bin_width/antinu_total_flux_SBND_antinuMode
#=================================================
#=========   Writing to CSV (** OK **)   =========
#=================================================
with open("../csv/Fluxes/NormProb/nu_NormProb_BNB_SBND_antinuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Neutrino Flux (AntiNuMode): {nu_total_flux_SBND_antinuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_SBND_antinuMode, Ehigh_SBND_antinuMode, nu_normProb_SBND_antinuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

with open("../csv/Fluxes/NormProb/antinu_NormProb_BNB_SBND_antinuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Antineutrino Flux (AntiNuMode): {antinu_total_flux_SBND_antinuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_SBND_antinuMode, Ehigh_SBND_antinuMode, antinu_normProb_SBND_antinuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

#=================================================
#======   Std Flux DUNE NuMode (** OK **)   ======
#=================================================
DUNEStd_nuMode_file = "../csv/Fluxes/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_neutrino_LBNEND_globes_flux.txt"
Energy_DUNEStd_nuMode = np.loadtxt(DUNEStd_nuMode_file, usecols=0, comments="#")
numu_flux_DUNEStd_nuMode = np.loadtxt(DUNEStd_nuMode_file, usecols=2, comments="#")
numubar_flux_DUNEStd_nuMode = np.loadtxt(DUNEStd_nuMode_file, usecols=5, comments="#")
#=================================================
Elow_DUNEStd_nuMode = Energy_DUNEStd_nuMode - 0.25/2	# bin width = 0.25 GeV
Ehigh_DUNEStd_nuMode = Energy_DUNEStd_nuMode + 0.25/2	# bin width = 0.25 GeV
#=================================================
nu_total_flux_DUNEStd_nuMode = simpson(numu_flux_DUNEStd_nuMode, Energy_DUNEStd_nuMode)
antinu_total_flux_DUNEStd_nuMode = simpson(numubar_flux_DUNEStd_nuMode, Energy_DUNEStd_nuMode)
nu_normProb_DUNEStd_nuMode = np.zeros_like(numu_flux_DUNEStd_nuMode)
antinu_normProb_DUNEStd_nuMode = np.zeros_like(numubar_flux_DUNEStd_nuMode)
#=================================================
for i in range(len(numu_flux_DUNEStd_nuMode)):
    bin_width = Ehigh_DUNEStd_nuMode[i] - Elow_DUNEStd_nuMode[i]
    nu_normProb_DUNEStd_nuMode[i] = numu_flux_DUNEStd_nuMode[i]*bin_width/nu_total_flux_DUNEStd_nuMode
    
for i in range(len(numubar_flux_DUNEStd_nuMode)): 
    bin_width = Ehigh_DUNEStd_nuMode[i] - Elow_DUNEStd_nuMode[i] 
    antinu_normProb_DUNEStd_nuMode[i] = numubar_flux_DUNEStd_nuMode[i]*bin_width/antinu_total_flux_DUNEStd_nuMode
#=================================================
#=========   Writing to CSV (** OK **)   =========
#=================================================
with open("../csv/Fluxes/NormProb/nu_NormProb_DUNEStd_nuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Neutrino Flux (NuMode): {nu_total_flux_DUNEStd_nuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_DUNEStd_nuMode, Ehigh_DUNEStd_nuMode, nu_normProb_DUNEStd_nuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

with open("../csv/Fluxes/NormProb/antinu_NormProb_DUNEStd_nuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Antineutrino Flux (NuMode): {antinu_total_flux_DUNEStd_nuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_DUNEStd_nuMode, Ehigh_DUNEStd_nuMode, antinu_normProb_DUNEStd_nuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

#=================================================
#====   Std Flux DUNE AntiNuMode (** OK **)   ====
#=================================================
DUNEStd_AntinuMode_file = "../csv/Fluxes/histos_g4lbne_v3r5p4_QGSP_BERT_OptimizedEngineeredNov2017_antineutrino_LBNEND_globes_flux.txt"
Energy_DUNEStd_AntinuMode = np.loadtxt(DUNEStd_AntinuMode_file, usecols=0, comments="#")
numu_flux_DUNEStd_AntinuMode = np.loadtxt(DUNEStd_AntinuMode_file, usecols=2, comments="#")
numubar_flux_DUNEStd_AntinuMode = np.loadtxt(DUNEStd_AntinuMode_file, usecols=5, comments="#")
#=================================================
Elow_DUNEStd_AntinuMode = Energy_DUNEStd_AntinuMode - 0.25/2	# bin width = 0.25 GeV
Ehigh_DUNEStd_AntinuMode = Energy_DUNEStd_AntinuMode + 0.25/2	# bin width = 0.25 GeV
#=================================================
nu_total_flux_DUNEStd_AntinuMode = simpson(numu_flux_DUNEStd_AntinuMode, Energy_DUNEStd_AntinuMode)
antinu_total_flux_DUNEStd_AntinuMode = simpson(numubar_flux_DUNEStd_AntinuMode, Energy_DUNEStd_AntinuMode)
nu_normProb_DUNEStd_AntinuMode = np.zeros_like(numu_flux_DUNEStd_AntinuMode)
antinu_normProb_DUNEStd_AntinuMode = np.zeros_like(numubar_flux_DUNEStd_AntinuMode)
#=================================================
for i in range(len(numu_flux_DUNEStd_AntinuMode)):
    bin_width = Ehigh_DUNEStd_AntinuMode[i] - Elow_DUNEStd_AntinuMode[i]
    nu_normProb_DUNEStd_AntinuMode[i] = numu_flux_DUNEStd_AntinuMode[i]*bin_width/nu_total_flux_DUNEStd_AntinuMode
    
for i in range(len(numubar_flux_DUNEStd_AntinuMode)): 
    bin_width = Ehigh_DUNEStd_AntinuMode[i] - Elow_DUNEStd_AntinuMode[i] 
    antinu_normProb_DUNEStd_AntinuMode[i] = numubar_flux_DUNEStd_AntinuMode[i]*bin_width/antinu_total_flux_DUNEStd_AntinuMode
#=================================================
#=========   Writing to CSV (** OK **)   =========
#=================================================
#with open("../csv/Fluxes/NormProb/nu_NormProb_DUNEStd_nuMode.txt", mode='w', newline="") as file:
#    file.write("########################################################################################## \n")
#    file.write(f"#####\t Integrated Neutrino Flux (NuMode): {nu_total_flux_DUNEStd_nuMode} m^-2 POT^-1 \n")
#    file.write("########################################################################################## \n")
#    for e_low, e_high, p_norm in (zip(Elow_DUNEStd_nuMode, Ehigh_DUNEStd_nuMode, nu_normProb_DUNEStd_nuMode)):
#        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

#with open("../csv/Fluxes/NormProb/antinu_NormProb_DUNEStd_nuMode.txt", mode='w', newline="") as file:
#    file.write("########################################################################################## \n")
#    file.write(f"#####\t Integrated Antineutrino Flux (NuMode): {antinu_total_flux_DUNEStd_nuMode} m^-2 POT^-1 \n")
#    file.write("########################################################################################## \n")
#    for e_low, e_high, p_norm in (zip(Elow_DUNEStd_nuMode, Ehigh_DUNEStd_nuMode, antinu_normProb_DUNEStd_nuMode)):
#        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

#=================================================
#====   tau-Opt Flux DUNE NuMode (** OK **)   ====
#=================================================
Energy_DUNETauOpt_nuMode = np.loadtxt("../csv/Fluxes/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_neutrino_LBNEND_globes_flux.txt", usecols=0, comments="#") 
numu_flux_DUNETauOpt_nuMode = np.loadtxt("../csv/Fluxes/histos_g4lbne_v3r5p4_QGSP_BERT_TauOptimized_neutrino_LBNEND_globes_flux.txt", usecols=2, comments="#")
#=================================================
Elow_DUNETauOpt_nuMode = Energy_DUNETauOpt_nuMode - 0.25/2		# bin width = 0.25 GeV
Ehigh_DUNETauOpt_nuMode = Energy_DUNETauOpt_nuMode + 0.25/2		# bin width = 0.25 GeV
#=================================================
nu_total_flux_DUNETauOpt_nuMode = simpson(numu_flux_DUNETauOpt_nuMode, Energy_DUNETauOpt_nuMode)
nu_normProb_DUNETauOpt_nuMode = np.zeros_like(numu_flux_DUNETauOpt_nuMode)
#=================================================
for i in range(len(numu_flux_DUNETauOpt_nuMode)):
    bin_width = Ehigh_DUNETauOpt_nuMode[i] - Elow_DUNETauOpt_nuMode[i]
    nu_normProb_DUNETauOpt_nuMode[i] = numu_flux_DUNETauOpt_nuMode[i]*bin_width/nu_total_flux_DUNETauOpt_nuMode
#=================================================
#=========   Writing to CSV (** OK **)   =========
#=================================================
with open("../csv/Fluxes/NormProb/nu_NormProb_DUNETauOpt_nuMode.txt", mode='w', newline="") as file:
    file.write("########################################################################################## \n")
    file.write(f"#####\t Integrated Neutrino Flux (NuMode): {nu_total_flux_DUNETauOpt_nuMode} m^-2 POT^-1 \n")
    file.write("########################################################################################## \n")
    for e_low, e_high, p_norm in (zip(Elow_DUNETauOpt_nuMode, Ehigh_DUNETauOpt_nuMode, nu_normProb_DUNETauOpt_nuMode)):
        file.write(f"\t {e_low:.4f}, \t {e_high:.4f}, \t {p_norm:.22f}, \n")

#=================================================
#=====   Loading Flux at Icarus (** OK **)   =====
#=================================================
Elow_ICARUS_nuMode = np.loadtxt("../csv/Fluxes/ICARUS_NuMode.txt", usecols=0, comments="#") 		# Lower bin edge (GeV)
Ehigh_ICARUS_nuMode = np.loadtxt("../csv/Fluxes/ICARUS_NuMode.txt", usecols=1, comments="#") 		# Upper bin edge (GeV)
numu_flux_ICARUS_nuMode = np.loadtxt("../csv/Fluxes/ICARUS_NuMode.txt", usecols=2, comments="#") 	# numu/m^2/POT/GeV
numubar_flux_ICARUS_nuMode = np.loadtxt("../csv/Fluxes/ICARUS_NuMode.txt", usecols=3, comments="#") # numubar/m^2/POT/GeV

Elow_ICARUS_antinuMode = np.loadtxt("../csv/Fluxes/ICARUS_AnuMode.txt", usecols=0, comments="#") 			# Lower bin edge (GeV)
Ehigh_ICARUS_antinuMode = np.loadtxt("../csv/Fluxes/ICARUS_AnuMode.txt", usecols=1, comments="#") 			# Upper bin edge (GeV)
numu_flux_ICARUS_antinuMode = np.loadtxt("../csv/Fluxes/ICARUS_AnuMode.txt", usecols=2, comments="#") 		# numu/m^2/POT/GeV
numubar_flux_ICARUS_antinuMode = np.loadtxt("../csv/Fluxes/ICARUS_AnuMode.txt", usecols=3, comments="#")	# numubar/m^2/POT/GeV

#=================================================
#=======   Plotting the Fluxes (** OK **)   ======
#=================================================
fig = plt.figure(figsize=(30,28))
gs = fig.add_gridspec(2,2,hspace=0.2, wspace=0.22)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[1, 0])
ax4 = fig.add_subplot(gs[1, 1])

#fig = plt.figure()
#ax = fig.add_subplot()

ax1.step(StepFunction(Elow_MiniBooNE_nuMode,Ehigh_MiniBooNE_nuMode,numu_flux_MiniBooNE_nuMode)[0],
		StepFunction(Elow_MiniBooNE_nuMode,Ehigh_MiniBooNE_nuMode,numu_flux_MiniBooNE_nuMode)[1],
		ls='-', color='#58094F', label=r'$\nu_\mu$ [$\nu$-Mode]', path_effects=[pe.Stroke(linewidth=4,foreground='k'), pe.Normal()])
ax1.step(StepFunction(Elow_MiniBooNE_nuMode,Ehigh_MiniBooNE_nuMode,numubar_flux_MiniBooNE_nuMode)[0],
		StepFunction(Elow_MiniBooNE_nuMode,Ehigh_MiniBooNE_nuMode,numubar_flux_MiniBooNE_nuMode)[1],
		ls='--', color='#58094F', label=r'$\overline{\nu}_\mu$ [$\nu$-Mode]', path_effects=[pe.Stroke(linewidth=4,foreground='k'), pe.Normal()])
ax1.step(StepFunction(Elow_MiniBooNE_antinuMode,Ehigh_MiniBooNE_antinuMode,numu_flux_MiniBooNE_antinuMode)[0],
		StepFunction(Elow_MiniBooNE_antinuMode,Ehigh_MiniBooNE_antinuMode,numu_flux_MiniBooNE_antinuMode)[1],
		ls='-', color='#295E11', label=r'$\nu_\mu$ [$\overline{\nu}$-Mode]', path_effects=[pe.Stroke(linewidth=4,foreground='k'), pe.Normal()])
ax1.step(StepFunction(Elow_MiniBooNE_antinuMode,Ehigh_MiniBooNE_antinuMode,numubar_flux_MiniBooNE_antinuMode)[0],
		StepFunction(Elow_MiniBooNE_antinuMode,Ehigh_MiniBooNE_antinuMode,numubar_flux_MiniBooNE_antinuMode)[1],
		ls='--', color='#295E11', label=r'$\overline{\nu}_\mu$ [$\overline{\nu}$-Mode]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])

ax1.set_xlim([0,7.5])
ax1.set_yscale('log')
ax1.set_ylim(top=1e-5)
ax1.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
ax1.yaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
ax1.tick_params(which='both', right=True, top=True)
ax1.set_title(r'{\bf Flux at MiniBooNE in} $\nu / \overline{\nu}$ {\bf Mode}', pad=16, fontsize=font)
ax1.set_xlabel(r'{\bf Neutrino Energy} $E_\nu$ [GeV]', fontsize=font)
ax1.set_ylabel(r'$\frac{\mathrm{d} \Phi}{\mathrm{d} E_\nu}$  [$\nu$/m$^2$/POT/GeV]', fontsize=font)
ax1.legend(fontsize=0.7*font, loc='upper right', frameon=False, framealpha=1, edgecolor='k', fancybox=True)
#plt.savefig('../Plots/Figures_PDF/MiniBooNE/MiniBooNE_Flux.pdf', transparent=False, bbox_inches='tight')
#plt.savefig('../Plots/Figures_PNG/MiniBooNE/MiniBooNE_Flux.png', transparent=False, bbox_inches='tight')

#=================================================
#=================================================
#fig = plt.figure()
#ax = fig.add_subplot()

ax2.step(StepFunction(Elow_SBND_nuMode,Ehigh_SBND_nuMode,numu_flux_SBND_nuMode)[0],
		StepFunction(Elow_SBND_nuMode,Ehigh_SBND_nuMode,numu_flux_SBND_nuMode)[1],
		ls='-', color='#58094F', label=r'$\nu_\mu$ [$\nu$-Mode]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.step(StepFunction(Elow_SBND_nuMode,Ehigh_SBND_nuMode,numubar_flux_SBND_nuMode)[0],
		StepFunction(Elow_SBND_nuMode,Ehigh_SBND_nuMode,numubar_flux_SBND_nuMode)[1],
		ls='--', color='#58094F', label=r'$\overline{\nu}_\mu$ [$\nu$-Mode]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.step(StepFunction(Elow_SBND_antinuMode,Ehigh_SBND_antinuMode,numu_flux_SBND_antinuMode)[0],
		StepFunction(Elow_SBND_antinuMode,Ehigh_SBND_antinuMode,numu_flux_SBND_antinuMode)[1],
		ls='-', color='#295E11', label=r'$\nu_\mu$ [$\overline{\nu}$-Mode]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])
ax2.step(StepFunction(Elow_SBND_antinuMode,Ehigh_SBND_antinuMode,numubar_flux_SBND_antinuMode)[0],
		StepFunction(Elow_SBND_antinuMode,Ehigh_SBND_antinuMode,numubar_flux_SBND_antinuMode)[1],
		ls='--', color='#295E11', label=r'$\overline{\nu}_\mu$ [$\overline{\nu}$-Mode]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])

ax2.set_xlim([0,7.5])
ax2.set_yscale('log')
ax2.set_ylim(top=3e-4)
ax2.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
ax2.yaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
ax2.tick_params(which='both', right=True, top=True)
ax2.set_title(r'{\bf Flux at SBND in} $\nu / \overline{\nu}$ {\bf Mode}', pad=16, fontsize=font)
ax2.set_xlabel(r'{\bf Neutrino Energy} $E_\nu$ [GeV]', fontsize=font)
ax2.set_ylabel(r'$\frac{\mathrm{d} \Phi}{\mathrm{d} E_\nu}$  [$\nu$/m$^2$/POT/GeV]', fontsize=font)
ax2.legend(fontsize=0.7*font, loc='upper right', frameon=False, framealpha=1, edgecolor='k', fancybox=True)
#plt.savefig('../Plots/Figures_PDF/SBND/SBND_Flux.pdf', transparent=False, bbox_inches='tight')
#plt.savefig('../Plots/Figures_PNG/SBND/SBND_Flux.png', transparent=False, bbox_inches='tight')

#=================================================
#=================================================
#fig = plt.figure()
#ax = fig.add_subplot()

ax3.step(StepFunction(Elow_DUNEStd_nuMode, Ehigh_DUNEStd_nuMode, numu_flux_DUNEStd_nuMode)[0],
		StepFunction(Elow_DUNEStd_nuMode, Ehigh_DUNEStd_nuMode, numu_flux_DUNEStd_nuMode)[1],
		color='#58094F', lw=4, ls='-', label=r'$\nu_\mu$ Std. [$\nu$-Mode]', path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
ax3.step(StepFunction(Elow_DUNEStd_nuMode, Ehigh_DUNEStd_nuMode, numubar_flux_DUNEStd_nuMode)[0], # ???????????
		StepFunction(Elow_DUNEStd_nuMode, Ehigh_DUNEStd_nuMode, numubar_flux_DUNEStd_nuMode)[1], # ???????????????
		color='#58094F', lw=4, ls='--', label=r'$\overline{\nu}_\mu$ Std. [$\nu$-Mode]', path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])
#ax3.step(StepFunction(Elow_DUNETauOpt_nuMode, Ehigh_DUNETauOpt_nuMode,numu_flux_DUNETauOpt_nuMode)[0],
		#StepFunction(Elow_DUNETauOpt_nuMode, Ehigh_DUNETauOpt_nuMode,numu_flux_DUNETauOpt_nuMode)[1],
		#color='#58094F', lw=4, ls='--', label=r'$\nu_\mu \, \tau$-Opt.', path_effects=[pe.Stroke(linewidth=5, foreground='k'), pe.Normal()])

ax3.set_yscale('log')
ax3.set_ylim(bottom=4e-8, top=1e-3)
ax3.set_xscale('log')
ax3.set_xlim(left=1e-1, right=1e2)
ax3.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
ax3.yaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
ax3.legend(fontsize=0.7*font, loc='upper right', frameon=False, framealpha=1, edgecolor='k', fancybox=True)
ax3.tick_params(which='both', right=True, top=True)
ax3.set_xlabel(r'{\bf Neutrino Energy} $E_\nu$ [GeV]', fontsize=font)
ax3.set_ylabel(r'$\frac{\mathrm{d} \Phi}{\mathrm{d} E_\nu}$  [$\nu$/m$^2$/POT/GeV]', fontsize=font)
ax3.set_title(r'{\bf Flux at DUNE in} $\nu$ {\bf Mode}', pad=16, fontsize=font)
#plt.savefig('../Plots/Figures_PDF/DUNE/FluxatDUNE.pdf', transparent=False, bbox_inches='tight')
#plt.savefig('../Plots/Figures_PNG/DUNE/FluxatDUNE.png', transparent=False, bbox_inches='tight')

#=================================================
#=================================================
ax4.step(StepFunction(Elow_ICARUS_nuMode,Ehigh_ICARUS_nuMode,numu_flux_ICARUS_nuMode)[0],
		StepFunction(Elow_ICARUS_nuMode,Ehigh_ICARUS_nuMode,numu_flux_ICARUS_nuMode)[1],
		ls='-', color='#58094F', label=r'$\nu_\mu$ [$\nu$-Mode]', path_effects=[pe.Stroke(linewidth=4,foreground='k'), pe.Normal()])
ax4.step(StepFunction(Elow_ICARUS_nuMode,Ehigh_ICARUS_nuMode,numubar_flux_ICARUS_nuMode)[0],
		StepFunction(Elow_ICARUS_nuMode,Ehigh_ICARUS_nuMode,numubar_flux_ICARUS_nuMode)[1],
		ls='--', color='#58094F', label=r'$\overline{\nu}_\mu$ [$\nu$-Mode]', path_effects=[pe.Stroke(linewidth=4,foreground='k'), pe.Normal()])
ax4.step(StepFunction(Elow_ICARUS_antinuMode,Ehigh_ICARUS_antinuMode,numu_flux_ICARUS_antinuMode)[0],
		StepFunction(Elow_ICARUS_antinuMode,Ehigh_ICARUS_antinuMode,numu_flux_ICARUS_antinuMode)[1],
		ls='-', color='#295E11', label=r'$\nu_\mu$ [$\overline{\nu}$-Mode]', path_effects=[pe.Stroke(linewidth=4,foreground='k'), pe.Normal()])
ax4.step(StepFunction(Elow_ICARUS_antinuMode,Ehigh_ICARUS_antinuMode,numubar_flux_ICARUS_antinuMode)[0],
		StepFunction(Elow_ICARUS_antinuMode,Ehigh_ICARUS_antinuMode,numubar_flux_ICARUS_antinuMode)[1],
		ls='--', color='#295E11', label=r'$\overline{\nu}_\mu$ [$\overline{\nu}$-Mode]', path_effects=[pe.Stroke(linewidth=4, foreground='k'), pe.Normal()])

ax4.set_xlim([0,6.0])
ax4.set_yscale('log')
ax4.set_ylim(bottom=2e-9, top=2e-5)
ax4.xaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
ax4.yaxis.grid(True, ls='-', which='major', color='gray', alpha=0.35)
ax4.tick_params(which='both', right=True, top=True)
ax4.set_title(r'{\bf NuMI Off-axis Flux at ICARUS in} $\nu / \overline{\nu}$ {\bf Mode}', pad=16, fontsize=font)
ax4.set_xlabel(r'{\bf Neutrino Energy} $E_\nu$ [GeV]', fontsize=font)
ax4.set_ylabel(r'$\frac{\mathrm{d} \Phi}{\mathrm{d} E_\nu}$  [$\nu$/m$^2$/POT/GeV]', fontsize=font)
ax4.legend(fontsize=0.7*font, loc='upper right', frameon=False, framealpha=1, edgecolor='k', fancybox=True)

plt.savefig('../Plots/Figures_PDF/Fluxall.pdf', transparent=False, bbox_inches='tight')
plt.savefig('../Plots/Figures_PNG/Fluxall.png', transparent=False, bbox_inches='tight')

plt.show()
