import numpy as np
from scipy.integrate import quad

#=================================================
#==========   Define constant values   ===========
#=================================================
aEM = 1/137.036 	# Electromagnetic Coupling
Mproton = 0.938272	# Proton mass in GeV
Mneutron = 0.939565 # Neutron mass in GeV
mup = 2.7928 		# magnetic moment of proton
mun = -1.913 		# magnetic moment of neutron

Ar40_A = 40		# Number of protons and neutrons in the Argon nuclei
Ar40_Z = 18		# Number of protons in the Argon nuclei
W184_A = 184	# Number of protons and neutrons in the Tungsten nuclei
W184_Z = 74		# Number of protons in the Tungsten nuclei

fm_to_GeVem1 = 1/0.1973 # fm to GeV^{-1}
GeVem1_to_fm = 0.1973 # GeV^{-1} to fm

def J1(x): #(** OK **)
	main = 1/3*np.ones_like(x)
	cond = x!= 0
	num = np.sin(x[cond]) - x[cond]*np.cos(x[cond])
	denom = x[cond]**3
	main[cond] = num/denom
	return main

def j1(x): #(** OK **)
	if x == 0:
		main = 1/3
	else:
		num = np.sin(x) - x*np.cos(x)
		denom = x**3
		main = num/denom
	return main

#=================================================
#=========   Argon [Ar40] Form Factors   =========
#=================================================
def Compute_Exp_FF(Energy, Theta, Xsection, A=Ar40_A, Z=Ar40_Z): # (** OK **)
    # This function compute the Mott Cross Section
    # Returns the momentum transfer (GeV) and the form factor
    M = Z*Mproton + (A-Z)*Mneutron 	# Target mass in GeV
    E0 = Energy
    Q = np.ones_like(Theta)
    FF = np.ones_like(Theta)
    Theta *= np.pi/180
    E = E0/( 1 + 2*E0*np.sin(Theta/2)**2/M )	# Energy in GeV
    Q2 = 4*E0*E*np.sin(Theta/2)**2 				# Momentum transfer squared
    Q = np.sqrt(Q2) 							# Momentum transfer in GeV
    Xsec_Mott = aEM**2 * np.cos(Theta/2)**2 / ( 4*E0**2*np.sin(Theta/2)**4 ) / ( 1 + 2*E0*np.sin(Theta/2)**2/M )
    Xsec_Mott *= Z**2 							# Mott Cross Section in GeV^{-2}
    Xsec_Mott *= (GeVem1_to_fm)**2 				# (fm)^2
    FF = np.sqrt( Xsection/Xsec_Mott )
    return Q, FF

#=================================================
#=========   Fourier Bessel  (** OK **)   ========
#=================================================
Ar40_FBSE_R = 9.0*fm_to_GeVem1 # 1/GeV
Ar40_FBSE_a = np.array([0.30451e-1, 0.55337e-1, 0.20203e-1, -0.16765e-1, -0.13578e-1, -0.43204e-4, 0.91988e-3, -0.41205e-3, 0.11971e-3, -0.19801e-4, -0.43204e-5, 0.61205e-5, -0.37803e-5, 0.18001e-5, -0.77407e-6], dtype=np.float64)

Ar40_FBSE_NORM = 0

for i in range(len(Ar40_FBSE_a)):
	n = i+1
	Ar40_FBSE_NORM -= Ar40_FBSE_a[i]*(-1)**n/n**2

def Ar40FBSE_FF(Q2):
	Q = np.sqrt(Q2)
	aux = np.ones_like(Q2)
	cond = Q!=0
	s = np.zeros_like(Q[cond])
	for i in range(len(Ar40_FBSE_a)):
		n = i + 1
		s += ( Ar40_FBSE_a[i]*(-1)**n/( Q[cond]*Ar40_FBSE_R + n*np.pi ) )*( np.sin(Q[cond]*Ar40_FBSE_R)/( Q[cond]*Ar40_FBSE_R - n*np.pi ) )
	s = s*np.pi**2/( Q[cond]*Ar40_FBSE_R )
	aux[cond] = s/Ar40_FBSE_NORM
	return aux

#=================================================
#=====   Three-parameter Fermi  (** OK **)   =====
#=================================================
Ar40_3pF_c = 3.73*fm_to_GeVem1 		# 1/GeV
Ar40_3pF_z = 0.62*fm_to_GeVem1 		# 1/GeV
Ar40_3pF_w = -0.19
Ar40_3pF_R = (3.48*fm_to_GeVem1)*5	# 1/GeV

def Ar40_RHO_FERMI(r):
	num = 1 + Ar40_3pF_w*r**2/Ar40_3pF_c**2
	denom = 1 + np.exp( (r-Ar40_3pF_c)/Ar40_3pF_z )
	return num/denom

def Ar40_RHO_FERMI_aux(r, Q2):
	Q = np.sqrt(Q2)
	if Q == 0:
		return r**2 * Ar40_RHO_FERMI(r)
	else:
		return r**2 * np.sin(Q*r)/(Q*r) * Ar40_RHO_FERMI(r)

Ar40_NORM_FERMI, _ = quad(Ar40_RHO_FERMI_aux, 0, Ar40_3pF_R, args=(0,), limit=1000)

def Ar403pF_FF(Q2):
	if Q2 == 0:
		main = Ar40_NORM_FERMI
	else:
		main, _ = quad(Ar40_RHO_FERMI_aux, 0, Ar40_3pF_R, args=(Q2,), limit=1000)
	return main/Ar40_NORM_FERMI

Ar403pF_FF = np.vectorize(Ar403pF_FF)

#=================================================
#========   Helm Form Factor  (** OK **)   =======
#=================================================
Ar40_HELM_r0 = 0.52	# fm
Ar40_HELM_s = 0.9 	# fm
Ar40_HELM_R0 = np.sqrt( ( 1.23*Ar40_A**(1/3) - 0.6 )**2 + 7*np.pi**2*Ar40_HELM_r0**2/3 - 5*Ar40_HELM_s**2 ) # fm
Ar40_HELM_R0 = Ar40_HELM_R0*fm_to_GeVem1	# 1/GeV
Ar40_HELM_s = Ar40_HELM_s*fm_to_GeVem1		# 1/GeV

def Ar40HELM_FF(Q2):
	Q = np.sqrt(Q2)
	return 3*J1(Q*Ar40_HELM_R0)*np.exp(-Q2*Ar40_HELM_s**2/2)

def Ar40HELM_ff(Q2):
	Q = np.sqrt(Q2)
	return 3*j1(Q*Ar40_HELM_R0)*np.exp(-Q2*Ar40_HELM_s**2/2)

#=================================================
#===   Klein-Nystrand Form Factor  (** OK **)   ==
#=================================================
Ar40_KN_RA = 1.23*Ar40_A**(1/3) 		# fm
Ar40_KN_RA = Ar40_KN_RA*fm_to_GeVem1	# 1/GeV

def Ar40KN_FF(Q2):
	Q = np.sqrt(Q2)
	a = 0.7*fm_to_GeVem1 # 1/GeV
	return 3*J1(Q*Ar40_KN_RA)/( 1 + Q2*a**2 )

def Ar40KN_ff(Q2):
	Q = np.sqrt(Q2)
	a = 0.7*fm_to_GeVem1 # 1/GeV
	return 3*j1(Q*Ar40_KN_RA)/( 1 + Q2*a**2 )

#=================================================
#=   (ad.)  Klein-Nystrand Form Factor  (** OK **)
#=================================================
Ar40_adKN_r0 = 3.427	# fm
Ar40_adKN_ak = 0.7		# fm
Ar40_adKN_RA = np.sqrt( 5*Ar40_adKN_r0**2/3 - 10*Ar40_adKN_ak**2 ) # fm
Ar40_adKN_RA = Ar40_adKN_RA*fm_to_GeVem1 # 1/GeV

def Ar40adKN_FF(Q2):
	Q = np.sqrt(Q2)
	a = 0.7*fm_to_GeVem1 # 1/GeV
	return 3*J1(Q*Ar40_adKN_RA)/( 1 + Q2*a**2 )

def Ar40adKN_ff(Q2):
	Q = np.sqrt(Q2)
	a = 0.7*fm_to_GeVem1 # 1/GeV
	return 3*j1(Q*Ar40_adKN_RA)/( 1 + Q2*a**2 )

#=================================================
#==   Tungsten [W184] Form Factors (** OK **)   ==
#=================================================
W184_WS_sigma = 0.523*fm_to_GeVem1 # 1/GeV
W184_WS_r = 1.126*W184_A**(1/3)*fm_to_GeVem1 # 1/GeV

def W184WS_FF(Q2):
    main = np.ones_like(Q2)
    Q = np.sqrt(Q2)
    cond = Q != 0
    num = np.pi*W184_WS_sigma*np.cosh(np.pi*Q[cond]*W184_WS_sigma)*np.sin(Q[cond]*W184_WS_r)/np.sinh(np.pi*Q[cond]*W184_WS_sigma) - W184_WS_r*np.cos(Q[cond]*W184_WS_r)
    num *= 3*np.pi*W184_WS_sigma
    denom = Q[cond]*W184_WS_r*np.sinh(np.pi*Q[cond]*W184_WS_sigma)
    denom *= (W184_WS_r**2 + np.pi**2*W184_WS_sigma**2)
    main[cond] = num/denom
    return main

#=================================================
#=====   Helm [v1] Form Factor  (** OK **)   =====
#=================================================
W184_HELMv1_c = ( 1.23*W184_A**(1/3) - 0.6 ) # fm
W184_HELMv1_r0 =  0.52 # fm
W184_HELMv1_s = 0.9 # fm
W184_HELMv1_R0 = np.sqrt( W184_HELMv1_c**2 + 7*np.pi**2*W184_HELMv1_r0**2/3 - 5*W184_HELMv1_s**2 ) # fm
W184_HELMv1_R0 = W184_HELMv1_R0*fm_to_GeVem1 # 1/GeV
W184_HELMv1_s = W184_HELMv1_s*fm_to_GeVem1 # 1/GeV

def W184HELMv1_FF(Q2):
	Q = np.sqrt(Q2)
	return 3*J1(Q*W184_HELMv1_R0)*np.exp(-Q2*W184_HELMv1_s**2/2)

def W184HELMv1_ff(Q2):
	Q = np.sqrt(Q2)
	return 3*j1(Q*W184_HELMv1_R0)*np.exp(-Q2*W184_HELMv1_s**2/2)

#=================================================
#=====   Helm [v2] Form Factor  (** OK **)   =====
#=================================================
W184_HELMv2_c = 6.51 # fm
W184_HELMv2_r0 =  0.535 # fm
W184_HELMv2_s = 0.9 # fm
W184_HELMv2_R0 = np.sqrt( W184_HELMv2_c**2 + 7*np.pi**2*W184_HELMv2_r0**2/3 - 5*W184_HELMv2_s**2 ) # fm
W184_HELMv2_R0 = W184_HELMv2_R0*fm_to_GeVem1 # 1/GeV
W184_HELMv2_s = W184_HELMv2_s*fm_to_GeVem1 # 1/GeV

def W184HELMv2_FF(Q2):
	Q = np.sqrt(Q2)
	return 3*J1(Q*W184_HELMv2_R0)*np.exp(-Q2*W184_HELMv2_s**2/2)

def W184HELMv2_ff(Q2):
	Q = np.sqrt(Q2)
	return 3*j1(Q*W184_HELMv2_R0)*np.exp(-Q2*W184_HELMv2_s**2/2)

#=================================================
#===   Klein-Nystrand Form Factor  (** OK **)   ==
#=================================================
W184_KN_RA = 1.23*W184_A**(1/3)			# fm
W184_KN_RA = W184_KN_RA*fm_to_GeVem1	# 1/GeV

def W184KN_FF(Q2):
	Q = np.sqrt(Q2)
	a = 0.7*fm_to_GeVem1 # 1/GeV
	return 3*J1(Q*W184_KN_RA)/( 1 + Q2*a**2 )

def W184KN_ff(Q2):
	Q = np.sqrt(Q2)
	a = 0.7*fm_to_GeVem1 # 1/GeV
	return 3*j1(Q*W184_KN_RA)/( 1 + Q2*a**2 )

#=================================================
#=   (ad.)  Klein-Nystrand Form Factor  (** OK **)
#=================================================
W184_adKN_r0 = 5.3658	# fm
W184_adKN_ak = 0.7		# fm
W184_adKN_RA = np.sqrt( 5*W184_adKN_r0**2/3 - 10*W184_adKN_ak**2 ) # fm
W184_adKN_RA = W184_adKN_RA*fm_to_GeVem1 # 1/GeV

def W184adKN_FF(Q2):
	Q = np.sqrt(Q2)
	a = 0.7*fm_to_GeVem1 # 1/GeV
	return 3*J1(Q*W184_adKN_RA)/( 1 + Q2*a**2 )

def W184adKN_ff(Q2):
	Q = np.sqrt(Q2)
	a = 0.7*fm_to_GeVem1 # 1/GeV
	return 3*j1(Q*W184_adKN_RA)/( 1 + Q2*a**2 )

#=================================================
#=======   Proton & Neutron Form Factors   =======
#=================================================

#=================================================
#===   Altmannshofer Form Factor  (** OK **)   ===
#=================================================
def PROTONGe_Alt(Q2):
    tau = Q2/(4*Mproton**2)
    num = 1 + (-0.19)*tau
    denom = 1 + (11.12)*tau + (15.16)*tau**2 + (21.25)*tau**3
    return num/denom

def PROTONGm_Alt(Q2):
    tau = Q2/(4*Mproton**2)
    num = 1 + (1.09)*tau
    denom = 1 + (12.31)*tau + (25.57)*tau**2 + (30.61)*tau**3
    return (mup*num)/denom

def NEUTRONGm_Alt(Q2):
    tau = Q2/(4*Mneutron**2)
    num = 1 + (8.28)*tau
    denom = 1 + (21.3)*tau + (77)*tau**2 + (238)*tau**3
    return (mun*num)/denom

def GD(Q2):
    return 1/(1 + Q2/0.71)**2

def NEUTRONGe_Alt(Q2):
    tau = Q2/(4*Mneutron**2)
    num = (1.68)*tau*GD(Q2)
    denom = 1 + (3.63)*tau
    return num/denom

#=================================================
#=======   Kelly Form Factor  (** OK **)   =======
#=================================================
def PROTONGe_Kelly(Q2):
    tau = Q2/(4*Mproton**2)
    num = 1 + (-0.24)*tau
    denom = 1 + (10.98)*tau + (12.82)*tau**2 + (21.97)*tau**3
    return num/denom

def PROTONGm_Kelly(Q2):
    tau = Q2/(4*Mproton**2)
    num = 1 + (0.12)*tau
    denom = 1 + (10.97)*tau + (18.86)*tau**2 + (6.55)*tau**3
    return (mup*num)/denom

def NEUTRONGm_Kelly(Q2):
    tau = Q2/(4*Mproton**2)
    num = 1 + (2.33)*tau
    denom = 1 + (14.72)*tau + (24.20)*tau**2 + (84.1)*tau**3
    return (mun*num)/denom

def NEUTRONGe_Kelly(Q2):
    tau = Q2/(4*Mproton**2)
    num = (1.70)*tau*GD(Q2)
    denom = 1 + (3.30)*tau
    return num/denom
