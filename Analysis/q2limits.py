import numpy as np
import matplotlib as plt
import matplotlib.pyplot as plt

plt.style.use('sty.mplstyle')

#=================================================
#========   Define constants (** OK **)   ========
#=================================================
me = 0.000511 	# electron mass (GeV)
mmu = 0.105658	# muon mass 	(GeV)
mtau = 1.778	# tau mass		(GeV)

MArgon = 39.95*0.9315		# Argon mass (GeV)
MIron = 55.85*0.9315		# Iron mass (GeV)
MTungsten = 183.84*0.9315	# Tungsten mass (GeV)
MCarbon = 12.01*0.9315		# Carbon mass (GeV)
MNitrogen = 14.01*0.9315	# Nitrogen mass (GeV)
MOxygen = 16.0*0.9315		# Oxygen mass (GeV)
MSilicon = 28.09*0.9315		# Silicon mass (GeV)

#=================================================
#========   Define Processes (** OK **)   ========
#=================================================
Process1 = r'$\nu_e \to \nu_e e^{+} e^{-}$'							# [1] nu_e -> nu_e e+ e-
Process2 = r'$\nu_e \to \nu_e \mu^{+} \mu^{-}$'						# [2] nu_e -> nu_e mu+ mu-
Process3 = r'$\nu_e \to \nu_\mu \mu^{+} e^{-}$'						# [3] nu_e -> nu_mu mu+ e-
Process4 = r'$\bar{\nu}_e \to \bar{\nu}_e e^{+} e^{-}$'				# [4] anti-nu_e -> anti-nu_e e+ e-
Process5 = r'$\bar{\nu}_e \to \bar{\nu}_e \mu^{+} \mu^{-}$'			# [5] anti-nu_e -> anti-nu_e mu+ mu-
Process6 = r'$\bar{\nu}_e \to \bar{\nu}_\mu e^{+} \mu^{-}$'			# [6] anti-nu_e -> anti-nu_mu e+ mu-
Process7 = r'$\nu_\mu \to \nu_\mu e^{+} e^{-}$'						# [7] nu_mu -> nu_mu e+ e-
Process8 = r'$\nu_\mu \to \nu_\mu \mu^{+} mu^{-}$'					# [8] nu_mu -> nu_mu mu+ mu-
Process9 = r'$\nu_\mu \to \nu_e e^{+} \mu^{-}$' 					# [9] nu_mu -> nu_e e+ mu-
Process10 = r'$\bar{\nu}_\mu \to \bar{\nu}_\mu e^{+} e^{-}$' 		# [10] anti-nu_mu -> anti-nu_mu e+ e-
Process11 = r'$\bar{\nu}_\mu \to \bar{\nu}_\mu \mu^{+} \mu^{-}$'	# [11] anti-nu_mu -> anti-nu_mu mu+ mu-
Process12 = r'$\bar{\nu}_\mu \to \bar{\nu}_\e \mu^{+} e^{-}$' 		# [12] anti-nu_mu -> anti-nu_e mu+ e-
Process13 = r'$\nu_\mu \to \nu_\mu \tau^{+} \tau^{-}$' 				# [13] nu_mu -> nu_mu tau+ tau-
Process14 = r'$\nu_\mu \to \nu_\tau \tau^{+} \mu^{-}$' 				# [14] nu_mu -> nu_tau tau+ mu-
Process15 = r'$\nu_e \to \nu_\tau \tau^{+} e^{-}$' 					# [15] nu_e -> nu_tau tau+ e-

m3 = np.array([me, mmu, mmu, me, mmu, me, me, mmu, me, me, mmu, mmu, mtau, mtau, mtau])
m4 = np.array([me, mmu, me, me, mmu, mmu, me, mmu, mmu, me, mmu, me, mtau, mmu, me])
labels = [Process1,Process2,Process3,Process4,Process5,Process6,Process7,Process8,Process9,Process10,Process11,Process12,Process13,Process14,Process15]

#=================================================
#=======   Computing Q2 limits (** OK **)   ======
#=================================================
def Q2limits(Emax, M=MArgon):
	X1max = []
	X1min = []
	Enu = []
	Enuth = ( (m3+m4+M)**2 - M**2 )/2/M
	for i in range(len(m3)):
		eps1 = np.linspace(Enuth[i], Emax, 10000)
		g = eps1/M
		mL = m3[i] + m4[i]
		coef = (mL/eps1)*(mL/eps1)/2
		x1max = 1 - coef*(1+g) + np.sqrt( (1-coef*(1+g))**2 - coef**2*(1+2*g) )
		x1max *= 2*eps1**2/(1+2*g)
		x1min = mL**4/(1+2*g)/x1max
		X1min.append(x1min)
		X1max.append(x1max)
		Enu.append(eps1)
	return Enu, X1min, X1max

#=================================================
#=======   Generating the data (** OK **)   ======
#=================================================
