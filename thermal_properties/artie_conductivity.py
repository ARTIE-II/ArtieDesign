import numpy as np
from matplotlib import pyplot as plt
from thermal import *

inner_r = 1.25
outer_r = 11.25
foam_r = np.linspace(0.1, 20, 1000)
length = 200
steel_width = 1.0
kapton_width = 0.00762

taus = []
for r in foam_r:
    artie = ArtieI(inner_r,outer_r,r,length)
    taus.append(tau_alpha(lar_rho,lar_alpha,artie.p_over_v(lar_boil,300)))

fig, axs = plt.subplots()
axs.plot(foam_r, taus)
axs.set_xlabel("Foam radius [cm]")
axs.set_ylabel(r"$\tau_{\alpha}$ [s]")
plt.show()

