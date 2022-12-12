import numpy as np
from matplotlib import pyplot as plt

magnitudes = [r for r in np.linspace(1.45,25.0,10000)]
def r_surface(
    r_steel_i,
    r_steel_o,
    r_foam_o,
    L,
    k
):
    f1 = 1.0/(k * r_steel_i)
    f2 = np.log(r_steel_o/r_steel_i)/.143
    f3 = np.log(r_foam_o/r_steel_o)/.0003
    f4 = 1.0/(k * r_foam_o)
    return (1.0/(2*np.pi*L)) * (f1 + f2 + f3 + f4)

vals_10 = [r_surface(1.25,1.45,r,168.0,10e-4) for r in magnitudes]
vals_100 = [r_surface(1.25,1.45,r,168.0,10e-1) for r in magnitudes]
val_10 = r_surface(1.25,1.45,11.25,168.0,10e-4)
val_100 = r_surface(1.25,1.45,11.25,168.0,10e-1)

fig, axs = plt.subplots()
axs.plot(magnitudes, vals_10, linestyle='-.',color='k', label=r"$k_{\alpha}=10^{-4}$")
axs.plot(magnitudes, vals_100, linestyle=':',color='b', label=r"$k_{\alpha}=10^{-1}$")
axs.scatter(11.25, val_10, color='r', marker='x', label=r"$R_{\mathrm{surface}}(r^o_{\mathrm{foam}}=11.25)$"+f"={val_10:.4f}")
axs.scatter(11.25, val_100, color='r', marker='x', label=r"$R_{\mathrm{surface}}(r^o_{\mathrm{foam}}=11.25)$"+f"={val_100:.4f}")
axs.set_xlabel(r"$r^o_{\mathrm{foam}}$" + " [cm]")
axs.set_ylabel(r"$R_{\mathrm{surface}}$")
plt.legend()
plt.tight_layout()
plt.savefig("r_surface.png")
plt.show()
