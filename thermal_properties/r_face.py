import numpy as np
from matplotlib import pyplot as plt

magnitudes = [pow(10,n) for n in np.linspace(-4,1,10000)]

def r_face(
    r_kap,
    r_steel,
    z_kap,
    z_steel,
    k
):
    c_kap = 1.0 / (np.pi * r_kap * r_kap)
    c_steel = 1.0 / (np.pi * (r_steel * r_steel) - np.pi * (r_kap * r_kap))
    f_kap = (z_kap/.003 + 2.0/k)
    f_steel = (z_steel/.143 + 2.0/k)
    return c_kap*f_kap + c_steel*f_steel

vals = [r_face(0.8, 1.25, 0.00762, 1.0,k) for k in magnitudes]
val_10 = r_face(0.8, 1.25, 0.00762, 1.0,10e-4)
val_100 = r_face(0.8, 1.25, 0.00762, 1.0,10e-1)


fig, axs = plt.subplots()
axs.plot(magnitudes, vals, linestyle='--',color='k')
axs.scatter(10e-4, val_10, color='r', marker='x', label=r"$R_{\mathrm{face}}(k=10^{-4})$"+f"={val_10:.4f}")
axs.scatter(10e-1, val_100, color='r', marker='x', label=r"$R_{\mathrm{face}}(k=10^{-1})$"+f"={val_100:.4f}")
axs.set_xscale("log")
axs.set_xlabel(r"$k_{\mathrm{air}} = k_{\mathrm{LAr}}$" + " [" + r"$W cm^{-2} K^{-1}$" + "]")
axs.set_ylabel(r"$R_{\mathrm{face}}$")
plt.legend()
plt.tight_layout()
plt.savefig("r_face.png")
plt.show()
