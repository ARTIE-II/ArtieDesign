import numpy as np
from matplotlib import pyplot as plt
magnitudes = [r for r in np.linspace(1.45,25.0,10000)]
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
 

vals_10 = [r_surface(1.25,1.45,r,168.0,10e-4) + r_face(0.8, 1.25, 0.00762, 1.0, 10e-4) for r in magnitudes]
vals_100 = [r_surface(1.25,1.45,r,168.0,10e-1) + r_face(0.8, 1.25, 0.00762, 1.0, 10e-1) for r in magnitudes]
val_10 = r_surface(1.25,1.45,11.25,168.0,10e-4) + r_face(0.8, 1.25, 0.00762, 1.0, 10e-4) 
val_100 = r_surface(1.25,1.45,11.25,168.0,10e-1) + r_face(0.8, 1.25, 0.00762, 1.0, 10e-1) 


factor = 3.43 * (1.25 * 1.25) * 168.0

t_alpha_10 = factor * np.array(vals_10)
t_alpha_100 = factor * np.array(vals_100)

t_10 = factor * val_10
t_100 = factor * val_100


fig, axs = plt.subplots()
#axs.plot(magnitudes, t_alpha_10, linestyle='-.',color='k', label=r"$k_{\alpha}=10^{-4}$")
axs.plot(magnitudes, t_alpha_100, linestyle=':',color='b', label=r"$k_{\alpha}=10^{-1}$")
#axs.scatter(11.25, t_10, color='r', marker='x', label=r"$\tau_{\alpha}(r^o_{\mathrm{foam}}=11.25)$"+f"={t_10:.4f}")
axs.scatter(11.25, t_100, color='r', marker='x', label=r"$\tau_{\alpha}(r^o_{\mathrm{foam}}=11.25)$"+f"={t_100:.4f}")
axs.set_xlabel(r"$r^o_{\mathrm{foam}}$" + " [cm]")
axs.set_ylabel(r"$\tau_{\alpha}$" + " [s]")
plt.legend()
plt.tight_layout()
plt.savefig("r_total.png")
plt.show()