"""
Class for calculating thermal properties of 
the ARTIE-II target.
"""
import numpy as np
from matplotlib import pyplot as plt

"""
Various thermal conductivities in units of W cm^-1 K^-1
and heat transfer coefficients in units of W cm^-2 K^-1
"""
sigma_polyurethane = 0.0003
sigma_steel = 0.143
sigma_aluminum = 2.37
sigma_kapton = 0.003
sigma_air = 0.00026
sigma_vacuum_rough = 0.000212
sigma_vacuum_high = 00.0000111
k_air = 1.0
k_LAr = 1.0

"""
Various thermal properties of LAr
"""
sigma_LAr = 1.26e-3     # W cm^-1 K^-1
gamma_LAr = 161.14      # J g^-1
c_p_LAr = 1.117         # J g^-1 K^-1
eta_LAr = 2.707e-5      # g cm^-1 s^-1
nu_LAr = 1.939e-5       # cm^2 s^-1
beta_LAr = 1.17e-2      # K^-1
alpha_LAr = 8.08e-4     # cm^2 s^-1
Pr_LAr = 1.71e-2        # no units
c_p_over_c_v_LAr = 1.67 # no units
triple_point_LAr = 83.8 # K
boiling_point_LAr = 87.3 # K
freezing_point_LAr = 83.8 # K
rho_triple_point_LAr = 1.417    # g cm^-3
rho_boiling_point_LAr = 1.396   # g cm^-3
rho_LAr = 1.3973

"""
Various constants
"""
g = 980.0   # cm s^-2
T_inf = 300 # K
Rayleigh_const = (rho_LAr * beta_LAr * g * (T_inf - 85.0))/(eta_LAr * alpha_LAr)

def inches_to_cm(inch):
    return inch * 2.54

ArtieII_config = {
    'r_LAr':    inches_to_cm(0.87/2),
    'r_steel_o':inches_to_cm((1.0)/2),
    'r_vac_o':  inches_to_cm((3.87/2)),
    'r_annulus_inner_i': inches_to_cm(4.0/2),
    'r_annulus_outer_i': inches_to_cm(5.87/2),
    'r_annulus_outer_o': inches_to_cm(6.0/2),
    'r_foam_o': inches_to_cm(14.0/2),
    'l_LAr': inches_to_cm(62.442),
    'l_annulus': inches_to_cm(54.808),
    'z_al_LAr': inches_to_cm(.079)
}

class ARTIE:
    """
    """
    def __init__(self,
        cfg:    dict=ArtieII_config
    ):
        self.cfg = cfg
        self.r_LAr = self.cfg['r_LAr']

        self.r_steel_o = self.cfg['r_steel_o']
        self.delta_r_steel = self.r_steel_o - self.r_LAr

        self.r_vac_o = self.cfg['r_vac_o']
        self.delta_r_vac = self.r_vac_o - self.r_steel_o
        
        self.r_annulus_inner_i = self.cfg['r_annulus_inner_i']
        self.delta_r_annulus_inner = self.r_annulus_inner_i - self.r_vac_o

        self.r_annulus_outer_i = self.cfg['r_annulus_outer_i']
        self.delta_r_annulus = self.r_annulus_outer_i - self.r_annulus_inner_i

        self.r_annulus_outer_o = self.cfg['r_annulus_outer_o']
        self.delta_r_annulus_outer = self.r_annulus_outer_o - self.r_annulus_outer_i

        self.r_foam_o = self.cfg['r_foam_o']
        self.delta_r_foam = self.r_foam_o - self.r_annulus_outer_o

        self.l_LAr = self.cfg['l_LAr']
        self.l_annulus = self.cfg['l_annulus']
        self.z_al_LAr = self.cfg['z_al_LAr']

        self.V_LAr = np.pi * self.r_LAr * self.r_LAr * self.l_LAr

    def Rayleigh_number(self):
        return Rayleigh_const * np.power(self.l_LAr, 3)
    
    def R_surf_LAr(self):
        return (1.0/(2*np.pi*self.l_LAr)) * (1.0/(k_LAr * self.r_LAr))
    
    def R_surf_steel(self):
        return (1.0/(2*np.pi*self.l_LAr)) * (np.log(self.r_steel_o/self.r_LAr)/sigma_steel)
    
    def R_surf_vac(self):
        return (1.0/(2*np.pi*self.l_LAr)) * (np.log(self.r_vac_o/self.r_steel_o)/sigma_vacuum_high)
    
    def R_surf_annulus_inner(self):
        return (1.0/(2*np.pi*self.l_LAr)) * (np.log(self.r_annulus_inner_i/self.r_vac_o)/sigma_steel)
    
    def R_surf_annulus_LAr(self):
        return (1.0/(2*np.pi*self.l_LAr)) * ((1.0/(k_LAr * self.r_annulus_inner_i)) + (1.0/(k_LAr * self.r_annulus_outer_i)))
    
    def R_surf_annulus_outer(self):
        return (1.0/(2*np.pi*self.l_LAr)) * (np.log(self.r_annulus_outer_o/self.r_annulus_outer_i)/sigma_steel)

    def R_surf_foam(self):
        return (1.0/(2*np.pi*self.l_LAr)) * (np.log(self.r_foam_o/self.r_annulus_outer_o)/sigma_polyurethane)
    
    def R_surf_air(self):
        return (1.0/(2*np.pi*self.l_LAr)) * (1.0/(k_air * self.r_foam_o))
    
    def R_surf_total(self):
        return (
            self.R_surf_LAr() + 
            self.R_surf_steel() + 
            self.R_surf_vac() + 
            self.R_surf_annulus_inner() + 
            self.R_surf_annulus_LAr() + 
            self.R_surf_annulus_outer() + 
            self.R_surf_foam() + 
            self.R_surf_air()
        )
    
    def R_face_LAr(self):
        return (1.0/(np.pi * self.r_LAr*self.r_LAr)) * (self.z_al_LAr / sigma_aluminum + 1.0/k_LAr + 1.0/k_air)
    
    def R_face_vac(self):
        return (1.0/(np.pi * (self.r_annulus_outer_o * self.r_annulus_outer_o - self.r_LAr*self.r_LAr))) * (self.z_al_LAr / sigma_aluminum)
    
    def R_face_foam(self):
        return (1.0/(np.pi * (self.r_foam_o * self.r_foam_o - self.r_annulus_outer_o * self.r_annulus_outer_o))) * (1.0/k_air)

    def R_face_total(self):
        return 2 * (
            self.R_face_LAr() + 
            self.R_face_vac() + 
            self.R_face_foam()
        )
    
    def tau_alpha(self):
        return (gamma_LAr * rho_LAr * np.pi * self.r_LAr*self.r_LAr * self.l_LAr/(300-85)) * (self.R_surf_total() + self.R_face_total())

    def tau_alpha_vac_r(self, r_val):
        self.r_vac_o = self.r_steel_o + r_val
        self.r_annulus_inner_i = self.r_vac_o + self.delta_r_annulus_inner
        self.r_annulus_outer_i = self.r_annulus_inner_i + self.delta_r_annulus
        self.r_annulus_outer_o = self.r_annulus_outer_i + self.delta_r_annulus_outer
        self.r_foam_o = self.r_annulus_outer_o + self.delta_r_foam
        tau = self.tau_alpha()
        self.r_vac_o = self.cfg['r_vac_o']       
        self.r_annulus_inner_i = self.cfg['r_annulus_inner_i']
        self.r_annulus_outer_i = self.cfg['r_annulus_outer_i']
        self.r_annulus_outer_o = self.cfg['r_annulus_outer_o']
        self.r_foam_o = self.cfg['r_foam_o']
        return tau

    def rate(self):
        return (artie.V_LAr / (artie.tau_alpha()/3600)) / 1000

    def plot_vac_r(self):
        steps = 10000
        delta_r = (self.r_vac_o + 10 - self.r_steel_o) / steps
        r_vals = [delta_r * ii for ii in range(steps)]

        v1_r_val = self.r_vac_o - self.r_steel_o
        v1_tau_alpha = self.tau_alpha()
        v1_rate = self.rate()

        v1_r_zero = 0
        v1_tau_alpha_zero = self.tau_alpha_vac_r(0)
        v1_rate_zero = self.V_LAr / (self.tau_alpha_vac_r(0)/3600) / 1000

        taus = []
        boil = []
        for r_val in r_vals:
            taus.append(self.tau_alpha_vac_r(r_val))
            boil.append(self.V_LAr / (self.tau_alpha_vac_r(r_val)/3600) / 1000)
        fig, axs = plt.subplots()
        axs.plot(
            r_vals, taus,
            linestyle='--',
            color='k'
        )
        axs.plot(
            v1_r_val, v1_tau_alpha,
            marker='x',
            linestyle='',
            color='r',
            label=r"$\tau_{\alpha}(r=$"+f"{self.r_vac_o:.2f})={v1_tau_alpha:.2f}"
        )
        axs.plot(
            v1_r_zero, v1_tau_alpha_zero,
            marker='x',
            linestyle='',
            color='g',
            label=r"$\tau_{\alpha}(r=0) = $" + f"{v1_tau_alpha_zero:.2f}"
        )
        axs.set_xlabel(r"$\Delta r$" + " [cm]")
        axs.set_ylabel(r"$\tau_{\alpha}$" + " [s]")
        axs.set_title(r"$\Delta r$" + " Vacuum shell vs. " + r"$\tau_{\alpha}$")
        plt.legend()
        plt.tight_layout()
        plt.show()



if __name__ == "__main__":
    artie = ARTIE()
    print(ArtieII_config)
    print(artie.R_surf_total())
    print(artie.tau_alpha())
    print(artie.tau_alpha()/3600)
    print(artie.rate())

    artie.plot_vac_r()