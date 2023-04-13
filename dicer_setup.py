import numpy as np
from matplotlib import pyplot as plt
import csv

if __name__ == "__main__":

    dicer_flux = []
    dicer_energy = []
    with open("data/dicer_dataset.csv", "r") as file:
        reader = csv.reader(file, delimiter=",")
        for row in reader:
            dicer_flux.append(float(row[1]))
            dicer_energy.append(float(row[0]))

    # fixed quantities
    J = 2.2e7
    E_T = 3.899e-3
    alpha_3 = 5.4876e2
    E_C = 2.07e-4
    C_1 = 1.5e6
    alpha_1 = 1
    C_2 = 1e4
    alpha_2 = 2.3

    entrance_r = 6.0    # cm
    rbb_r = 0.4 # cm
    bc_r = 0.05 # cm
    as_r = .6   # cm

    entrance_area = np.pi*entrance_r*entrance_r
    rbb_area = 2*np.pi*rbb_r*rbb_r
    bc_area = 2*np.pi*bc_r*bc_r
    as_area = 2*np.pi*as_r*as_r

    rbb_flux_red = rbb_area / entrance_area
    bc_flux_red = bc_area / rbb_area
    as_area_red = as_area / bc_area

    def x_factor(energy):
        return alpha_3*(energy - E_C)

    def flux(energy):
        x = x_factor(energy)
        first = (J * energy/E_T)*np.exp(-energy/E_T)
        second = (1 - np.exp(-x)*(1 + x + x*x/2.0))
        third = C_1*np.power(energy, -alpha_1) + C_2*np.power(energy, -alpha_2)
        return first + second*third

    energies = np.linspace(10e-4,10e5,100000)
    entrance_fluxes = flux(energies)
    rbb_fluxes = rbb_flux_red * entrance_fluxes
    bc_fluxes = bc_flux_red * rbb_fluxes

    fig, axs = plt.subplots()
    axs.plot(energies, entrance_fluxes, label='7.9m')
    axs.plot(energies, rbb_fluxes, label='14.65m (rbb)')
    axs.plot(energies, bc_fluxes, label='15.15m (bc)')
    axs.scatter(dicer_energy, dicer_flux, marker='x', color='r', label='data 7.9m')
    axs.set_xscale("log")
    axs.set_yscale("log")
    axs.set_xlim(10e-4,10e5)
    plt.legend()
    plt.tight_layout()
    plt.show()

    energy = 57000.0
    entrance_flux = flux(energy)
    diameters = np.linspace(0.05, 0.4)
    areas = 2*np.pi*diameters*diameters
    reductions = 3600 * areas * areas / entrance_area
    neutrons_hr = 3600 * areas * areas / entrance_area

    fig, axs = plt.subplots()
    axs.plot(diameters, entrance_flux * reductions, linestyle='--', color='k')
    axs.plot(0.05, entrance_flux*reductions[0], marker='x', color='r', label="0.05 cm (BC)")
    axs.plot(0.4, entrance_flux*reductions[-1], marker='x', color='r', label="0.4 cm (RBB)")
    # axs2 = axs.twinx()
    # axs2.plot(diameters, entrance_flux * neutrons_hr)
    # axs2.set_ylabel("N [neutron/hr]")

    axs.set_xlabel("BC radius [cm]")
    axs.set_ylabel(r"$\Phi$" + "(E) [neutrons/hr]")
    axs.set_title("Neutrons/hr at 57 keV vs. BC radius")
    plt.tight_layout()
    plt.show()

    delta_d = 6.0
    
    d_rbb = 1.6
    l_rbb = 30.5
    z_rbb = 1435

    d_as = 2.2
    l_as = 30
    z_as = 1850

    l_target = 200

    # calculate tolerance on BC placement
    delta_rbb_max = delta_d - d_rbb
    delta_z_bc_max = (l_rbb/d_rbb) * delta_rbb_max
    
    bc_rs = np.linspace(0,2.9,10000)
    bc_l_mins = np.zeros(bc_rs.shape)
    for ii in range(len(bc_rs)):
        d_max = (delta_d - bc_rs[ii])/2.0
        bc_l_mins[ii] = (bc_rs[ii]/d_max)*l_target
    d_bc_30 = delta_d / (2*l_target/30 + 1)
    l_bc_rbb = 2 * 1.6 * l_target / (delta_d - 1.6)
    
    fig, axs = plt.subplots()
    axs.plot(bc_rs, bc_l_mins, linestyle='--',color='k')
    axs.plot(d_bc_30, 30, marker='x', color='r', label=f"(r,l)=({d_bc_30:.2f}, {30})")
    axs.plot(1.6, l_bc_rbb, marker='x', color='c', label=f"(r,l)=({1.6}, {l_bc_rbb:.2f})")
    axs.set_xlabel(r"$d_{\mathrm{BC}}$ [cm]")
    axs.set_ylabel(r"$\ell^{\mathrm{min}}_{\mathrm{BC}}$ [cm]")
    axs.set_title("Minimum non-crossing length vs. hole diameter")
    plt.legend()
    plt.tight_layout()
    plt.show()