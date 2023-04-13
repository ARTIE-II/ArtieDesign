import numpy as np
from matplotlib import pyplot as plt
import csv

m_ar = 66.337316
endf_file = "../data/endf_natar.csv"

energy = []
cross_section = []

with open(endf_file, "r") as file:
    reader = csv.reader(file, delimiter=",")
    for row in reader:
        energy.append(float(row[0])*1000)
        cross_section.append(float(row[1]))


def compute_transmission(density, target_length):
    print(density * target_length / m_ar)
    transmissions = [
        np.exp(-density * target_length * cs / m_ar )
        for cs in cross_section
    ]
    return transmissions

def plot_lengths(density, target_lengths):
    fig, axs = plt.subplots(figsize=(10,6))
    for length in target_lengths:
        transmissions = compute_transmission(density, length)
        axs.plot(energy, transmissions, label=f"d = {length} cm", linestyle='--')
    axs.set_xlabel("Energy [eV]", fontsize=12)
    axs.set_ylabel("Transmission", fontsize=12)
    axs.set_title("Transmission vs. Energy [eV] " + r"$\rho$" + f" = {density}")
    axs.set_xlim(2e4,2e5)
    labels = axs.get_xticks() / 1e6
    axs.set_xscale("log")
    #axs.set_xticklabels(labels, rotation=45, ha='right', fontsize=12)
    plt.legend(fontsize=12)
    plt.tight_layout()
    plt.savefig("transmission_energy.png")
    plt.show()

if __name__ == "__main__":
    target_lengths = [1, 15, 200]
    rho_liquid = 1.3954
    rho_gas = 0.0017984
    plot_lengths(rho_liquid, target_lengths)