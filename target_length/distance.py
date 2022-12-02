import numpy as np
from matplotlib import pyplot as plt
import csv

m_ar = 66.337316
endf_file = "../data/endf_ar.csv"

energy = []
cross_section = []

with open(endf_file, "r") as file:
    reader = csv.reader(file, delimiter=",")
    for row in reader:
        energy.append(float(row[0]))
        cross_section.append(float(row[1]))


def compute_transmission(density, target_length):
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
    axs.set_xlabel("Energy [eV]")
    axs.set_ylabel("Transmission")
    axs.set_title("Transmission vs. Energy [MeV] " + r"$\rho$" + f" = {density}")
    axs.set_xlim(10000,200000)
    plt.legend()
    plt.tight_layout()
    plt.savefig("transmission_energy.png")
    plt.show()

if __name__ == "__main__":
    target_lengths = [10, 50, 100, 200]
    plot_lengths(1.3954, target_lengths)