import numpy as np
import matplotlib.pyplot as plt

hbar = 1.0545718e-34  # constante de Planck réduite (J·s)
m = 9.10938356e-31    # masse de l'électron (kg)
eV = 1.60218e-19       # 1 eV en Joules

# Paramètres du potentiel
V0 = 10 * eV           # Hauteur du puits de potentiel (J)
a = 1e-10              # Largeur du puits (m)

# Énergies de 0.01 eV à 20 eV
E = np.linspace(0.01, 20, 1000) * eV

# Transmission à travers un puits carré (modèle simplifié)
def transmission(E):
    T = np.zeros_like(E)

    # Cas 1 : E < V0 (barrière)
    mask1 = E < V0
    kappa = np.sqrt(2 * m * (V0 - E[mask1])) / hbar
    T[mask1] = 1 / (1 + (V0**2 * np.sinh(kappa * a)**2) / (4 * E[mask1] * (V0 - E[mask1])))

    # Cas 2 : E >= V0 (puits)
    mask2 = E >= V0
    k = np.sqrt(2 * m * (E[mask2] - V0)) / hbar
    T[mask2] = 1 / (1 + (V0**2 * np.sin(k * a)**2) / (4 * E[mask2] * (E[mask2] - V0)))

    return T

T_vals = transmission(E)

# Tracer
plt.figure(figsize=(10, 6))
plt.plot(E / eV, T_vals, color='blue')

plt.title("Effet Ramsauer-Townsend : la transmission en fonction de l'énergie")

plt.ylim(0,1)
plt.xlim(0,100)
plt.xlabel("Energie E (J)")
plt.ylabel("Transimission T")

plt.axhline(1, color='gray', linestyle='--', linewidth=0.5)
plt.tight_layout()
#plt.legend()
plt.grid(True)
plt.show()