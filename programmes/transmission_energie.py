import numpy as np
import matplotlib.pyplot as plt

# Constantes physiques
hbar = 1.0545718e-34  # constante de planck réduite
m = 9.10938356e-31    # masse de l'électron

# Paramètres du puits
V0 = 5.0          # profondeur du puits

a = 1E-18          # largeur du puits

# Énergies
E= np.linspace(0.1, 30, 1000)


# Vecteurs d'onde de la région de l'intérieur du puits
k2 = np.sqrt(2 * m * (E + V0)) / hbar


# Formule de transmission
T = 1 / (1 + (V0**2 / (4 * E * (E + V0))) * np.sin(k2 * a)**2)

# Plot
plt.figure(figsize=(10,6))
plt.plot(E, T, label="Transmission en fonction de l'énergie")
plt.xlabel("Énergie")
plt.ylabel("Transmission")
plt.title("Effet Ramsauer-Townsend")
plt.grid(True)
plt.legend()
plt.ylim(0, 1.1)
plt.show()