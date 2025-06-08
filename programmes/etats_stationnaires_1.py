import numpy as np
import matplotlib.pyplot as plt
from scipy.linalg import eigh_tridiagonal


dx = 0.001
nx = int(1 / dx) * 2  # Nombre de point spatiaux
x = np.linspace(0, (nx - 1) * dx, nx) #Vecteur des positions x de taille nx
v0 = -4000 #valeur du potentiel
V = np.zeros(nx)
V[(x >= 0.8) & (x <= 0.9)] = v0  # définition du potentiel pour x appartenant à [0.8, 0.9]

#construction de l'hamiltoniens
diag_principale = 1 / dx**2 + V           # Diagonale principale ("V(n)" + 1/dx**2)
hors_diag = -0.5 / dx**2 * np.ones(nx - 1)  # Terme (-1/2dx**2) hors de la diagonale corresponde au psi(i+-1) 

# Résolution du problème H*psi = E*psi
vp_Energie, vectp_onde = eigh_tridiagonal(diag_principale, hors_diag) #fonction qui renvoie les valeurs propre et vecteur propre importé par scipy

# création des 3 premiers états propres
plt.figure(figsize=(10, 6))
for n in range(3):
    psi_n = vectp_onde[:, n]
    psi_n /= np.sqrt(np.sum(psi_n**2) * dx)  # Normalisation
    plt.plot(x, psi_n + vp_Energie[n], label=f"n={n}, E={vp_Energie[n]:.1f}") # affichage des fonction avec des décalage pour ne pas les superposé


#affichage   
plt.plot(x, V, 'k--', label='Potentiel V(x)')
plt.title("États stationnaires dans un puits rectangulaire")
plt.xlabel("x")
plt.ylabel("ψₙ(x) + énergie")
plt.legend()
plt.grid(True)
plt.show()