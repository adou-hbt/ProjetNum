import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

def init():     # ça c'est pour initialiser une liste ou un truc
    line.set_data([], [])
    return line,

def animate(j):     # Pour animer
    line.set_data(x_array, final_density[j,:]) #Crée un graphique pour chaque densite sauvegarde
    return line,


# INITIALISATION DES DONNEES-------------------------
dt=1E-7
dx=0.01            
nx=int(1/dx)*2      # Nbr de points dans l'espace
nt=90000            # "Nombre"

n_frames=int(nt/1000)+1     #nombre d image dans notre animation
s=dt/(dx**2)        # Coef de la dérivée seconde

v0=-4000        # Profondeur du potentiel
e=5             # Valeur du rapport E/V0
E=e*v0          # Energie du paquet d'onde
k=math.sqrt(2*abs(E))

x_array = np.linspace(0, (nx - 1) * dx, nx)
V_potential = np.zeros(nx)
V_potential[:] = 0


#CALCUL DES DONNEES----------------------------------
xc=0.6  # Centre du paquet d'onde
sigma=0.05  # Largeur du paquet d'onde

normalisation=1/(math.sqrt(sigma*math.sqrt(math.pi)))
wp_gauss = normalisation * np.exp(1j * k * x_array - ((x_array - xc) ** 2) / (2 * (sigma ** 2)))        #Initialisation de l'équation d'onde

#wave packet Real part 
wp_re=np.zeros(nx)
wp_re[:]=np.real(wp_gauss[:])

#wave packet Imaginary part 
wp_im=np.zeros(nx)
wp_im[:]=np.imag(wp_gauss[:])

density = np.zeros((nt,nx))
density[0,:] = np.absolute(wp_gauss[:]) ** 2

final_density =np.zeros((n_frames,nx))


# FONCTION D'ONDE DANS LE TEMPS----------------
psi = wp_gauss.copy()
for t in range(1, nt): # Boucle dans le temps
    psi_new = psi.copy()
    for i in range(1, nx - 1): # Boucle dans l'espace
        laplacian = psi[i+1] - 2*psi[i] + psi[i-1]
        potential_term = V_potential[i] * psi[i]
        psi_new[i] = psi[i] + 1j * s * laplacian - 1j * dt * potential_term

    psi = psi_new
    
    
    density[t, :] = np.abs(psi) ** 2

    # animation : une image tous les 1000 pas
    if t % 1000 == 0:
        frame_index = t // 1000
        final_density[frame_index, :] = density[t, :]

#Algo devant retourner la densité de probabilité de présence de la particule à différents instants



# ANIMATION ----------------------------------------
plot_title = "E/Vo="+str(e)

fig = plt.figure() # initialise la figure principale
line, = plt.plot([], [])
plt.ylim(-1,10)
plt.xlim(0,2)
plt.plot(x_array,V_potential,label="Potentiel")
plt.title(plot_title)
plt.xlabel("x")
plt.ylabel("Densité de probabilité de présence")
plt.legend() #Permet de faire apparaitre la legende

ani = animation.FuncAnimation(fig,animate,init_func=init, frames=n_frames, blit=False, interval=100, repeat=False)
file_name = 'paquet_onde_e='+str(e)+'.mp4'
ani.save(file_name, writer = animation.FFMpegWriter(fps=120, bitrate=5000))
plt.show()
