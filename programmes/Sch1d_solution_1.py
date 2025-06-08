import numpy as np
import math
import matplotlib.pyplot as plt
import matplotlib.animation as animation

#~CALCUL ETATS STATIONNAIRES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
def Energie(n_energies,m,L):
    hbar = 1.055e-34       # constante de Planck réduite
    eV = 1.602e-19
    
    energie = np.zeros(n_energies)
    for n in range(1 ,n_energies+1):
        En = pow((hbar * n * np.pi),2) / (2 * m * pow(L*pow(10,-9),2))
        En_eV = En/eV
        energie[n-1] = En_eV + v0
    return energie

#~INITIALISATION DES VARIABLES~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
#Variables de l'environnement (espace, temps et animation)
dt=1E-7
dx=0.001
nx=int(1/dx)*2
nt=90000 
nd=int(nt/1000)+1
n_frame = nd

#Variables du paquet d'ondes gaussien
s=dt/(dx**2)
xc=0.6
sigma=0.05
A=1/(math.sqrt(sigma*math.sqrt(math.pi)))

#Variables de l'équation de schrodinger
v0=-4000
e=5#Valeur du rapport E/V0
E=e*v0
m = 9.109e-31          # masse d’un éléctron
L = 0.6

# Initialisation des tableaux
o = np.linspace(0, (nx - 1) * dx, nx)
V = np.zeros(nx)
#V[o >= 1] = v0  # Potentiel
V[(o >= 0.8) & (o<=0.9)] = v0  # Potentiel

Etats = Energie(1,m,L)

#~ANIMATION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
plot_title = "Marche Ascendante avec E/Vo="+str(e)
fig = plt.figure() # initialise la figure principale
colors = ['red', 'green', 'blue', 'orange', 'purple']
lines = [plt.plot([], [], label=f"État {i+1}", color=colors[i % len(colors)])[0] for i in range(len(Etats))]
plt.ylim(-10,13)
plt.xlim(0,2)

plt.plot(o,V,label="Potentiel")

for i in range(len(Etats)):
    plt.plot([0,2],[Etats[i],Etats[i]], color=colors[i])

plt.title(plot_title)
plt.xlabel("x")
plt.ylabel("Densité de probabilité de présence")
plt.legend() #Permet de faire apparaitre la legende

ani = []
all_final_densites = []

for j in Etats:   
    k=math.sqrt(2*abs(j))   #nombre d'éléctrons
    #Equation Schrödinger
    cpt = A * np.exp(1j * k * o - ((o - xc) ** 2) / (2 * (sigma ** 2)))
    re=np.zeros(nx)
    re[:]=np.real(cpt[:])
    b=np.zeros(nx)
    im=np.zeros(nx)
    im[:]=np.imag(cpt[:])

    #Densité
    densite=np.zeros((nt,nx))
    densite[0,:] = np.absolute(cpt[:]) ** 2
    final_densite=np.zeros((n_frame,nx))
    all_final_densites.append(final_densite)

    #~MODELISATION FONCTION~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    it=0
    for i in range(1, nt):
        if i % 2 != 0:
            b[1:-1]=im[1:-1]
            im[1:-1] = im[1:-1] + s * (re[2:] + re[:-2]) - 2 * re[1:-1] * (s + V[1:-1] * dt)
            densite[i,1:-1] = re[1:-1]*re[1:-1] + im[1:-1]*b[1:-1]
        else:
            re[1:-1] = re[1:-1] - s * (im[2:] + im[:-2]) + 2 * im[1:-1] * (s + V[1:-1] * dt)

    for i in range(1,nt):
        if((i-1)%1000==0):
            it+=1
            final_densite[it][:]=densite[i][:]


def init():
    for line in lines:
        line.set_data([], [])
    return lines

def animate(j):
    for i, line in enumerate(lines):
        line.set_data(o, all_final_densites[i][j, :])
    return lines

#file_name = 'paquet_onde_e='+str(e)+'.mp4'
ani.append(animation.FuncAnimation(fig,animate,init_func=init, frames=nd, blit=False, interval=100, repeat=False))
#ani.save(file_name, writer = animation.FFMpegWriter(fps=120, bitrate=5000))
plt.legend()
plt.grid(True)
plt.show()
