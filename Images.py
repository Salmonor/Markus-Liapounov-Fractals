#Code pour l'affichage des fractales de Liapounov

import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap
import numpy as np
import math
import random

def f(r, x):
    '''(r,float),(x,float)
    Fonction logistique (r doit être un réel de [0,4] et x de [0,1])
    (r,float),(x,float)-> float
    '''
    return r * x * (1 - x)

def liapunov_sab( s, f , a , b , N_e, N):
    '''(s:string,f:fun , a:float , b:float ,N:int)
    Calcule à partir d'une séquence de type AB / BA / AAB / BAABBA 
    et d'une fonction itératrice l'exposant de liapunov associé à
    cette séquence avec N désignant le nombre d'itérations du motif s 
    et N_e le nombre d'itérations du régime transitoire.
    (s:string,f:fun , a:float , b:float ,N:int)->float
    '''
    #Préconditions
    assert isinstance(a,(float,int)) and isinstance(b,(float,int))
    assert isinstance(N,int) and isinstance(s,str)
    
    #initialisation
    u=.5
    lp_e=.0

    for n in range(N_e): #Régime transitoire
        for i in s:
            if i=='A':
                u=f(a,u)
            else:
                u=f(b,u)

    for n in range(N):
        for i in s:
            if i=='A':
                u=f(a,u)
                if a*(1-2*u) != 0:
                    lp_e+=math.log(abs(a*(1-2*u)))
            else:
                u=f(b,u)
                if b*(1-2*u) != 0:
                    lp_e+=math.log(abs(b*(1-2*u)))

    return lp_e/(N*len(s))

def fractale_Markus_Liapunov(s, N, Ne, Ni,x1,x2,y1,y2):
    plt.figure(figsize=(9,6))

    r_a = np.linspace(y1, y2, N)
    r_b = np.linspace(x1, x2, N)
    L = np.zeros((N, N))
    S = np.zeros((N, N))

    for i, a in enumerate(r_a):
        for j, b in enumerate(r_b):
            L[i, j] = liapunov_sab(s, f, a, b, Ne, Ni)
            if L[i,j] < 0:
                S[i,j] = -1

    colors_pos = [(0, 'black'),(0.6, 'blue'), (1, 'white')]
    cmap_pos = LinearSegmentedColormap.from_list('custom_colormap_pos', colors_pos)

    colors_neg = [(0, 'white'), (0.5, 'black'), (1, 'orange')]
    cmap_neg = LinearSegmentedColormap.from_list('custom_colormap_neg', colors_neg)

    plt.imshow(L, origin='lower', cmap=cmap_pos, vmin=0, vmax=np.max(L),
    aspect='auto', extent=(x1, x2, y1, y2))
    plt.imshow(np.ma.masked_where(S != -1, L), origin='lower', cmap=cmap_neg,
     vmin=np.min(L), vmax=0, aspect='auto', extent=(x1, x2, y1, y2))

    plt.xlabel('b')
    plt.ylabel('a')
    plt.title('Fractale de Markus-Liapounov pour le motif ' + s)
    plt.savefig("AB.png", dpi=2000)
    plt.show()

fractale_Markus_Liapunov('AB', 2000, 400, 200,2,4,2,4)

