from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

# Carpeta donde se guardan las figuras (miniproyecto1-cuantica/figuras)
FIG = Path(__file__).resolve().parent.parent / "figuras"
FIG.mkdir(exist_ok=True)

def energias_permisibles(mumax: float,*,eps=1e-12,a=1,hbar=1,M=1/8 * np.pi**2) -> np.ndarray:
    '''
    Se calcula las energías permisibles para el pozo finito. Se sabe que la función tan(pi/2 * x),
    tiene polos en los valores pares de x y siempre es positivo en el intervalo 2n, 2n+1.

    Args:
        mumax:      the maximum value of mumax. (sqrt(U0/V0))
        eps:        a small tolerance parameter
        a:          width of the finite well
        hbar:       the value of Planck reduced constant
        M:          the mass of the particle
    returns:
        Tuple(array,array) with even and odd energies respectively.

    '''
    f = lambda x: x*np.sin(np.pi*x/2) - np.sqrt(mumax**2 - x**2)*np.cos(np.pi*x/2)
    g = lambda x:  np.sqrt(mumax**2 - x**2)*np.sin(np.pi*x/2) +x*np.cos(np.pi*x/2)
    raicespar = []
    raicesimpar = []
    n = 0
    while 2*n < mumax:
        #estados pares
        low = 2*n + 1e-12
        high = min(2*n + 1-eps, mumax)
        if f(low)*f(high) < 0:
            raicespar.append(brentq(f, low, high))
        #estados impares
        low+=1; high= min(high+1,mumax)
        if low>mumax:
            n+=1
            break
        if g(low)*g(high) < 0:
            raicesimpar.append(brentq(g,low,high))
        n += 1
    E1 = (hbar*np.pi)**2/(2*M*(2*a)**2)    
    x_par, x_imp = np.array(raicespar),np.array(raicesimpar)
    energias = (E1*x_par**2, E1*x_imp**2)
    return energias





def plot_energies(energias,*,a=1,ax=None): 
    '''
    Plots allowed energie levels in horizontal lines.
    
    Args:
        energias:       array(s) that contain the energies
        a:              width of the well
        ax:             if included, the axis for graph
    Raises:
        Value error if energies is empty.
    '''
    energias =  np.sort(np.concatenate(energias))
    if len(energias) ==0:
        raise ValueError("Energías must at least have 1 value.")
    if ax is None:
        _, ax = plt.subplots()

    y_prev = -np.inf
    for (index,e) in enumerate(energias,start=1):
        y = max(e, y_prev + 6)
        ax.text(1, y, f"  n={index}", va="center_baseline")
        y_prev = y
    for e in energias:
        plt.hlines(e,-a,a,lw=0.6)
    plt.ylabel("Energía") 
    ax = graph_well(ax)
    plt.title("Valor de las Energías permisibles")
    plt.savefig(FIG / "energiaspermisibles.png", dpi=200)
    plt.show()

def graph_well(ax,*,a=1,U0=110):
    def potential(x: np.ndarray, *, a: float, U0: float) -> np.ndarray:
        """
        Helper function that returns 0 for |x| <= a and U0 otherwise.

        Args:
            x: Positions at which to evaluate the potential.
            a: Half-width of the well.
            U0: Potential height outside the well.

        Returns:
            Array of potential values with the same shape as x.
        """
        return np.where(np.abs(x) <= a, 0.0, U0)
    steps = 2001
    x = np.linspace(-2 * a, 2 * a, steps)
    y = potential(x, a=a, U0=U0)

    ax.set_xlim(-2 * a, 2 * a)
    ax.set_ylim(-0.1 * U0, 1.3 * U0)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$U(x)$")
    ax.grid(alpha=0.3)
    ax.plot(x,y,lw=1)
    return ax

def plot_eigenstates(energias,M,hbar,U0,a,*, n=8):
    '''
    Function that plots eigenstates in a 4x2 grid. To graph, it also calculates the eigenstates.

    Args:
        energias:   tuple that contains even and odd energies in arrays.
        M:          mass of the particle
        hbar:       value of Planck reduced constant
        U0:         value of the potential outside the well
        a:          width of the well
        n:          number of eigenstates to graph (max 8)

    Raises:
        ValueError if there is not enough energies or n>8
    ''' 

    if (len(energias[1]) < n/2):
        raise ValueError("No hay suficientes energías.")

    if (n>8):
        raise ValueError(f"El valor de $n$ es muy grande. El máximo es 8")
    knpar = np.sqrt(2*M*energias[0])/hbar
    qnpar = np.sqrt(2*M*(U0-energias[0]))/hbar
    knimpar = np.sqrt(2*M*energias[1])/hbar
    qnimpar = np.sqrt(2*M*(U0-energias[1]))/hbar
    Anpar = np.sqrt(qnpar/(a*qnpar+1))
    Bnpar = np.cos(knpar*a)/np.exp(-qnpar*a)*Anpar
    Animpar = np.sqrt(knimpar/(a*knimpar-np.tan(a*knimpar)))
    Bnimpar = np.sin(knimpar*a) / np.exp(-qnimpar*a) * Animpar
    steps= 2001
    fig, axes = plt.subplots(4, 2, figsize=(10, 12), sharex=True)
    x1 = np.linspace(-2*a,-a,steps)
    x2 = np.linspace(-a,a,steps)
    x3 = np.linspace(a,2*a,steps)
    count = 0
    for count,ax in zip(range(n),axes.flat):
        i = count // 2
        if count % 2 == 0:
            # ESTADO PAR
            psi = np.concatenate([Bnpar[i]*np.exp(qnpar[i]*x1),
                                  Anpar[i]*np.cos(knpar[i]*x2),
                                  Bnpar[i]*np.exp(-qnpar[i]*x3)])
        else:
            # ESTADO IMPAR
            psi = np.concatenate([-Bnimpar[i]*np.exp(qnimpar[i]*x1),
                                  Animpar[i]*np.sin(knimpar[i]*x2),
                                  Bnimpar[i]*np.exp(-qnimpar[i]*x3)])
        ax.plot(np.concatenate([x1, x2, x3]), psi, label=f"n={count+1}")
        ax.set_title(f"n={count+1}")
        # well line
        ax.axvline(-a, color="red", lw=0.4, ls="--")
        ax.axvline(a, color="red", lw=0.4, ls="--")
    fig.suptitle(f"Primeros {n} Eigenestados del pozo finito", fontsize=14)
    fig.tight_layout()
    for ax in axes[-1]:
        ax.set_xlabel("x")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\psi_n(x)$")

    plt.legend()
    plt.savefig(FIG / "eigenestados.png", dpi=200)
    plt.show()


def eigenestados(
    energias: tuple[np.ndarray, np.ndarray],
    M: float,
    hbar: float,
    U0: float,
    a: float,
    *,
    n: int = 8,
    x_ext: float = 3,
    steps: int = 6001,
) -> tuple[np.ndarray, np.ndarray]:
    """
    Computes eigenstate for a finite potential well and returns them in an increasing energy order.

    Args:
        energias:   tuple that contains even and odd energies in arrays.
        M:          mass of the particle
        hbar:       value of Planck reduced constant
        U0:         value of the potential outside the well
        a:          width of the well
        n:          number of eigenstates to graph (max 8)
        x_ext:      the máximum value of exterior that the eigenstates takes
        steps:      the total amount of points evaluated between -x_ext and x_ext

    Returns:
        A tuple x,psi where x is an array of length steps and psi are the eigenstates in a shape (n,steps).

    Raises:
        ValueError: If there are not enough energies for n states, or if
            any energy lies outside (0, U0).
    """
    E_par, E_imp = (np.asarray(e, dtype=float) for e in energias)
    if len(E_par) < (n + 1) // 2 or len(E_imp) < n // 2:
        raise ValueError("No hay suficientes energías.")
    if np.any((E_par <= 0) | (E_par >= U0)) or np.any((E_imp <= 0) | (E_imp >= U0)):
        raise ValueError("Las energías deben cumplir 0 < E < U0.")

    x = np.linspace(-x_ext * a, x_ext * a, steps)
    dentro = np.abs(x) <= a
    psi = np.empty((n, steps))

    for j in range(n):
        i = j // 2
        E = E_par[i] if j % 2 == 0 else E_imp[i]
        k = np.sqrt(2 * M * E) / hbar
        q = np.sqrt(2 * M * (U0 - E)) / hbar
        A = np.sqrt(q / (a * q + 1))
        if j % 2 == 0:
            B = A * np.cos(k * a) * np.exp(q * a)
            psi[j] = np.where(dentro, A * np.cos(k * x), B * np.exp(-q * np.abs(x)))
        else:
            B = A * np.sin(k * a) * np.exp(q * a)
            psi[j] = np.where(dentro, A * np.sin(k * x),
                              np.sign(x) * B * np.exp(-q * np.abs(x)))

    return x, psi

def is_ortonormal(energias,M,hbar,U0,a,*,n=11):
    ''' 
    Function that shows a plot of S, the overlap matrix, that is used in quantum mechanics to measure 
    the orthogonality of the functions. The S_ij element is given by the internal product of psi_i and psi_j.

    Args:
        energias:   tuple that contains even and odd energies in arrays.
        M:          mass of the particle
        hbar:       value of Planck reduced constant
        U0:         value of the potential outside the well
        a:          width of the well
        n:          number of eigenstates to includei n matrix
    
    
    '''
    x, psi = eigenestados(energias,M,hbar,U0,a,n=n)
    S = np.trapezoid(psi[:, None, :] * psi[None, :, :], x, axis=-1)
    fig, ax = plt.subplots()
    im = ax.imshow(S, cmap="Greys", vmin=-1, vmax=1)
    fig.colorbar(im, ax=ax, label=r"$\langle\psi_i|\psi_j\rangle$")
    ax.set_xticks(range(len(S)), labels=range(1, len(S)+1))
    ax.set_yticks(range(len(S)), labels=range(1, len(S)+1))
    ax.set_xlabel("j")
    ax.set_ylabel("i")
    ax.set_title(f"Matriz de traslape $S$")
    plt.savefig(FIG / "overlapS.png", dpi=200)
    plt.show()


# DATA:
hbar = 1
a = 1
M = 1/8 * np.pi**2
V0 = (hbar*np.pi)**2/(2*M*(2*a)**2)    
U0 = 110
mumax= np.sqrt(U0/V0)

# RESULTS:
energias = energias_permisibles(mumax)
print("Energías pares:  ", np.round(energias[0], 4))
print("Energías impares:", np.round(energias[1], 4))
plot_energies(energias)
plot_eigenstates(energias,M,hbar,U0,a)
is_ortonormal(energias,M,hbar,U0,a)
