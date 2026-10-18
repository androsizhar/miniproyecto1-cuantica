from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt

# Compatibilidad: np.trapezoid (numpy >= 2.0) o np.trapz (versiones anteriores)
trapezoid = getattr(np, "trapezoid", None) or np.trapz

# Carpeta donde se guardan las figuras (miniproyecto1-cuantica/figuras)
FIG = Path(__file__).resolve().parent.parent / "figuras"
FIG.mkdir(exist_ok=True)

# =============================================================================
# DATA (mismos valores que el inciso 2)
# =============================================================================
hbar = 1
a = 1
M = 1/8 * np.pi**2
V0 = (hbar*np.pi)**2/(2*M*(2*a)**2)
U0 = 110
mumax = np.sqrt(U0/V0)


# =============================================================================
# Funciones del inciso 2
# =============================================================================
def energias_permisibles(mumax: float, *, eps=1e-12, a=1, hbar=1, M=1/8*np.pi**2):
    '''
    Calcula las energías permisibles del pozo finito simétrico resolviendo las
    ecuaciones (14) y (17). Los estados pares tienen mu en (2n, 2n+1) y los
    impares en (2n+1, 2n+2).

    Returns:
        Tuple(array, array) con las energías pares e impares respectivamente.
    '''
    f = lambda x: x*np.sin(np.pi*x/2) - np.sqrt(mumax**2 - x**2)*np.cos(np.pi*x/2)
    g = lambda x: np.sqrt(mumax**2 - x**2)*np.sin(np.pi*x/2) + x*np.cos(np.pi*x/2)
    raicespar, raicesimpar = [], []
    n = 0
    while 2*n < mumax:
        low, high = 2*n + eps, min(2*n + 1 - eps, mumax)
        if f(low)*f(high) < 0:
            raicespar.append(brentq(f, low, high))
        low += 1; high = min(high + 1, mumax)
        if low > mumax:
            break
        if g(low)*g(high) < 0:
            raicesimpar.append(brentq(g, low, high))
        n += 1
    E1 = (hbar*np.pi)**2/(2*M*(2*a)**2)
    return E1*np.array(raicespar)**2, E1*np.array(raicesimpar)**2


def eigenestados(energias, M, hbar, U0, a, *, n=8, x_ext=3, steps=6001):
    '''
    Evalúa los primeros n eigenestados normalizados (amplitudes analíticas)
    ordenados por energía creciente.

    Returns:
        x (steps,), psi (n, steps), E (n,)
    '''
    E_par, E_imp = energias
    x = np.linspace(-x_ext*a, x_ext*a, steps)
    dentro = np.abs(x) <= a
    psi = np.empty((n, steps))
    E = np.empty(n)
    for j in range(n):
        i = j // 2
        E[j] = E_par[i] if j % 2 == 0 else E_imp[i]
        k = np.sqrt(2*M*E[j])/hbar
        q = np.sqrt(2*M*(U0 - E[j]))/hbar
        A = np.sqrt(q/(a*q + 1))
        if j % 2 == 0:
            B = A*np.cos(k*a)*np.exp(q*a)
            psi[j] = np.where(dentro, A*np.cos(k*x), B*np.exp(-q*np.abs(x)))
        else:
            B = A*np.sin(k*a)*np.exp(q*a)
            psi[j] = np.where(dentro, A*np.sin(k*x), np.sign(x)*B*np.exp(-q*np.abs(x)))
    return x, psi, E


# =============================================================================
# Funciones del inciso 3
# =============================================================================
def evolucion(psi, cn, E, t, *, hbar=1):
    '''
    phi(x,t) = sum_n c_n phi_n(x) exp(-i E_n t / hbar)  (Ec. 2).

    Returns:
        Arreglo complejo de forma (len(t), len(x)).
    '''
    phi = np.zeros((len(t), psi.shape[1]), dtype=complex)
    for n in range(len(E)):
        if cn[n] != 0:
            phi += cn[n]*np.exp(-1j*E[n]*t/hbar)[:, None]*psi[n][None, :]
    return phi


def valor_esperado_x(x, phi):
    '''
    <x>(t) = int x |phi|^2 dx / int |phi|^2 dx
    '''
    rho = np.abs(phi)**2
    return trapezoid(x*rho, x, axis=1)/trapezoid(rho, x, axis=1)


def mapa_densidad(x, t, phi, titulo, archivo, *, a=1):
    '''
    Gráfica bidimensional de |phi(x,t)|^2.
    '''
    fig, ax = plt.subplots(figsize=(7, 4.5))
    im = ax.pcolormesh(x, t, np.abs(phi)**2, shading="auto", cmap="viridis")
    fig.colorbar(im, ax=ax, label=r"$|\phi(x,t)|^2$")
    ax.axvline(-a, color="red", lw=0.8, ls="--")
    ax.axvline(a, color="red", lw=0.8, ls="--")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$t$")
    ax.set_title(titulo)
    fig.tight_layout()
    fig.savefig(archivo, dpi=200)
    plt.show()


def grafica_x(t, xm, titulo, archivo, *, a=1):
    '''
    Gráfica del valor esperado de la posición en función del tiempo.
    '''
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(t, xm, lw=1)
    ax.axhline(-a, color="red", lw=0.6, ls="--")
    ax.axhline(a, color="red", lw=0.6, ls="--")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle x\rangle$")
    ax.set_title(titulo)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(archivo, dpi=200)
    plt.show()


# =============================================================================
# RESULTS: inciso 3(a)
# =============================================================================
if __name__ == "__main__":
    energias = energias_permisibles(mumax)
    N = len(energias[0]) + len(energias[1])          # N = 11
    x, psi, E = eigenestados(energias, M, hbar, U0, a, n=N)

    # Coeficientes de la Ec. (4) y constante de normalización
    coef = np.array([1, 0.7, -1.2j, 2 - 1j])
    A = 1/np.sqrt(np.sum(np.abs(coef)**2))          # A = 1/sqrt(7.93)
    cn = np.zeros(N, dtype=complex)
    cn[:4] = A*coef

    t = np.linspace(0, 10, 1201)
    phi = evolucion(psi, cn, E, t)

    print(f"A = {A:.4f}")
    print(f"Norma en t=0: {trapezoid(np.abs(phi[0])**2, x):.10f}")

    mapa_densidad(x, t, phi, r"Inciso 3(a): $|\phi(x,t)|^2$", FIG / "p1_3a_mapa.png")
    grafica_x(t, valor_esperado_x(x, phi), r"Inciso 3(a): $\langle x\rangle(t)$", FIG / "p1_3a_x.png")
