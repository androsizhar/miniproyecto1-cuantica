import numpy as np
from scipy.optimize import brentq

hbar = 1
a = 1
M = 1/8 * np.pi**2
V0 = (hbar*np.pi)**2/(2*M*(2*a)**2)
U0 = 110
mumax = np.sqrt(U0/V0)


def energias_permisibles(mumax: float, *, eps=1e-12, a=1, hbar=1, M=1/8*np.pi**2):
    '''
    Calcula las energías permisibles del pozo finito simétrico resolviendo las
    ecuaciones (14) y (17).

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


def energias_unificadas(mumax: float, *, V0=1):
    '''
    Resuelve la ecuación característica unificada
        mu = n - (2/pi) arcsin(mu/mumax),   n = 1, 2, ...
    Hay una raíz por cada n mientras F(mumax) = mumax + 1 > n.

    Returns:
        Array con las energías E_n = V0 mu_n^2 en orden creciente.
    '''
    F = lambda mu, n: mu + (2/np.pi)*np.arcsin(mu/mumax) - n
    mus = []
    n = 1
    while F(mumax, n) > 0:
        mus.append(brentq(F, 0, mumax, args=(n,), xtol=1e-14))
        n += 1
    return V0*np.array(mus)**2


if __name__ == "__main__":
    E_par, E_imp = energias_permisibles(mumax)
    E_paridad = np.sort(np.concatenate((E_par, E_imp)))
    E_unif = energias_unificadas(mumax, V0=V0)

    print(f"N (paridad) = {len(E_paridad)},  N (unificada) = {len(E_unif)},  ceil(mumax) = {int(np.ceil(mumax))}")
    print(f"{'n':>3} {'Ecs. 14 y 17':>14} {'Unificada':>12} {'Diferencia':>12}")
    for n, (e1, e2) in enumerate(zip(E_paridad, E_unif), start=1):
        print(f"{n:>3} {e1:>14.4f} {e2:>12.4f} {abs(e1 - e2):>12.1e}")
    print(f"Diferencia máxima: {np.max(np.abs(E_paridad - E_unif)):.1e}")
