"""
Part 1, Question 4(b) — Unified characteristic equation.

The round-trip phase condition 4ka + 2*gamma = 2*pi*m gives a single equation
for even and odd states:
    mu_n = n - (2/pi) arcsin(mu_n / mu_max),   n = 1, 2, ...
This script solves it and compares with the parity equations (14) and (17).
"""
import numpy as np
from scipy.optimize import brentq

from q2_eigenstates import V0, MU_MAX, allowed_energies


def unified_energies(mu_max, *, v0=V0):
    """Roots of F(mu) = mu + (2/pi) arcsin(mu/mu_max) - n; one root per n while F(mu_max) > n."""
    F = lambda mu, n: mu + (2 / np.pi) * np.arcsin(mu / mu_max) - n
    mus = []
    n = 1
    while F(mu_max, n) > 0:
        mus.append(brentq(F, 0, mu_max, args=(n,), xtol=1e-14))
        n += 1
    return v0 * np.array(mus) ** 2


if __name__ == "__main__":
    e_even, e_odd = allowed_energies(MU_MAX)
    e_parity = np.sort(np.concatenate((e_even, e_odd)))
    e_unified = unified_energies(MU_MAX)

    print(f"N (parity) = {len(e_parity)},  N (unified) = {len(e_unified)},  ceil(mu_max) = {int(np.ceil(MU_MAX))}")
    print(f"{'n':>3} {'Eqs. 14 & 17':>14} {'Unified':>12} {'Difference':>12}")
    for n, (e1, e2) in enumerate(zip(e_parity, e_unified), start=1):
        print(f"{n:>3} {e1:>14.4f} {e2:>12.4f} {abs(e1 - e2):>12.1e}")
    print(f"Max difference: {np.max(np.abs(e_parity - e_unified)):.1e}")
