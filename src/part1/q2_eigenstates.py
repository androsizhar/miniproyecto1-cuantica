"""
Part 1, Question 2 — Eigenenergies and eigenstates of the symmetric finite well.

Finds every bound-state energy of the well |x| < a with walls of height U0,
plots the energy levels, the first eight normalized eigenstates (analytic
amplitudes A_n, B_n) and the overlap matrix S_ij = <phi_i|phi_j>.

Outputs (figures/part1/):
    energy_levels.png, eigenstates.png, overlap_matrix.png
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import brentq

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures" / "part1"
FIG.mkdir(parents=True, exist_ok=True)

trapezoid = getattr(np, "trapezoid", None) or np.trapz

# Units: hbar = 1, a = 1, M = pi^2/8  ->  V0 = 1 (ground energy of the infinite well of width 2a)
HBAR = 1.0
A_HALF = 1.0
MASS = np.pi**2 / 8
V0 = (HBAR * np.pi) ** 2 / (2 * MASS * (2 * A_HALF) ** 2)
U0 = 110.0
MU_MAX = np.sqrt(U0 / V0)


def allowed_energies(mu_max, *, eps=1e-12, a=A_HALF, hbar=HBAR, m=MASS):
    """
    Solve the even (Eq. 14) and odd (Eq. 17) characteristic equations in terms of
    mu = (2a/pi) k. Even roots lie in (2n, 2n+1) and odd roots in (2n+1, 2n+2).

    Returns:
        (even_energies, odd_energies) as arrays.
    """
    f_even = lambda x: x * np.sin(np.pi * x / 2) - np.sqrt(mu_max**2 - x**2) * np.cos(np.pi * x / 2)
    f_odd = lambda x: np.sqrt(mu_max**2 - x**2) * np.sin(np.pi * x / 2) + x * np.cos(np.pi * x / 2)
    even_roots, odd_roots = [], []
    n = 0
    while 2 * n < mu_max:
        low, high = 2 * n + eps, min(2 * n + 1 - eps, mu_max)
        if f_even(low) * f_even(high) < 0:
            even_roots.append(brentq(f_even, low, high))
        low += 1
        high = min(high + 1, mu_max)
        if low > mu_max:
            break
        if f_odd(low) * f_odd(high) < 0:
            odd_roots.append(brentq(f_odd, low, high))
        n += 1
    e1 = (hbar * np.pi) ** 2 / (2 * m * (2 * a) ** 2)
    return e1 * np.array(even_roots) ** 2, e1 * np.array(odd_roots) ** 2


def eigenstates(energies, *, n=None, a=A_HALF, u0=U0, hbar=HBAR, m=MASS, x_ext=3.0, steps=6001):
    """
    Evaluate the normalized eigenstates (analytic amplitudes of Eqs. 18-19),
    ordered by increasing energy.

    Returns:
        x (steps,), psi (n, steps), E (n,)
    """
    e_even, e_odd = energies
    total = len(e_even) + len(e_odd)
    n = total if n is None else n
    x = np.linspace(-x_ext * a, x_ext * a, steps)
    inside = np.abs(x) <= a
    psi = np.empty((n, steps))
    E = np.empty(n)
    for j in range(n):
        i = j // 2
        E[j] = e_even[i] if j % 2 == 0 else e_odd[i]
        k = np.sqrt(2 * m * E[j]) / hbar
        q = np.sqrt(2 * m * (u0 - E[j])) / hbar
        amp = np.sqrt(q / (a * q + 1))
        if j % 2 == 0:
            b = amp * np.cos(k * a) * np.exp(q * a)
            psi[j] = np.where(inside, amp * np.cos(k * x), b * np.exp(-q * np.abs(x)))
        else:
            b = amp * np.sin(k * a) * np.exp(q * a)
            psi[j] = np.where(inside, amp * np.sin(k * x), np.sign(x) * b * np.exp(-q * np.abs(x)))
    return x, psi, E


def plot_energy_levels(E, path, *, a=A_HALF, u0=U0):
    """Energy-level diagram with leader lines so that close levels stay readable."""
    fig, ax = plt.subplots(figsize=(7, 5))
    x = np.linspace(-2 * a, 2 * a, 2001)
    ax.plot(x, np.where(np.abs(x) <= a, 0.0, u0), lw=1)
    y_prev = -np.inf
    for n, e in enumerate(E, start=1):
        ax.hlines(e, -a, a, lw=0.6)
        y = max(e, y_prev + 7.5)
        ax.plot([a, a + 0.12, a + 0.22], [e, e, y], color="gray", lw=0.5)
        ax.text(a + 0.25, y, f"n={n}  ($E$={e:.2f})", va="center", fontsize=8)
        y_prev = y
    ax.set_xlim(-2 * a, 2 * a)
    ax.set_ylim(-0.1 * u0, 1.3 * u0)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$U(x)$")
    ax.set_title("Allowed energies of the symmetric finite well")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_eigenstates(x, psi, path, *, a=A_HALF, n_plot=8):
    """First n_plot eigenstates in a 4x2 grid."""
    fig, axes = plt.subplots(4, 2, figsize=(10, 12), sharex=True)
    sel = np.abs(x) <= 2 * a
    for j, ax in zip(range(n_plot), axes.flat):
        ax.plot(x[sel], psi[j, sel], lw=1.2)
        ax.axvline(-a, color="red", lw=0.4, ls="--")
        ax.axvline(a, color="red", lw=0.4, ls="--")
        ax.set_title(f"n={j + 1}")
        ax.grid(alpha=0.3)
    for ax in axes[-1]:
        ax.set_xlabel(r"$x$")
    for ax in axes[:, 0]:
        ax.set_ylabel(r"$\varphi_n(x)$")
    fig.suptitle(f"First {n_plot} eigenstates of the finite well", fontsize=14)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_overlap_matrix(x, psi, path):
    """Overlap matrix S_ij = int phi_i phi_j dx (trapezoidal rule)."""
    S = trapezoid(psi[:, None, :] * psi[None, :, :], x, axis=-1)
    fig, ax = plt.subplots()
    im = ax.imshow(S, cmap="Greys", vmin=-1, vmax=1)
    fig.colorbar(im, ax=ax, label=r"$\langle\varphi_i|\varphi_j\rangle$")
    ax.set_xticks(range(len(S)), labels=range(1, len(S) + 1))
    ax.set_yticks(range(len(S)), labels=range(1, len(S) + 1))
    ax.set_xlabel("j")
    ax.set_ylabel("i")
    ax.set_title("Overlap matrix $S$")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)
    return S


if __name__ == "__main__":
    energies = allowed_energies(MU_MAX)
    x, psi, E = eigenstates(energies)
    print(f"N = {len(E)} bound states")
    for n, e in enumerate(E, start=1):
        print(f"  E_{n:<2d} = {e:9.4f}")

    plot_energy_levels(E, FIG / "energy_levels.png")
    plot_eigenstates(x, psi, FIG / "eigenstates.png")
    S = plot_overlap_matrix(x, psi, FIG / "overlap_matrix.png")
    print(f"max |S - I| = {np.max(np.abs(S - np.eye(len(S)))):.1e}")
    print(f"Figures saved in {FIG.relative_to(ROOT)}/")
