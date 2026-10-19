"""
Part 1, Question 3(a) — Time evolution of a four-state superposition.

phi(x,0) = A [phi_1 + 0.7 phi_2 - 1.2i phi_3 + (2 - i) phi_4],  A = 1/sqrt(7.93)

The state is evolved with the eigenstate expansion
    phi(x,t) = sum_n c_n phi_n(x) exp(-i E_n t / hbar).

Outputs (figures/part1/):
    superposition_density.png, superposition_mean_position.png
"""
import numpy as np
import matplotlib.pyplot as plt

from q2_eigenstates import FIG, ROOT, HBAR, A_HALF, MU_MAX, allowed_energies, eigenstates

trapezoid = getattr(np, "trapezoid", None) or np.trapz


def evolve(psi, cn, E, t, *, hbar=HBAR):
    """phi(x,t) for every time in t; returns a complex array of shape (len(t), len(x))."""
    phi = np.zeros((len(t), psi.shape[1]), dtype=complex)
    for n in range(len(E)):
        if cn[n] != 0:
            phi += cn[n] * np.exp(-1j * E[n] * t / hbar)[:, None] * psi[n][None, :]
    return phi


def mean_position(x, phi):
    """<x>(t) = int x |phi|^2 dx / int |phi|^2 dx."""
    rho = np.abs(phi) ** 2
    return trapezoid(x * rho, x, axis=1) / trapezoid(rho, x, axis=1)


def plot_density(x, t, phi, title, path, *, a=A_HALF):
    """2D map of |phi(x,t)|^2."""
    fig, ax = plt.subplots(figsize=(7, 4.5))
    im = ax.pcolormesh(x, t, np.abs(phi) ** 2, shading="auto", cmap="viridis")
    fig.colorbar(im, ax=ax, label=r"$|\phi(x,t)|^2$")
    ax.axvline(-a, color="red", lw=0.8, ls="--")
    ax.axvline(a, color="red", lw=0.8, ls="--")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$t$")
    ax.set_title(title)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_mean_position(t, xm, title, path, *, a=A_HALF):
    """<x> as a function of time, with the walls as reference."""
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(t, xm, lw=1)
    ax.axhline(-a, color="red", lw=0.6, ls="--")
    ax.axhline(a, color="red", lw=0.6, ls="--")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle x\rangle$")
    ax.set_title(title)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


if __name__ == "__main__":
    x, psi, E = eigenstates(allowed_energies(MU_MAX))

    coeffs = np.array([1, 0.7, -1.2j, 2 - 1j])
    A = 1 / np.sqrt(np.sum(np.abs(coeffs) ** 2))
    cn = np.zeros(len(E), dtype=complex)
    cn[:4] = A * coeffs

    t = np.linspace(0, 10, 1201)
    phi = evolve(psi, cn, E, t)
    print(f"A = {A:.4f}")
    print(f"Norm at t=0: {trapezoid(np.abs(phi[0]) ** 2, x):.10f}")

    plot_density(x, t, phi, r"Four-state superposition: $|\phi(x,t)|^2$",
                 FIG / "superposition_density.png")
    plot_mean_position(t, mean_position(x, phi), r"Four-state superposition: $\langle x\rangle(t)$",
                       FIG / "superposition_mean_position.png")
    print(f"Figures saved in {FIG.relative_to(ROOT)}/")
