"""
Part 2, Question 6 — Asymmetric finite well (Python version).

Potential: U1 for x < 0, 0 for 0 <= x < L, U2 = c*U1 for x >= L.
Parameters: hbar = 1, L = 2, M = pi^2/8, U1 = 85, c = 2.

Same method and parameters as the original Julia code (q6_asymmetric_well.jl):
bound-state energies from the characteristic equation, normalized eigenstates,
and the evolution of a Gaussian packet (b = 0.2, centre 0.4, v0 = 0, t in [0, 5]).

Outputs:
    figures/part2/asymmetric_well_potential.png
    figures/part2/energy_roots.png
    figures/part2/eigenstates.png
    figures/part2/expansion_coefficients.png
    figures/part2/gaussian_density.png
    figures/part2/gaussian_mean_position.png
    videos/gaussian_packet_asymmetric_well.mp4
"""
from pathlib import Path

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter
from scipy.optimize import brentq

try:
    import imageio_ffmpeg
    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError:
    pass

ROOT = Path(__file__).resolve().parents[2]
FIG = ROOT / "figures" / "part2"
VID = ROOT / "videos"
FIG.mkdir(parents=True, exist_ok=True)
VID.mkdir(exist_ok=True)

trapezoid = getattr(np, "trapezoid", None) or np.trapz

HBAR = 1.0
L = 2.0
MASS = np.pi**2 / 8
U1 = 85.0
C_ASYM = 2.0
U2 = C_ASYM * U1


def bound_energies(u1=U1, u2=U2, *, length=L, m=MASS, hbar=HBAR):
    """
    Solve the characteristic equation written as a phase condition,
        w L - arctan(kappa/w)' ... = (n - 1) pi,
    i.e.  g(E) = (L/hbar) sqrt(2 M E) - arctan(sqrt((U1-E)/E)) - arctan(sqrt((U2-E)/E)) = (n-1) pi.
    g(E) is increasing on (0, U1), so there is one root per n while g(U1) > (n-1) pi.
    """
    g = lambda e: (length / hbar) * np.sqrt(2 * m * e) - np.arctan(np.sqrt((u1 - e) / e)) - np.arctan(np.sqrt((u2 - e) / e))
    e_max = u1 * (1 - 1e-12)
    roots = []
    n = 1
    while g(e_max) > (n - 1) * np.pi:
        roots.append(brentq(lambda e: g(e) - (n - 1) * np.pi, 1e-9, e_max, xtol=1e-13))
        n += 1
    return np.array(roots)


def eigenstate(x, E, u1=U1, u2=U2, *, length=L, m=MASS, hbar=HBAR):
    """
    Normalized eigenstate:
        x < 0:       C exp(kappa x)
        0 <= x < L:  C [cos(w x) + (kappa/w) sin(w x)]
        x >= L:      D exp(-gamma (x - L))
    with the interior amplitude sqrt(2 kappa gamma / (kappa gamma L + kappa + gamma)).
    """
    w = np.sqrt(2 * m * E) / hbar
    kap = np.sqrt(2 * m * (u1 - E)) / hbar
    gam = np.sqrt(2 * m * (u2 - E)) / hbar
    amp = np.sqrt(2 * kap * gam / (kap * gam * length + kap + gam))
    C = (w / np.sqrt(w**2 + kap**2)) * amp
    D = C * (np.cos(w * length) + (kap / w) * np.sin(w * length))
    inside = C * (np.cos(w * x) + (kap / w) * np.sin(w * x))
    return np.where(x < 0, C * np.exp(kap * np.minimum(x, 0)),
                    np.where(x < length, inside, D * np.exp(-gam * (np.maximum(x, length) - length))))


def potential(x, u1=U1, u2=U2, *, length=L):
    return np.where(x < 0, u1, np.where(x < length, 0.0, u2))


def plot_potential_diagram(path):
    """Schematic of the asymmetric well (Fig. 1 of the assignment)."""
    fig, ax = plt.subplots(figsize=(6, 3.2))
    x = np.linspace(-1.2, L + 1.2, 1000)
    ax.plot(x, potential(x, 1.0, 2.0), color="firebrick", lw=2)
    ax.axhline(0.7, color="gray", lw=0.8, ls=":")
    ax.text(-1.15, 0.75, r"$E$", fontsize=12)
    ax.text(-1.15, 1.07, r"$U_1$", fontsize=12)
    ax.text(L + 0.6, 2.07, r"$U_2 = c\,U_1$", fontsize=12)
    ax.set_xticks([0, L], ["0", r"$L$"])
    ax.set_yticks([])
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$U(x)$")
    ax.set_ylim(-0.2, 2.5)
    for side in ("top", "right"):
        ax.spines[side].set_visible(False)
    ax.set_title("Asymmetric finite well")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_energy_roots(E_n, path):
    """Graphical solution: tan(wL) against the right-hand side of the characteristic equation."""
    e = np.linspace(1e-3, U1 - 1e-3, 6000)
    lhs = np.tan((L / HBAR) * np.sqrt(2 * MASS * e))
    rhs = np.sqrt(e) * (np.sqrt(U1 - e) + np.sqrt(U2 - e)) / (e - np.sqrt((U1 - e) * (U2 - e)))
    lhs[np.abs(np.diff(lhs, append=lhs[-1])) > 160] = np.nan
    rhs[np.abs(np.diff(rhs, append=rhs[-1])) > 20] = np.nan
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.plot(e, lhs, lw=1, label=r"$\tan(\omega L)$")
    ax.plot(e, rhs, lw=1, label=r"$\frac{\omega(\kappa+\gamma)}{\omega^2-\kappa\gamma}$")
    y_n = np.tan((L / HBAR) * np.sqrt(2 * MASS * E_n))
    ax.scatter(E_n, y_n, color="red", s=20, zorder=3, label="Eigenenergies")
    for n, (en, yn) in enumerate(zip(E_n, y_n), start=1):
        dy = 14 if n % 2 else -18
        ax.annotate(f"$E_{n}$", (en, yn), textcoords="offset points", xytext=(0, dy),
                    ha="center", fontsize=8, arrowprops=dict(arrowstyle="-", lw=0.4, color="gray"))
    ax.set_ylim(-160, 160)
    ax.set_xlabel(r"Energy $E$")
    ax.set_title(f"Graphical solution of the characteristic equation ($c$ = {C_ASYM:g})")
    ax.legend(loc="lower left", fontsize=8)
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_eigenstates(x, psi, E_n, path, *, scale=3.0):
    """Eigenstates drawn at the height of their energy inside the potential."""
    fig, ax = plt.subplots(figsize=(8.5, 6))
    ax.plot(x, potential(x), color="black", lw=1.2, label=r"$U(x)$")
    for n, (p, e) in enumerate(zip(psi, E_n), start=1):
        ax.hlines(e, x[0], x[-1], color="gray", lw=0.3, ls=":")
        ax.plot(x, e + scale * p, lw=1, label=f"n = {n}")
    ax.set_ylim(-5, U1 + 15)
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"Energy")
    ax.set_title(r"Eigenstates of the asymmetric well ($\varphi_n$ shifted by $E_n$)")
    ax.legend(fontsize=7, loc="upper left", ncol=2)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_coefficients(cn, path):
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.bar(np.arange(1, len(cn) + 1), np.abs(cn) ** 2)
    ax.set_xticks(range(1, len(cn) + 1))
    ax.set_xlabel(r"$n$")
    ax.set_ylabel(r"$|c_n|^2$")
    ax.set_title(rf"Expansion weights, $\sum|c_n|^2$ = {np.sum(np.abs(cn) ** 2):.4f}")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def evolve(psi, cn, E, t):
    phi = np.zeros((len(t), psi.shape[1]), dtype=complex)
    for n in range(len(E)):
        phi += cn[n] * np.exp(-1j * E[n] * t / HBAR)[:, None] * psi[n][None, :]
    return phi


def plot_density(x, t, phi, path):
    fig, ax = plt.subplots(figsize=(7, 4.5))
    im = ax.pcolormesh(x, t, np.abs(phi) ** 2, shading="auto", cmap="viridis")
    fig.colorbar(im, ax=ax, label=r"$|g(x,t)|^2$")
    for xw in (0, L):
        ax.axvline(xw, color="red", lw=0.8, ls="--")
    ax.set_xlabel(r"$x$")
    ax.set_ylabel(r"$t$")
    ax.set_title(r"Gaussian packet, asymmetric well: $|g(x,t)|^2$")
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def plot_mean_position(t, xm, path):
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(t, xm, lw=1)
    for xw in (0, L):
        ax.axhline(xw, color="red", lw=0.6, ls="--")
    ax.set_xlabel(r"$t$")
    ax.set_ylabel(r"$\langle x\rangle$")
    ax.set_title(r"Gaussian packet, asymmetric well: $\langle x\rangle(t)$")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def make_video(x, t, phi, path, *, fps=30):
    rho = np.abs(phi) ** 2
    ymax = 1.1 * rho.max()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    ax.fill_between(x, 0, potential(x) / U2 * ymax * 0.9, color="red", alpha=0.12, label=r"$U(x)$ (scaled)")
    line, = ax.plot(x, rho[0], lw=1.2, label=r"$|g(x,t)|^2$")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(0, ymax)
    ax.set_xlabel(r"$x$")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
    label = ax.set_title("")

    def update(i):
        line.set_ydata(rho[i])
        label.set_text(f"Asymmetric well   t = {t[i]:.3f}")
        return line, label

    FuncAnimation(fig, update, frames=len(t), blit=False).save(path, writer=FFMpegWriter(fps=fps, bitrate=2400))
    plt.close(fig)


if __name__ == "__main__":
    E_n = bound_energies()
    print(f"c = {C_ASYM:g}:  N = {len(E_n)} bound states")
    for n, e in enumerate(E_n, start=1):
        print(f"  E_{n} = {e:.4f}")

    e_sym = bound_energies(110.0, 110.0)
    print(f"Check c = 1 (U1 = U2 = 110): N = {len(e_sym)}, E_1 = {e_sym[0]:.4f}, E_11 = {e_sym[-1]:.4f}")

    x = np.linspace(-1.5, L + 1.5, 4001)
    psi = np.array([eigenstate(x, e) for e in E_n])
    S = trapezoid(psi[:, None, :] * psi[None, :, :], x, axis=-1)
    print(f"max |S - I| = {np.max(np.abs(S - np.eye(len(E_n)))):.1e}")

    b, centre, v0 = 0.2, 0.4, 0.0
    A = (2 / (np.pi * b**2)) ** 0.25
    f = A * np.exp(-((x - centre) ** 2) / b**2) * np.exp(1j * MASS * v0 * x / HBAR)
    cn = np.array([trapezoid(p * f, x) for p in psi])
    print("c_n =", np.round(cn.real, 4))
    print(f"sum |c_n|^2 = {np.sum(np.abs(cn) ** 2):.4f}")

    plot_potential_diagram(FIG / "asymmetric_well_potential.png")
    plot_energy_roots(E_n, FIG / "energy_roots.png")
    plot_eigenstates(x, psi, E_n, FIG / "eigenstates.png")
    plot_coefficients(cn, FIG / "expansion_coefficients.png")

    t = np.linspace(0, 5, 1001)
    phi = evolve(psi, cn, E_n, t)
    rho = np.abs(phi) ** 2
    xm = trapezoid(x * rho, x, axis=1) / trapezoid(rho, x, axis=1)
    plot_density(x, t, phi, FIG / "gaussian_density.png")
    plot_mean_position(t, xm, FIG / "gaussian_mean_position.png")
    print(f"Figures saved in {FIG.relative_to(ROOT)}/")

    print("Rendering video...")
    t_v = np.linspace(0, 5, 301)
    sel = (x >= -0.5) & (x <= L + 0.5)
    make_video(x[sel], t_v, evolve(psi[:, sel], cn, E_n, t_v), VID / "gaussian_packet_asymmetric_well.mp4")
    print(f"Video saved in {VID.relative_to(ROOT)}/gaussian_packet_asymmetric_well.mp4")
