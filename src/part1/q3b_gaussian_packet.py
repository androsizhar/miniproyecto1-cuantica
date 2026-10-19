"""
Part 1, Question 3(b) — Gaussian wave packet in the symmetric finite well.

f(x) = A exp[-(x-c)^2/b^2] exp(i M v0 x / hbar),  A = (2/(pi b^2))^(1/4)
with b = 0.2, c = 0.4 and v0 = 8. Only the N = 11 bound states are used, so the
packet is reconstructed as g(x) = sum_n c_n phi_n(x) with sum |c_n|^2 ~ 0.91.

Outputs:
    figures/part1/gaussian_reconstruction.png
    figures/part1/gaussian_density.png
    figures/part1/gaussian_mean_position.png
    videos/gaussian_packet_symmetric_well.mp4
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

from q2_eigenstates import FIG, ROOT, HBAR, MASS, A_HALF, U0, MU_MAX, allowed_energies, eigenstates
from q3a_superposition import evolve, mean_position, plot_density, plot_mean_position

try:
    import imageio_ffmpeg
    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError:
    pass

trapezoid = getattr(np, "trapezoid", None) or np.trapz
VID = ROOT / "videos"
VID.mkdir(exist_ok=True)


def gaussian(x, b, c, v0, *, m=MASS, hbar=HBAR):
    """Normalized Gaussian packet of Eq. (5)."""
    A = (2 / (np.pi * b**2)) ** 0.25
    return A * np.exp(-((x - c) ** 2) / b**2) * np.exp(1j * m * v0 * x / hbar)


def expansion_coefficients(x, psi, f):
    """c_n = int phi_n^*(x) f(x) dx (Eq. 3, trapezoidal rule)."""
    return np.array([trapezoid(np.conj(p) * f, x) for p in psi])


def plot_reconstruction(x, f, g, path, *, a=A_HALF):
    """Initial packet |f|^2 versus its bound-state reconstruction |g|^2."""
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(x, np.abs(f) ** 2, lw=1, label=r"$|f(x)|^2$")
    ax.plot(x, np.abs(g) ** 2, lw=1, ls="--", label=r"$|g(x)|^2$")
    ax.axvline(-a, color="red", lw=0.6, ls="--")
    ax.axvline(a, color="red", lw=0.6, ls="--")
    ax.set_xlim(-2 * a, 2 * a)
    ax.set_xlabel(r"$x$")
    ax.set_title(r"Initial packet $f(x)$ and its reconstruction $g(x)$")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(path, dpi=200)
    plt.close(fig)


def make_video(x, t, phi, path, *, a=A_HALF, u0=U0, fps=30, title=""):
    """mp4 of |phi(x,t)|^2 together with the (scaled) potential."""
    rho = np.abs(phi) ** 2
    ymax = 1.1 * rho.max()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    pot = np.where(np.abs(x) <= a, 0, u0) / u0 * ymax * 0.9
    ax.fill_between(x, 0, pot, color="red", alpha=0.12, label=r"$U(x)$ (scaled)")
    line, = ax.plot(x, rho[0], lw=1.2, label=r"$|\phi(x,t)|^2$")
    re, = ax.plot(x, np.real(phi[0]), lw=0.6, alpha=0.6, label=r"Re$\,\phi$")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-np.sqrt(ymax), ymax)
    ax.set_xlabel(r"$x$")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
    label = ax.set_title("")

    def update(i):
        line.set_ydata(rho[i])
        re.set_ydata(np.real(phi[i]))
        label.set_text(f"{title}   t = {t[i]:.3f}")
        return line, re, label

    FuncAnimation(fig, update, frames=len(t), blit=False).save(path, writer=FFMpegWriter(fps=fps, bitrate=2400))
    plt.close(fig)


if __name__ == "__main__":
    x, psi, E = eigenstates(allowed_energies(MU_MAX))

    b, c, v0 = 0.2, 0.4, 8
    f = gaussian(x, b, c, v0)
    cn = expansion_coefficients(x, psi, f)
    g = sum(cn[n] * psi[n] for n in range(len(E)))
    P = np.sum(np.abs(cn) ** 2)
    print(f"A = {(2 / (np.pi * b**2)) ** 0.25:.4f}")
    print(f"sum |c_n|^2 = {P:.4f}")
    print(f"<E> of g = {np.sum(np.abs(cn) ** 2 * E) / P:.2f}")

    plot_reconstruction(x, f, g, FIG / "gaussian_reconstruction.png")
    t = np.linspace(0, 3, 1201)
    phi = evolve(psi, cn, E, t)
    plot_density(x, t, phi, r"Gaussian packet, symmetric well: $|\phi(x,t)|^2$",
                 FIG / "gaussian_density.png")
    plot_mean_position(t, mean_position(x, phi), r"Gaussian packet, symmetric well: $\langle x\rangle(t)$",
                       FIG / "gaussian_mean_position.png")
    print(f"Figures saved in {FIG.relative_to(ROOT)}/")

    print("Rendering video (about one minute)...")
    t_v = np.linspace(0, 3, 451)
    sel = np.abs(x) <= 2 * A_HALF
    make_video(x[sel], t_v, evolve(psi[:, sel], cn, E, t_v),
               VID / "gaussian_packet_symmetric_well.mp4", title="Symmetric well")
    print(f"Video saved in {VID.relative_to(ROOT)}/gaussian_packet_symmetric_well.mp4")
