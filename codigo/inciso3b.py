from pathlib import Path
import numpy as np
from scipy.optimize import brentq
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, FFMpegWriter

trapezoid = getattr(np, "trapezoid", None) or np.trapz

FIG = Path(__file__).resolve().parent.parent / "figuras"
FIG.mkdir(exist_ok=True)
VID = Path(__file__).resolve().parent.parent / "videos"
VID.mkdir(exist_ok=True)

try:
    import imageio_ffmpeg
    plt.rcParams["animation.ffmpeg_path"] = imageio_ffmpeg.get_ffmpeg_exe()
except ImportError:
    pass

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


def eigenestados(energias, M, hbar, U0, a, *, n=8, x_ext=3, steps=6001):
    '''
    Evalúa los primeros n eigenestados normalizados ordenados por energía.

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


def evolucion(psi, cn, E, t, *, hbar=1):
    '''
    phi(x,t) = sum_n c_n phi_n(x) exp(-i E_n t / hbar)  (Ec. 2).
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


def coeficientes(x, psi, f):
    '''
    c_n = int phi_n^*(x) f(x) dx  (Ec. 3) con la regla del trapecio.
    '''
    return np.array([trapezoid(np.conj(p)*f, x) for p in psi])


def gaussiano(x, b, c, v0, *, M=1/8*np.pi**2, hbar=1):
    '''
    Pulso Gaussiano de la Ec. (5) con A = (2/(pi b^2))^(1/4).
    '''
    A = (2/(np.pi*b**2))**0.25
    return A*np.exp(-(x - c)**2/b**2)*np.exp(1j*M*v0*x/hbar)


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


def grafica_fg(x, f, g, archivo, *, a=1):
    '''
    Compara el pulso inicial f(x) con su reconstrucción g(x) usando los N estados.
    '''
    fig, ax = plt.subplots(figsize=(7, 3.5))
    ax.plot(x, np.abs(f)**2, lw=1, label=r"$|f(x)|^2$")
    ax.plot(x, np.abs(g)**2, lw=1, ls="--", label=r"$|g(x)|^2$")
    ax.axvline(-a, color="red", lw=0.6, ls="--")
    ax.axvline(a, color="red", lw=0.6, ls="--")
    ax.set_xlim(-2*a, 2*a)
    ax.set_xlabel(r"$x$")
    ax.set_title(r"Pulso inicial $f(x)$ y su reconstrucción $g(x)$")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    fig.savefig(archivo, dpi=200)
    plt.show()


def video(x, t, phi, archivo, *, a=1, U0=110, fps=30, titulo=""):
    '''
    Video mp4 de la evolución de |phi(x,t)|^2 junto con el potencial.
    '''
    rho = np.abs(phi)**2
    ymax = 1.1*rho.max()
    fig, ax = plt.subplots(figsize=(8, 4.5))
    pot = np.where(np.abs(x) <= a, 0, U0)/U0*ymax*0.9
    ax.fill_between(x, 0, pot, color="red", alpha=0.12, label=r"$U(x)$ (escalado)")
    linea, = ax.plot(x, rho[0], lw=1.2, label=r"$|\phi(x,t)|^2$")
    re, = ax.plot(x, np.real(phi[0]), lw=0.6, alpha=0.6, label=r"Re$\,\phi$")
    ax.set_xlim(x[0], x[-1])
    ax.set_ylim(-np.sqrt(ymax), ymax)
    ax.set_xlabel(r"$x$")
    ax.legend(loc="upper right")
    ax.grid(alpha=0.3)
    texto = ax.set_title("")

    def update(i):
        linea.set_ydata(rho[i])
        re.set_ydata(np.real(phi[i]))
        texto.set_text(f"{titulo}   t = {t[i]:.3f}")
        return linea, re, texto

    anim = FuncAnimation(fig, update, frames=len(t), blit=False)
    anim.save(archivo, writer=FFMpegWriter(fps=fps, bitrate=2400))
    plt.close(fig)


if __name__ == "__main__":
    energias = energias_permisibles(mumax)
    N = len(energias[0]) + len(energias[1])
    x, psi, E = eigenestados(energias, M, hbar, U0, a, n=N)

    b, c, v0 = 0.2, 0.4, 8
    f = gaussiano(x, b, c, v0)
    cn = coeficientes(x, psi, f)
    g = sum(cn[n]*psi[n] for n in range(N))
    P = np.sum(np.abs(cn)**2)

    print(f"A = {(2/(np.pi*b**2))**0.25:.4f}")
    print(f"Norma de f: {trapezoid(np.abs(f)**2, x):.6f}")
    print("|c_n|^2 =", np.round(np.abs(cn)**2, 4))
    print(f"Suma |c_n|^2 = {P:.4f}")
    print(f"<E> de g = {np.sum(np.abs(cn)**2*E)/P:.2f}")

    grafica_fg(x, f, g, FIG / "p1_3b_fg.png")

    t = np.linspace(0, 3, 1201)
    phi = evolucion(psi, cn, E, t)
    mapa_densidad(x, t, phi, r"Inciso 3(b): $|\phi(x,t)|^2$ del pulso Gaussiano", FIG / "p1_3b_mapa.png")
    grafica_x(t, valor_esperado_x(x, phi), r"Inciso 3(b): $\langle x\rangle(t)$ del pulso Gaussiano", FIG / "p1_3b_x.png")

    print("Generando video (tarda alrededor de un minuto)...")
    t_v = np.linspace(0, 3, 451)
    sel = np.abs(x) <= 2*a
    video(x[sel], t_v, evolucion(psi[:, sel], cn, E, t_v), VID / "pulso_pozo_simetrico.mp4",
          titulo="Pozo simétrico")
    print("Video guardado en videos/pulso_pozo_simetrico.mp4")
