# Mini-proyecto 1 — Mecánica Cuántica

**Sistemas simples en 1D: pozo finito de potencial simétrico y asimétrico**

Tecnológico de Monterrey · Fundamentación de la electrodinámica (Gpo 601)
Profesor: Julio César Gutiérrez Vega

**Equipo**
- Luis Daniel González Alcocer — A01286903
- Andrés Saul Pérez Martínez — A01199376
- Andros Izhar Alvarado Sánchez — A00842107

**Entrega:** lunes 19 de octubre de 2026, antes de las 08:59 AM.

---

## Estructura del repositorio

```
miniproyecto1-cuantica/
├── README.md
├── requirements.txt
├── .gitignore
├── codigo/       ← un script de Python por inciso
├── figuras/      ← las gráficas que generan los scripts (.png)
├── videos/       ← animaciones de la evolución del pulso (.mp4)
└── reporte/      ← proyecto LaTeX del reporte (miniproyecto1.tex)
```

## Cómo correr el código

1. Instalar Python 3.10 o más reciente.
2. Instalar las librerías (una sola vez):
   ```bash
   pip install -r requirements.txt
   ```
3. Para los videos se necesita además **ffmpeg**:
   - Mac: `brew install ffmpeg`
   - Windows: `winget install ffmpeg`
4. Abrir la carpeta del repositorio en VS Code y correr cualquier script de `codigo/` con ▶ (Run Python File).
   También se puede desde la terminal:
   ```bash
   python codigo/inciso2.py
   ```

Cada script es independiente y guarda sus figuras en `figuras/` automáticamente, sin importar desde qué carpeta se corra.

**Parámetros numéricos** (los mismos en todos los scripts): ħ = 1, M = π²/8, a = 1 (L = 2), lo que da V₀ = 1.

## Avance

### Parte 1: Pozo finito simétrico
| Inciso | Script | Estado |
|---|---|---|
| 1. Amplitudes Aₙ, Bₙ (analítico) | — | ✅ En el reporte |
| 2. Eigenenergías, eigenestados, ortonormalidad | `inciso2.py` | ✅ |
| 3(a). Evolución de la superposición de 4 estados | `inciso3a.py` | ✅ |
| 3(b). Pulso Gaussiano + video | `inciso3b.py` | ⬜ |
| 4(a). Ecuación característica unificada (analítico) | — | ⬜ |
| 4(b). Verificación numérica de la ecuación unificada | `inciso4b.py` | ⬜ |

### Parte 2: Pozo finito asimétrico
| Inciso | Script | Estado |
|---|---|---|
| 1. Forma de las eigenfunciones | — | ✅ En el reporte |
| 2. Ecuación característica | — | ✅ En el reporte |
| 3. Número de estados N | — | ⬜ |
| 4. Condición de corte | — | ⬜ |
| 5. Amplitudes (normalización) | — | ⬜ Corregir: hay que integrar \|φ\|² |
| 6(a). Energías con c = 2 | `p2_inciso6a.py` | ⬜ |
| 6(b). Eigenestados | `p2_inciso6b.py` | ⬜ |
| 6(c). Pulso Gaussiano + video | `p2_inciso6c.py` | ⬜ |

### Discusión final
Individual: cada integrante escribe la suya.

## Cómo trabajar en equipo

1. **Antes de empezar**, traer los cambios de los demás:
   ```bash
   git pull
   ```
2. **Repartirse los incisos** para no editar el mismo archivo al mismo tiempo. Cada inciso tiene su propio script, así que casi no hay conflictos.
3. **Al terminar algo**, subirlo:
   ```bash
   git add .
   git commit -m "Inciso 3(b): pulso gaussiano y video"
   git push
   ```
4. Si `git push` falla porque alguien subió algo antes, hacer `git pull` y luego `git push` otra vez.
5. Para el reporte, cuidado con editar `reporte/mini1/parte1.tex` o `parte2.tex` al mismo tiempo que otra persona. Avisen en el chat del equipo quién está editando qué.

## Notas

- El repositorio debe ser **privado** hasta después de la fecha de entrega.
- Los videos pesan poco (unos 4 MB). Si alguno pasa de 100 MB, GitHub no lo acepta: súbanlo a YouTube o Drive y pongan el link aquí y en el reporte.
