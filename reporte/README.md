# Mini-proyecto de Cuántica

Informe LaTeX del pozo finito simétrico (Parte 1) y asimétrico (Parte 2).

## Estructura

- `miniproyecto1.tex`: documento principal; reúne las secciones y el apéndice.
- `mini1/parte1.tex` y `mini1/parte2.tex`: desarrollo de las partes del informe.
- `mini1/appendix.tex`: apéndice con material complementario.
- `figuras/`, `mini1/` y la raíz: imágenes referenciadas por el documento.

Se mantienen las rutas actuales para no cambiar el contenido del informe ni romper referencias a figuras.

## Compilar

Desde la raíz del proyecto, ejecutar:

```sh
latexmk -pdf -interaction=nonstopmode -halt-on-error miniproyecto1.tex
```

Requiere una distribución LaTeX con los paquetes usados por el documento y `latexmk`. La compilación genera archivos auxiliares junto a las fuentes; `.gitignore` los excluye.

## Colaborar

- Edita la sección correspondiente en `mini1/` y conserva el documento principal como punto de entrada.
- Compila desde la raíz: varias rutas de imágenes se resuelven desde ese directorio.
- Si agregas o cambias una figura, comprueba que su nombre y ruta coincidan exactamente con `\includegraphics` en los `.tex`.
- Antes de compartir cambios, revisa `git status` y asegúrate de incluir las imágenes necesarias en el repositorio; los archivos que Git no rastrea no llegarán a las demás personas.
- Mantén los productos de compilación fuera de los cambios y describe en tus commits qué sección o figuras modificaste.

## Flujo sugerido

1. Crea una rama para cada cambio o sección.
2. Edita y compila el documento desde la raíz.
3. Revisa el PDF y `git status` antes de compartir la rama.
4. Abre una solicitud de revisión con un resumen breve de los cambios.