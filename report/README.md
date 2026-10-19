# Report

LaTeX sources of the written report (in Spanish). The compiled version is [`report.pdf`](report.pdf).

```
report/
├── main.tex              preamble, title page and document structure
├── sections/
│   ├── part1.tex         Part 1: symmetric finite well
│   ├── part2.tex         Part 2: asymmetric finite well
│   └── appendix.tex      source code listings
├── reflections/          individual final reflections
└── figures/              images used in the report (p1_* Part 1, p2_* Part 2)
```

The Part 1 figures were produced by the scripts in [`src/part1/`](../src/part1/) (with the plot labels in Spanish) and the Part 2 figures by [`src/part2/q6_asymmetric_well.jl`](../src/part2/q6_asymmetric_well.jl). `p2_asymmetric_well_assignment.png` is the potential diagram from the assignment statement.

## Building

From this folder:

```sh
latexmk -pdf main.tex
```
