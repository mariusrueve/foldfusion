# Thesis Build Guide

This directory contains the LaTeX sources for the accompanying master’s thesis. All writing must use American English and a professional, direct style. See `LANGUAGE.md` for rules.

## Build

- Preferred: `make thesis` (uses `latexmk -pdf`; biber configured in `.latexmkrc`)
- Clean aux files: `make thesis-clean` (uses `latexmk -C`)
- Artifact: `thesis/main.pdf`

Direct command (equivalent):

```
cd thesis
latexmk -pdf -interaction=nonstopmode -halt-on-error main.tex
```

## Prerequisites

Install a TeX distribution with: `latexmk`, `pdflatex`, `biber`, and common packages (e.g., `biblatex`, `csquotes`, `microtype`, `siunitx`, `hyperref`, `cleveref`, `geometry`, `graphicx`).

## Common Pitfalls

- Biber missing: install `biber` (TeX Live/MiKTeX or system package manager).
- “Please (re)run Biber”: run `make thesis` again; ensure the `.bib` file exists and is referenced correctly.
- Missing LaTeX package: install via TeX package manager (e.g., `tlmgr install <pkg>`).
- Stale aux files or odd references: `make thesis-clean` then rebuild.
- Figure not found: verify relative paths and file extensions; update `\graphicspath{}` if used.
- Encoding issues: ensure sources and bibliography are UTF‑8; avoid non-ASCII in labels/keys.

Keep content modular and consistent; prefer clarity to complexity.

## Acronyms & Glossary

- Definitions live in `thesis/acronyms.tex` via `\newacronym{label}{SHORT}{Long form}`.
- Use in text: `\gls{label}` (first use: Long form (SHORT); later uses: SHORT). Capitalized: `\Gls{}`.
- Print list: automatically inserted after the Table of Contents.
- Add entries as needed (e.g., `md`, `rmsd`, `pdb`, `af`).
- Build: no extra steps — `latexmk` with `\makenoidxglossaries` handles it.
