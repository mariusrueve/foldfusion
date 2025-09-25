# Latexmk configuration for the thesis build

# Produce PDF via pdflatex
$pdf_mode = 1;

# Use pdflatex with strict, non-interactive options
$pdflatex = 'pdflatex -interaction=nonstopmode -halt-on-error %O %S';

# Always use biber for bibliography processing
# Use biber explicitly and tell latexmk to use it
$biber = 'biber %O %B';
$bibtex_use = 2;  # 2 = use biber (not bibtex)
