from setuptools import setup, find_packages

setup(
    name="foldfusion",
    version="0.1.0",
    packages=find_packages(),
    description="FoldFusion: Ligand Transplantation and Optimization in Protein Models",
    author="Marius Rueve",
    author_email="",
    install_requires=["setuptools", "numpy", "pandas", "requests"],
    python_requires=">=3.12",
)
