import os

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

HERE = os.path.dirname(os.path.abspath(__file__))
INCLUDE = os.path.join(HERE, "..", "include")

ext_modules = [
    Pybind11Extension(
        "opal._opal",
        ["bindings.cpp"],
        include_dirs=[INCLUDE],
        cxx_std=17,
    ),
]

setup(
    name="opal-pricing",
    version="0.1.0",
    description="Institutional option pricing library (C++ core, Python API)",
    long_description=(
        "Opal prices and risk-manages equity and interest rate options: "
        "Black-Scholes-Merton, Black-76, Bachelier, Heston, SABR, Hull-White; "
        "binomial/trinomial trees, Crank-Nicolson PDE and Monte Carlo engines; "
        "vanilla, digital, barrier, Asian and lookback payoffs; caps/floors "
        "and swaptions."
    ),
    packages=["opal"],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    python_requires=">=3.8",
    zip_safe=False,
)
