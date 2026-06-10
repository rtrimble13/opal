import os
import shutil

from pybind11.setup_helpers import Pybind11Extension, build_ext
from setuptools import setup

HERE = os.path.dirname(os.path.abspath(__file__))

# The C++ headers live at the repo root (../include). Bundle a copy inside
# the package directory so sdists and cibuildwheel builds are self-contained;
# refresh the copy whenever the repo headers are available.
LOCAL_INCLUDE = os.path.join(HERE, "include")
REPO_INCLUDE = os.path.normpath(os.path.join(HERE, "..", "include"))
if os.path.isdir(REPO_INCLUDE):
    if os.path.isdir(LOCAL_INCLUDE):
        shutil.rmtree(LOCAL_INCLUDE)
    shutil.copytree(REPO_INCLUDE, LOCAL_INCLUDE)
if not os.path.isdir(LOCAL_INCLUDE):
    raise RuntimeError(
        "opal headers not found: expected ./include (sdist) or ../include (repo)"
    )

ext_modules = [
    Pybind11Extension(
        "opal._opal",
        ["bindings.cpp"],
        include_dirs=[LOCAL_INCLUDE],
        cxx_std=17,
    ),
]

setup(
    name="opal-pricing",
    version="0.2.0",
    description="Institutional option pricing library (C++ core, Python API)",
    long_description=(
        "Opal prices and risk-manages equity and interest rate options: "
        "Black-Scholes-Merton (incl. discrete cash dividends), Black-76, "
        "Bachelier, Heston, SABR, Hull-White; binomial/trinomial trees, "
        "Crank-Nicolson PDE, Monte Carlo and Longstaff-Schwartz American "
        "Monte Carlo engines; vanilla, digital, barrier, Asian and lookback "
        "payoffs; caps/floors and swaptions with multi-curve OIS discounting."
    ),
    long_description_content_type="text/plain",
    author="rtrimble13",
    url="https://github.com/rtrimble13/opal",
    license="MIT",
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Financial and Insurance Industry",
        "License :: OSI Approved :: MIT License",
        "Programming Language :: C++",
        "Programming Language :: Python :: 3",
        "Topic :: Office/Business :: Financial :: Investment",
    ],
    packages=["opal"],
    ext_modules=ext_modules,
    cmdclass={"build_ext": build_ext},
    python_requires=">=3.8",
    zip_safe=False,
)
