import os
import re
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


def read_version():
    """Derive the package version from the single source of truth,
    include/opal/version.hpp, so the C++ macros, CMake and the wheel never
    disagree (issue #19). The header is bundled into the package above, so it
    is present both in the repo and inside an sdist."""
    for inc in (LOCAL_INCLUDE, REPO_INCLUDE):
        header = os.path.join(inc, "opal", "version.hpp")
        if os.path.isfile(header):
            with open(header) as fh:
                m = re.search(r'OPAL_VERSION_STRING\s+"([^"]+)"', fh.read())
            if m:
                return m.group(1)
    raise RuntimeError("could not parse OPAL_VERSION_STRING from version.hpp")


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
    version=read_version(),
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
