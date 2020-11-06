from pathlib import Path
from setuptools import setup, find_packages

SRC_DIR = Path(__file__).parent


def _get_version():
    env = {}
    with (SRC_DIR / "hamplicons" / "__init__.py").open() as handle:
        exec(handle.read(), env)

    return env["__version__"]


def _get_requirements():
    lines = []
    with (SRC_DIR / "requirements.txt").open() as handle:
        for line in handle:
            line = line.strip()
            if line:
                lines.append(line)

    return lines


setup(
    name="hamplicons",
    version=_get_version(),
    packages=find_packages(),
    description="Estimation of indel sizes in amplicons using Hamming distances",
    install_requires=_get_requirements(),
    entry_points={"console_scripts": ["hamplicons=hamplicons.main:entry_point"]},
    zip_safe=True,
)
