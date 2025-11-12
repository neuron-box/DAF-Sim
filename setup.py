"""Setup configuration for DAF Test Bench Framework."""

from setuptools import setup, find_packages
from pathlib import Path

# Read the README file
readme_file = Path(__file__).parent / "daf_test_bench" / "README.md"
long_description = readme_file.read_text(encoding="utf-8") if readme_file.exists() else ""

setup(
    name="daf-test-bench",
    version="1.0.0",
    author="DAF-Sim Development Team",
    author_email="",
    description="Unified testing and benchmarking framework for DAF simulation engines",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/neuron-box/DAF-Sim",
    packages=find_packages(include=["daf_test_bench", "daf_test_bench.*"]),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Chemistry",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Programming Language :: Python :: 3.10",
        "Programming Language :: Python :: 3.11",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
    install_requires=[
        "numpy>=1.22.0",
        "scipy>=1.10.0",
        "pandas>=1.3.0",
        "psutil>=5.9.0",
    ],
    extras_require={
        "engines": [
            "pillar3-physics-model",
            "floc-kinetics-pbm",
        ],
        "dev": [
            "pytest>=6.0",
            "pytest-cov>=2.12",
            "black>=21.0",
            "flake8>=3.9",
        ],
    },
    entry_points={
        "console_scripts": [
            "daf-benchmark=daf_benchmark:main",
        ],
    },
    include_package_data=True,
    package_data={
        "daf_test_bench": [
            "test_cases/*/*.json",
            "test_cases/*/*/*.csv",
        ],
    },
    zip_safe=False,
)
