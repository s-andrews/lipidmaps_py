from setuptools import setup, find_packages

setup(
    name="lipidmaps_py",  # Umbrella package name
    version="0.1.0",
    author="LIPID MAPS",
    description="LIPID MAPS Python API suite for data input, normalization, data processing and LIPID MAPS reactions",
    packages=find_packages(where="src"),  # This will find ALL packages under src/
    package_dir={"": "src"},
    install_requires=[
        "pydantic",
        "pandas",
        "numpy",
        "requests",
        "networkx",
        "matplotlib",
    ],
    extras_require={
        "dev": [
            "pytest>=7.0",
            "pytest-cov",
            "black",
            "flake8",
            "mypy",
        ],
    },
    classifiers=[
        "Programming Language :: Python :: 3",
        "License :: OSI Approved :: MIT License",
        "Operating System :: OS Independent",
    ],
    python_requires=">=3.8",
)
