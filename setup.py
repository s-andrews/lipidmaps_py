from setuptools import setup, find_packages

setup(
    name='lipidmaps-suite',  # Umbrella package name
    version='0.1.0',
    author='LIPID MAPS',
    description='LIPID MAPS Python API suite including reactions, BioPAN, and data processing tools',
    packages=find_packages(where='src'),  # This will find ALL packages under src/
    package_dir={'': 'src'},
    install_requires=[
        'pandas',
        'numpy',
        'requests',
    ],
    extras_require={
        'dev': [
            'pytest>=7.0',
            'pytest-cov',
            'black',
            'flake8',
            'mypy',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
