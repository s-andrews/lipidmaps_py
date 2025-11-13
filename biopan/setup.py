from setuptools import setup, find_packages

setup(
    name='biopan-framework',
    version='0.1.0',
    author='Chetin Baloglu',
    author_email='baloglu@gmail.com',
    description='A backend tool for processing lipid data in BioPAN.',
    packages=find_packages(where='src'),
    package_dir={'': 'src'},
    install_requires=[
        'pandas',
        'numpy',
        'requests',
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)