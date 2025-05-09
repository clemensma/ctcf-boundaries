# setup.py

from setuptools import setup, find_packages

setup(
    name='ctcf-boundaries',
    version='0.1',
    description='Command-line tool to detect CTCF boundaries from .cool and BED files',
    author='Clemens Mauksch',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'cooler',
        'cooltools',
        'bioframe'
    ],
    entry_points={
        'console_scripts': [
            'ctcf-boundaries=ctcf_boundaries.cli:main',
        ],
    },
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
    ],
    python_requires='>=3.8',
)