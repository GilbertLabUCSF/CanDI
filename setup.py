from setuptools import setup, find_packages
from CanDI.__version__ import version

from pathlib import Path

this_directory = Path(__file__).parent

setup(
    name='CanDI',
    description='A cancer data integration package',
    version=version,
    packages=find_packages(),
    python_requires='>=3.9',
    install_requires=[
        "pandas",
        "configparser",
        "requests",
        "tqdm",
    ],
    url = 'https://github.com/GilbertLabUCSF/CanDI',
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
)