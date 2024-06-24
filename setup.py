from setuptools import setup, find_packages
from CanDI.__version__ import version

from pathlib import Path

this_directory = Path(__file__).parent
long_description = (this_directory / "README.rst").read_text()

setup(
    name='PyCanDI',
    description='A cancer data integration package',
    version=version,
    packages=find_packages(),
    long_description=long_description,
    long_description_content_type='text/x-rst',
    python_requires='>=3.11,<4.0',
    install_requires=[
        "pandas",
        "configparser",
        "requests",
        "tqdm",
    ],
    url = 'https://github.com/GilbertLabUCSF/CanDI',
    entry_points={
        'console_scripts': [
            'candi-install = CanDI.setup.install:main',
        ],
    },
    classifiers=[
        'Intended Audience :: Science/Research',
        'License :: OSI Approved :: MIT License',
        'Programming Language :: Python :: 3',
    ],
)
