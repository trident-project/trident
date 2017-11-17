import os
from setuptools import setup
from setuptools import find_packages

# Utility function to read the README file.
# Used for the long_description.  It's nice, because now 1) we have a top level
# README file and 2) it's easier to type in the README file than to put a raw
# string in below ...
def read(fname):
    return open(os.path.join(os.path.dirname(__file__), fname)).read()

setup(
    name = "trident",
    version = "1.1.0",
    author = "Cameron Hummels, Devin Silvia, Britton Smith",
    author_email = "trident-project-users@googlegroups.com",
    description = ("Spectrum generator for astrophysical simulation data"),
    long_description = ("""Used for generating synthetic absorption-line 
spectra from astrophysical hydrodynamical data"""),
    license = "BSD",
    keywords = ["simulation", "spectra", "astronomy", "astrophysics"], 
    url = "https://trident-project.org",
    packages=find_packages(),
    include_package_data=True,
    classifiers=[
        "Development Status :: 5 - Production/Stable",
        "Intended Audience :: Science/Research",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
    install_requires=['yt>=3.4.0', 'h5py', 'numpy', 'matplotlib', 'astropy', 'requests', 'setuptools', 'future'],
)
