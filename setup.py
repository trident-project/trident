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
    version = "0.4",
    author = "Cameron Hummels, Devin Silvia, Britton Smith",
    author_email = "chummels@gmail.com",
    description = ("Spectrum generator for astrophysical simulation data"),
    long_description = ("""Used for generating synthetic absorption-line 
spectra from astrophysical hydrodynamical data"""),
    license = "BSD",
    keywords = ["simulation", "spectra"], 
    url = "https://bitbucket.org/trident-project/trident",
    packages=find_packages(),
    classifiers=[
        "Development Status :: 4 - Beta",
        "Intended Audience :: Science/Research",
        "Topic :: Utilities",
        "License :: OSI Approved :: BSD License",
        "Natural Language :: English",
        "Programming Language :: Python :: 2.7",
    ],
    install_requires=['yt', 'h5py', 'numpy', 'matplotlib'],
)
