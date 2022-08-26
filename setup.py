from setuptools import setup
from setuptools import find_packages

def get_version(filename):
    """
    Get version from a file.

    Inspired by https://github.mabuchilab/QNET/.
    """
    with open(filename) as f:
        for line in f.readlines():
            if line.startswith("__version__"):
                return line.split("=")[1].strip()[1:-1]
    raise RuntimeError(
        "Could not get version from %s." % filename)


VERSION = get_version("trident/__init__.py")

with open('README.md') as f:
    long_description = f.read()

dev_requirements = [
    'coveralls', 'flake8', 'pytest>=3.6', 'pytest-cov', 'twine', 'wheel',
    'sphinx', 'sphinx_rtd_theme', 'gitpython']

setup(
    name = "trident",
    version = VERSION,
    description = ("Spectrum generator for astrophysical simulation data"),
    long_description=long_description,
    long_description_content_type='text/markdown',
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
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.7",
        "Programming Language :: Python :: 3.8",
        "Programming Language :: Python :: 3.9",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
      extras_require={
          'dev': dev_requirements,
      },
    install_requires=[
        'astropy',
        'h5py',
        'matplotlib',
        'numpy',
        'requests',
        'scipy',
        'yt>=4.0.1',
        'yt_astro_analysis>=1.1.0'],
      python_requires='>=3.7'
)
