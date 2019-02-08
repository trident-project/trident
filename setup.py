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
    'sphinx', 'scipy', 'sphinx_rtd_theme', 'gitpython']

setup(
    name = "trident",
    version = VERSION,
    author = "Cameron Hummels, Devin Silvia, Britton Smith",
    author_email = "trident-project-users@googlegroups.com",
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
        "Programming Language :: Python :: 2",
        "Programming Language :: Python :: 2.7",
        "Programming Language :: Python :: 3",
        "Programming Language :: Python :: 3.5",
        "Programming Language :: Python :: 3.6",
        "Topic :: Scientific/Engineering :: Astronomy",
        "Topic :: Scientific/Engineering :: Physics",
        "Topic :: Scientific/Engineering :: Visualization",
    ],
      extras_require={
          'dev': dev_requirements,
      },
    install_requires=[
        'astropy',
        'future',
        'h5py',
        'matplotlib',
        'numpy!=1.14.0',
        'requests',
        'yt>=3.4.0',
        'yt_astro_analysis'],
      python_requires='>=2.7,!=3.0.*,!=3.1.*,!=3.2.*,!=3.3.*,!=3.4.*'
)
