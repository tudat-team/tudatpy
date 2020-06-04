
[![Build Status](https://img.shields.io/circleci/project/github/ggarrett13/tudatpy/master.svg?style=for-the-badge&logo=circleci)](https://circleci.com/gh/ggarrett13/tudatpy)
[![Conda](https://img.shields.io/conda/pn/ggarrett13/tudatpy?color=orange&logo=anaconda&style=for-the-badge)](https://anaconda.org/ggarrett13/tudatpy)


# TudatPy

TU Delft Astrodynamics Toolbox in Python, or **TudatPy**, is a library that primarily exposes the powerful set of C++ 
libraries, [Tudat](https://tudat.tudelft.nl/). TudatPy aims at accelerating the implementation of Tudat simulations,
providing and interface between Tudat and popular machine learning frameworks and establishing a platform to provide 
quality education in the field of astrodynamics. See the [documentation](https://ggarrett13.github.io/tudatpy/) for more.

Installation
===================

Installing `tudatpy` from the `ggarrett13` channel can be achieved by adding `ggarrett13` to your channels with:

```
conda config --add channels ggarrett13
```

Once the `ggarrett13` channel has been enabled, `tudatpy` can be installed with:

```
conda install tudatpy
```

It is possible to list all of the versions of `tudatpy` available on your platform with:

```
conda search tudatpy --channel ggarrett13
```

## Creating an issue on GitHub

post an issue on the [TudatPy](https://github.com/ggarrett13/tudatpy)
GitHub page using the issue template ``<tudatpy>/docs/issue_template.md``.