[![Build Status](https://img.shields.io/circleci/project/github/tudat-team/tudatpy/master.svg?style=for-the-badge&logo=circleci)](https://circleci.com/gh/tudat-team/tudatpy)
[![Conda](https://img.shields.io/conda/pn/tudat-team/tudatpy?color=orange&logo=anaconda&style=for-the-badge)](https://anaconda.org/tudat-team/tudatpy) [![Join the chat at https://gitter.im/tudat-space/tudatpy](https://badges.gitter.im/tudat-space/tudatpy.svg)](https://gitter.im/tudat-space/tudatpy?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)

# tudatpy

TU Delft Astrodynamics Toolbox in Python, or **tudatpy**, is a library that primarily exposes the powerful set of C++ 
libraries, [Tudat](https://tudat.tudelft.nl/). TudatPy aims at accelerating the implementation of Tudat simulations,
providing an interface between Tudat and popular machine learning frameworks and establishing a platform to provide 
quality education in the field of astrodynamics. See the [documentation](https://tudat-space.readthedocs.io) for more.

Installation
===================

Installing `tudatpy` from the `tudat-team` channel can be achieved by adding `tudat-team` to your channels with:

```
conda config --add channels tudat-team
```

Once the `tudat-team` channel has been enabled, `tudatpy` can be installed with:

```
conda install tudatpy
```

It is possible to list all of the versions of `tudatpy` available on your platform with:

```
conda search tudatpy --channel tudat-team
```
