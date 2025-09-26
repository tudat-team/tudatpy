/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GAS_PROPERTIES_H
#define TUDAT_GAS_PROPERTIES_H

namespace tudat
{

namespace aerodynamics
{

//! Molecular collision diameter of Argon in m
constexpr static double DIAMETER_ARGON = 340E-12;

//! Molecular collision diameter of Atomic Hydrogen in m
constexpr static double DIAMETER_ATOMIC_HYDROGEN = 260E-12;

//! Molecular collision diameter of Helium in m
constexpr static double DIAMETER_HELIUM = 256E-12;

//! Molecular collision diameter of Nitrogen in m
constexpr static double DIAMETER_NITROGEN = 370E-12;

//! Molecular collision diameter of Oxygen in m
constexpr static double DIAMETER_OXYGEN = 358E-12;

//! Molecular collision diameter of Atomic Nitrogen in m
constexpr static double DIAMETER_ATOMIC_NITROGEN = 290E-12;

//! Molecular collision diameter of Atomic Oxygen in m
constexpr static double DIAMETER_ATOMIC_OXYGEN = 280E-12;

//! Molar mass of Argon in kg/mole
constexpr static double MOLAR_MASS_ARGON = 39.948E-3;

//! Molar mass of Atomic Hydrogen in kg/mole
constexpr static double MOLAR_MASS_ATOMIC_HYDROGEN = 1.008E-3;

//! Molar mass of Helium in kg/mole
constexpr static double MOLAR_MASS_HELIUM = 4.002602E-3;

//! Molar mass of Nitrogen in kg/mole
constexpr static double MOLAR_MASS_NITROGEN = 2.0 * 14.007E-3;

//! Molar mass of Oxygen in kg/mole
constexpr static double MOLAR_MASS_OXYGEN = 2.0 * 15.999E-3;

//! Molar mass of Atomic Nitrogen in kg/mole
constexpr static double MOLAR_MASS_ATOMIC_NITROGEN = 14.007E-3;

//! Molar mass of Atomic Oxygen in kg/mole
constexpr static double MOLAR_MASS_ATOMIC_OXYGEN = 15.999E-3;

}  // namespace aerodynamics
}  // namespace tudat

#endif  // TUDAT_GAS_PROPERTIES_H_
