/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_DEFAULT_H
#define TUDATPY_DEFAULT_H

#include "../../expose_astro/prototype/frames.h"
#include "simpleBody.h"
#include <map>
namespace tudat {

namespace bodies {



namespace solar_system {

std::shared_ptr<SimpleBody> Sun = std::make_shared<SimpleBody>(
    1.32712440018E20,// gravitationalParameter
    "Sun",           // name
    "\u2609",        // symbol
    0.,              // radius
    0.,              // polarRadius
    0.,              // meanRadius
    0.,              // rotationalPeriod
    0.,              // J2
    0.,              // J3
    0.,              // mass
    0.,              // atmosphereDistance
    0.               // flybyLimit
);

std::shared_ptr<SimpleBody> Mercury = std::make_shared<SimpleBody>(
    Sun,       // parent
    2.20329E13,// gravitationalParameter
    "Mercury", // name
    "Mercury", // symbol
    0.,        // radius
    0.,        // polarRadius
    0.,        // meanRadius
    0.,        // rotationalPeriod
    0.,        // J2
    0.,        // J3
    0.,        // mass
    0.,        // atmosphereDistance
    0.         // flybyLimit
);

std::shared_ptr<SimpleBody> Venus = std::make_shared<SimpleBody>(
    Sun,        // parent
    3.248599E14,// gravitationalParameter
    "Venus",    // name
    "Venus",    // symbol
    0.,         // radius
    0.,         // polarRadius
    0.,         // meanRadius
    0.,         // rotationalPeriod
    0.,         // J2
    0.,         // J3
    0.,         // mass
    0.,         // atmosphereDistance
    0.          // flybyLimit
);

std::shared_ptr<SimpleBody> Earth = std::make_shared<SimpleBody>(
    Sun,            // parent
    3.9860044188E14,// gravitationalParameter
    "Earth",        // name
    "Earth",        // symbol
    0.,             // radius
    0.,             // polarRadius
    0.,             // meanRadius
    0.,             // rotationalPeriod
    0.,             // J2
    0.,             // J3
    0.,             // mass
    0.,             // atmosphereDistance
    0.              // flybyLimit

);

std::shared_ptr<SimpleBody> Mars = std::make_shared<SimpleBody>(
    Sun,         // parent
    4.2828372E13,// gravitationalParameter
    "Mars",      // name
    "Mars",      // symbol
    0.,          // radius
    0.,          // polarRadius
    0.,          // meanRadius
    0.,          // rotationalPeriod
    0.,          // J2
    0.,          // J3
    0.,          // mass
    0.,          // atmosphereDistance
    0.           // flybyLimit

);

SimpleSystemOfBodies SolarSystem = SimpleSystemOfBodies(
    Sun,
    std::make_shared<ReferenceFrame>(Sun->getName(), "ECLIPJ2000"),
    {Mercury, Venus, Earth, Mars});

}// namespace solar_system



}// namespace bodies

}// namespace tudat
#endif//TUDATPY_DEFAULT_H
