/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EIHPARTIALS_H
#define TUDAT_EIHPARTIALS_H

#include "tudat/astro/relativity/einsteinInfeldHoffmannEquations.h"

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

bool

//! Class to calculate the partials of the central gravitational acceleration w.r.t. parameters and states.
class EihEquationsPartials
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param gravitationalAcceleration Central gravitational acceleration w.r.t. which partials are to be taken.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     */
    EihEquationsPartials(
            const std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations );

    void update( const double currentTime )
    {

    }

protected:

    std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations_;


    std::vector< std::vector< Eigen::Matrix3d > > currentSingleAccelerationsWrtExertingPosition_;

    std::vector< std::vector< Eigen::Vector3d > > currentSingleAccelerationsWrtGravitationalParameter_;



    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentLocalPotentialsWrtPosition_;

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentLocalPotentialsWrtGravitationalParameter_;


    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentScalarTerm6WrtExertingPosition_;



    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentLocalPotentialsWrtGravitationalParameter_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_EIHPARTIALS_H
