/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_NBODYCOWELLSTATEDERIVATIVE_H
#define TUDAT_NBODYCOWELLSTATEDERIVATIVE_H

#include "Tudat/Astrodynamics/Propagators/nBodyStateDerivative.h"

namespace tudat
{

namespace propagators
{

template< typename StateScalarType = double, typename TimeType = double >
class NBodyCowellStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:
    NBodyCowellStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                                const boost::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                                const std::vector< std::string >& bodiesToIntegrate ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, cowell, bodiesToIntegrate ){ }

    ~NBodyCowellStateDerivative( ){ }

    //! Calculates the state derivative of the translational motion of the system.
    /*!
     * Calculates the state derivative (velocity+acceleration of each body) of the translational
     * motion of the system at the given time and position/velocity of bodies.
     * \param time Time (TDB seconds since J2000) at which the system is to be updated.
     * \param stateOfSystemToBeIntegrated List of 6 * bodiesToBeIntegratedNumerically_.size( ),
     * containing Caartesian position/velocity of the bodies being integrated. The order of the
     * values is defined by the order of bodies in bodiesToBeIntegratedNumerically_
     * \return Current state derivative (velocity+acceleration) of system of bodies integrated numerically.
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated )
    {
        return ( this->sumStateDerivativeContributions(
            stateOfSystemToBeIntegrated.template cast< double >( ),
            time ) ).template cast< StateScalarType >( );
    }

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
        convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic,
            Eigen::Dynamic >& cartesianSolution, const TimeType& time )
    {
        return cartesianSolution;
    }

    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
        convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic,
            Eigen::Dynamic >& internalSolution, const TimeType& time )
    {
        return internalSolution;
    }
};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_NBODYCOWELLSTATEDERIVATIVE_H
