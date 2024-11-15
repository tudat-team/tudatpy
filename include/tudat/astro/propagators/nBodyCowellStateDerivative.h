/*    Copyright (c) 2010-2019, Delft University of Technology
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

#include "tudat/astro/propagators/nBodyStateDerivative.h"

namespace tudat
{

namespace propagators
{

//! Class for computing the state derivative of translational motion of N bodies, using a Cowell propagator.
template< typename StateScalarType = double, typename TimeType = double >
class NBodyCowellStateDerivative: public NBodyStateDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     *  \param accelerationModelsPerBody A map containing the list of accelerations acting on each
     *  body, identifying the body being acted on and the body acted on by an acceleration. The map
     *  has as key a string denoting the name of the body the list of accelerations, provided as the
     *  value corresponding to a key, is acting on.  This map-value is again a map with string as
     *  key, denoting the body exerting the acceleration, and as value a pointer to an acceleration
     *  model.
     *  \param centralBodyData Object responsible for providing the current integration origins from
     *  the global origins.
     *  \param bodiesToIntegrate List of names of bodies that are to be integrated numerically.
     */
    NBodyCowellStateDerivative( const basic_astrodynamics::AccelerationMap& accelerationModelsPerBody,
                                const std::shared_ptr< CentralBodyData< StateScalarType, TimeType > > centralBodyData,
                                const std::vector< std::string >& bodiesToIntegrate ):
        NBodyStateDerivative< StateScalarType, TimeType >(
            accelerationModelsPerBody, centralBodyData, cowell, bodiesToIntegrate ){ }

    //! Destructor
    ~NBodyCowellStateDerivative( ){ }

    //! Calculates the state derivative of the translational motion of the system.
    /*!
     * Calculates the state derivative (velocity+acceleration of each body) of the translational motion of the system
     * at the given time and position/velocity of bodies.
     *  \param time Time (TDB seconds since J2000) at which the system is to be updated.
     *  \param stateOfSystemToBeIntegrated List of 6 * bodiesToBeIntegratedNumerically_.size( ), containing Caartesian
     *  position/velocity of the bodies being integrated. The order of the values is defined by the order of bodies in
     *  bodiesToBeIntegratedNumerically_
     *  \param stateDerivative Current state derivative (velocity+acceleration) of system of bodies integrated numerically
     *  (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time, const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative.setZero( );
        this->sumStateDerivativeContributions( stateOfSystemToBeIntegrated, stateDerivative, true );
    }

    //! Function to convert the state in the conventional form to the propagator-specific form.
    /*!
     * Function to convert the state in the conventional form to the propagator-specific form. For the Cowell propagator,
     * the two are equivalent, and this function returns the input state.
     * \param cartesianSolution State in 'conventional form'
     * \param time Current time at which the state is valid (not used in this class).
     * \return State (outputSolution), converted to the 'propagator-specific form' (which is equal to outputSolution).
     */
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >
        convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& cartesianSolution, const TimeType& time )
    {
        return cartesianSolution;
    }

    //! Function to convert the propagator-specific form of the state to the conventional form.
    /*!
     * Function to convert the propagator-specific form of the state to the conventional form. For the Cowell propagator,
     * the two are equivalent, and this function returns the input state.
     * In contrast to the convertCurrentStateToGlobalRepresentation function, this
     * function does not provide the state in the inertial frame, but instead provides it in the
     * frame in which it is propagated.
     * \param internalSolution State in propagator-specific form (i.e. form that is used in
     * numerical integration, equal to conventional form for this class).
     * \param time Current time at which the state is valid (not used in this class).
     * \param currentCartesianLocalSoluton State (internalSolution), converted to the 'conventional form',
     *  which is equal to outputSolution for this class (returned by reference).
     */
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        currentCartesianLocalSoluton = internalSolution;
    }

};

//extern template class NBodyCowellStateDerivative< double, double >;

} // namespace propagators

} // namespace tudat

#endif // TUDAT_NBODYCOWELLSTATEDERIVATIVE_H
