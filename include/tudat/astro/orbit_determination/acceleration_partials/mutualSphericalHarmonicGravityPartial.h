/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MUTUALSPHERICALHARMONICGRAVITYPARTIAL_H
#define TUDAT_MUTUALSPHERICALHARMONICGRAVITYPARTIAL_H

#include "tudat/astro/gravitation/mutualSphericalHarmonicGravityModel.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"


namespace tudat
{

namespace acceleration_partials
{
//! Class for calculating partial derivatives of a mutual spherical harmonic gravitational acceleration.
/*!
 *  Class for calculating partial derivatives of a mutual spherical harmonic gravitational acceleration, as calculated by the
 *  MutualSphericalHarmonicsGravitationalAccelerationModel class. It incorporate the extended body - point mass and point mass -
 *  extended body interactions.
 */
class MutualSphericalHarmonicsGravityPartial: public AccelerationPartial
{
public:

    //! Constructor
    /*!
     * Constructor from partial objects of the two spherical harmonic gravity accelerations that nuild up the mutual acceleration.
     * Note that the C_{00} coefficient for one of the bodies must be 0, and that their gravitational parameters must be equal.
     * \param accelerationPartialOfShExpansionOfBodyExertingAcceleration Partial for spherical harmonic acceleration due to body
     * exerting acceleration
     * \param accelerationPartialOfShExpansionOfBodyUndergoingAcceleration Partial for spherical harmonic acceleration due to body
     * undergoing acceleration
     * \param acceleratedBody Name of body undergoing acceleration
     * \param acceleratingBody Name of body exerting acceleration
     * \param accelerationUsesMutualAttraction Variable denoting whether point mass attraction from body undergoing acceleration
     * on  body exerting acceleration is included (i.e. whether gravitational parameter is the property
     * of the body exerting the acceleration, if variable is false, or the sum of the gravitational parameters,
     * if the variable is true.
     */
    MutualSphericalHarmonicsGravityPartial(
            const std::shared_ptr< SphericalHarmonicsGravityPartial >
            accelerationPartialOfShExpansionOfBodyExertingAcceleration,
            const std::shared_ptr< SphericalHarmonicsGravityPartial >
            accelerationPartialOfShExpansionOfBodyUndergoingAcceleration,
            const std::string& acceleratedBody, const std::string& acceleratingBody, const bool accelerationUsesMutualAttraction ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::mutual_spherical_harmonic_gravity ),
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_(
            accelerationPartialOfShExpansionOfBodyExertingAcceleration ),
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_(
            accelerationPartialOfShExpansionOfBodyUndergoingAcceleration ),
        accelerationUsesMutualAttraction_( accelerationUsesMutualAttraction ){ }

    //! Destructor
    ~MutualSphericalHarmonicsGravityPartial( ) { }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratedBody(
                Eigen::Block< Eigen::MatrixXd > partialMatrix,
                const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->wrtPositionOfAcceleratedBody(
                    partialMatrix, addContribution, startRow, startColumn );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->wrtPositionOfAcceleratingBody(
                    partialMatrix, !addContribution, startRow, startColumn );
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->wrtPositionOfAcceleratingBody(
                    partialMatrix, addContribution, startRow, startColumn );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->wrtPositionOfAcceleratedBody(
                    partialMatrix, !addContribution, startRow, startColumn );
    }

    //! Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
     *  and adding it to the existing partial block. Function calls constituent spherical harmonic model functions
     *  \param partialMatrix Block of partial derivatives of where current partial is to be added.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which partial is to be computed.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     */
    void wrtNonTranslationalStateOfAdditionalBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType,
            const bool addContribution )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->
                        wrtNonTranslationalStateOfAdditionalBody(
                            partialMatrix, stateReferencePoint, integratedStateType, true );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->
                        wrtNonTranslationalStateOfAdditionalBody(
                            partialMatrix, stateReferencePoint, integratedStateType, false );
    }

    //! Function for determining if the acceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  No dependency is implemented, but a warning is provided if partial w.r.t. mass of body exerting acceleration
     *  (and undergoing acceleration if mutual attraction is used) is requested.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
                const std::pair< std::string, std::string >& stateReferencePoint,
                const propagators::IntegratedStateType integratedStateType )
    {
        if( accelerationPartialOfShExpansionOfBodyExertingAcceleration_->
                isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) ||
                accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->
                isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
        {
            return true;
        }
        else
        {
            return false;
        }
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function<void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );


    //! Function to calculate the partial wrt the gravitational parameter.
    /*!
     *  Function to calculate the partial wrt the gravitational parameter of the central body. Note that in the case of
     *  mutual attraction (see SphericalHarmonicsGravitationalAccelerationModel), the partial wrt the gravitational
     *  parameter of the body exerting acceleration is equal to the partial wrt to the body undergoing the acceleration,
     *  which is zero otherwise.
     *  \return Partial wrt the gravitational parameter
     */
    void wrtGravitationalParameter( Eigen::MatrixXd& partialMatrix )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->wrtGravitationalParameterOfCentralBody(
                    partialMatrix, 0 );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->wrtGravitationalParameterOfCentralBody(
                    partialMatrix, -1 );
    }

    //! Function for updating the partial object to current state and time.
    /*!
     *  Function for updating the partial object to current state and time. Calculates the variables that are
     *  used for the calculation of multple partials, to prevent multiple calculations of same function.
     *  \param currentTime Time to which object is to be updated
     */
    void update( const double currentTime )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->update( currentTime );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->update( currentTime );

        currentTime_ = currentTime;

    }

    //! Function to set a dependency of this partial object w.r.t. a given double parameter.
    /*!
     * Function to set a dependency of this partial object w.r.t. a given double parameter. If a dependency exists, the given
     * partial is recomputed on every call of updateParameterPartials.
     * \param parameter Partial w.r.t. which dependency is to be checked and set.
     * \return Size (number of columns) of parameter partial. Zero if no dependency, 1 otherwise.
     */
    int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function to set a dependency of this partial object w.r.t. a given vector parameter.
    /*!
     * Function to set a dependency of this partial object w.r.t. a given vector parameter. If a dependency exists, the given
     * partial is recomputed on every call of updateParameterPartials.
     * \param parameter Partial w.r.t. which dependency is to be checked and set.
     * \return Size (number of columns) of parameter partial. Zero if no dependency, size of parameter otherwise.
     */
    int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

protected:

    //! Function to reset the constituent SphericalHarmonicsGravityPartial objects to the current time.
    void resetCurrentTimeOfMemberObjects( )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->resetCurrentTime( );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->resetCurrentTime( );
    }

    //! Function to update the parameter partials of the constituent SphericalHarmonicsGravityPartial objects.
    void updateParameterPartialsOfMemberObjects( )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->updateParameterPartials( );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->updateParameterPartials( );
    }

    //! Partial for spherical harmonic acceleration due to body exerting acceleration
    std::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialOfShExpansionOfBodyExertingAcceleration_;

    //!  Partial for spherical harmonic acceleration due to body undergoing acceleration
    std::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_;

    //! Boolean denoting whether the mutual point mass attraction between the bodies is taken into account
    /*!
     *  Boolean denoting whether the mutual point mass  attraction between the bodies is taken into account, as is the case when
     *  integrting in a (non-rotating) frame centered on the body exerting the acceleration, in which case the
     *  gravitational acceleration of the body undergoing the acceleration on that exerting the acceleration must be taken
     *  into account as an inertial effect.
     */
    bool accelerationUsesMutualAttraction_;

};

}

}

#endif // TUDAT_MUTUALSPHERICALHARMONICGRAVITYPARTIAL_H
