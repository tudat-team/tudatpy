/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_THIRDBODYGRAVITYPARTIAL_H
#define TUDAT_THIRDBODYGRAVITYPARTIAL_H

#include <memory>

#include "tudat/basics/tudatTypeTraits.h"

#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/mutualSphericalHarmonicGravityPartial.h"

#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Function to get the third body acceleration type from the direct acceleration partial object.
/*!
 *  Function to get the third body acceleration type from the direct acceleration partial object.
 *  \param directGravityPartial Partial derivative of direct acceleration.
 *  \return Type of acceleration for third body acceleration for which direct acceleration partial is given by
 *  directGravityPartial (i.e. third_body_point_mass_gravity if input is of type CentralGravitationPartial).
 */
template< typename DirectGravityPartial >
basic_astrodynamics::AvailableAcceleration getAccelerationTypeOfThirdBodyGravity(
        const std::shared_ptr< DirectGravityPartial > directGravityPartial )
{
    using namespace basic_astrodynamics;
    AvailableAcceleration accelerationType;

    // Check type of direct partial derivative.
    if( std::dynamic_pointer_cast< CentralGravitationPartial >( directGravityPartial ) != nullptr )
    {
        accelerationType = third_body_point_mass_gravity;
    }
    else if( std::dynamic_pointer_cast< SphericalHarmonicsGravityPartial >( directGravityPartial ) != nullptr )
    {
        accelerationType = third_body_spherical_harmonic_gravity;
    }
    else if( std::dynamic_pointer_cast< MutualSphericalHarmonicsGravityPartial >( directGravityPartial ) != nullptr )
    {
        accelerationType = third_body_mutual_spherical_harmonic_gravity;
    }

    else
    {
        throw std::runtime_error( "Error when getting third body partial type, type not identified" );
    }
    return accelerationType;
}

//! Class to calculate the partials of a third-body gravitational acceleration w.r.t. parameters and states.
/*!
 *  Class to calculate the partials of a third-body gravitational acceleration w.r.t. parameters and states. This class may
 *  be used for any direct gravitational acceleration (central, spherical harmonic, mutual spherical harmonic, etc.),
 *  providin a generic third-body partial interface. The template parameter is the derived class of AccelerationPartial
 *  for the associated direct acceleration partial.
 */
template< typename DirectGravityPartial,
                    typename std::enable_if< is_direct_gravity_partial< DirectGravityPartial >::value, int >::type = 0 >
class ThirdBodyGravityPartial: public AccelerationPartial
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param partialOfDirectGravityOnBodyUndergoingAcceleration Partial derivative of direct acceleration from
     * acceleratingBody on acceleratedBody.
     * \param partialOfDirectGravityOnCentralBody Partial derivative of direct acceleration from
     * centralBodyName on acceleratedBody.
     * \param acceleratedBody Name of body undergoing acceleration
     * \param acceleratingBody Name of body exerting acceleration (third-body)
     * \param centralBodyName Name of central body w.r.t. which the acceleration is computed.
     */
    ThirdBodyGravityPartial(
            const std::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnBodyUndergoingAcceleration,
            const std::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnCentralBody,
            const std::string& acceleratedBody, const std::string& acceleratingBody,
            const std::string& centralBodyName ):
        AccelerationPartial( acceleratedBody, acceleratingBody, getAccelerationTypeOfThirdBodyGravity(
                                 partialOfDirectGravityOnBodyUndergoingAcceleration ) ),
        partialOfDirectGravityOnBodyUndergoingAcceleration_( partialOfDirectGravityOnBodyUndergoingAcceleration ),
        partialOfDirectGravityOnCentralBody_( partialOfDirectGravityOnCentralBody ), centralBodyName_( centralBodyName ){ }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body undergoing acceleration.
     */
    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAcceleratedBody(
                    partialMatrix, addContribution, startRow, startColumn );

        // Check if acceleration on central body is dependent on acceleratedBody_
        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonnullptr( acceleratedBody_ ) == 1 )
        {
            partialOfDirectGravityOnCentralBody_->wrtPositionOfAdditionalBody(
                        acceleratedBody_, partialMatrix, addContribution, startRow, startColumn   );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. velocity of body undergoing acceleration.
     */
    void wrtVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                       const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAcceleratedBody(
                    partialMatrix, addContribution, startRow, startColumn  );

        // Check if acceleration on central body is dependent on acceleratedBody_
        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonnullptr( acceleratedBody_ ) == 1 )
        {
            partialOfDirectGravityOnCentralBody_->wrtVelocityOfAdditionalBody(
                        acceleratedBody_, partialMatrix, addContribution, startRow, startColumn  );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. position of body exerting acceleration.
     */
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )

    {
        // Add partials for both direct acceleration and acceleration on central body.
        partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAcceleratingBody(
                    partialMatrix, addContribution, startRow, startColumn  );
        partialOfDirectGravityOnCentralBody_->wrtPositionOfAcceleratingBody(
                    partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn  );

    }


    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \return Partial derivative of acceleration w.r.t. velocity of body exerting acceleration.
     */
    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        // Add partials for both direct acceleration and acceleration on central body.
        partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAcceleratingBody(
                    partialMatrix, addContribution, startRow, startColumn );
        partialOfDirectGravityOnCentralBody_->wrtVelocityOfAcceleratingBody(
                    partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn );
    }

    //! Function for calculating the partial of the acceleration w.r.t. the position of an additiona body.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of an additional body (in addition
     *  to the body undergoing and exerting the acceleration) and adding it to the existing partial block.
     *  This function check if the requested additional body equals the central body name. Also, it checks whether either
     *  the direct or indirect acceleration depend on this additional body. The partial computation is then updated
     *  accordingly
     *  \param bodyName Name of additional body.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of third body where
     *  current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtPositionOfAdditionalBody( const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                      const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( bodyName == centralBodyName_ )
        {
            partialOfDirectGravityOnCentralBody_->wrtPositionOfAcceleratedBody(
                        partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn  );
        }

        if( partialOfDirectGravityOnBodyUndergoingAcceleration_->isAccelerationPartialWrtAdditionalBodyNonnullptr( bodyName ) )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtPositionOfAdditionalBody(
                        bodyName, partialMatrix, addContribution, startRow, startColumn  );
        }

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonnullptr( bodyName ) )
        {
            partialOfDirectGravityOnCentralBody_->wrtPositionOfAdditionalBody(
                        bodyName, partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn );
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of an additiona body.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of an additional body (in addition
     *  to the body undergoing and exerting the acceleration) and adding it to the existing partial block.
     *  This function check if the requested additional body equals the central body name. Also, it checks whether either
     *  the direct or indirect acceleration depend on this additional body. The partial computation is then updated
     *  accordingly
     *  \param bodyName Name of additional body.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of third body where
     *  current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAdditionalBody( const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                      const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        if( bodyName == centralBodyName_ )
        {
            partialOfDirectGravityOnCentralBody_->wrtVelocityOfAcceleratedBody(
                        partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn  );
        }

        if( partialOfDirectGravityOnBodyUndergoingAcceleration_->isAccelerationPartialWrtAdditionalBodyNonnullptr( bodyName ) )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->wrtVelocityOfAdditionalBody(
                        bodyName, partialMatrix, addContribution, startRow, startColumn );
        }

        if( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonnullptr( bodyName ) )
        {
            partialOfDirectGravityOnCentralBody_->wrtVelocityOfAdditionalBody(
                        bodyName, partialMatrix, ( ( addContribution == true ) ? ( false ) : ( true ) ), startRow, startColumn );
        }
    }

    //! Function for checking whether the partial of the acceleration w.r.t. the state of an additional body is non-zero.
    /*!
     *  Function for checking whether the partial of the acceleration w.r.t. the state of an additional body is non-zero
     *  \param bodyName Name of additional body.
     *  \return True if dependency exists
     */
    bool isAccelerationPartialWrtAdditionalBodyNonnullptr( const std::string& bodyName )
    {
        bool isAccelerationDependentOnBody = 0;
        if( bodyName == centralBodyName_ )
        {
            isAccelerationDependentOnBody = 1;
        }

        if( ( partialOfDirectGravityOnCentralBody_->isAccelerationPartialWrtAdditionalBodyNonnullptr( bodyName ) == 1 ) ||
                ( partialOfDirectGravityOnBodyUndergoingAcceleration_->isAccelerationPartialWrtAdditionalBodyNonnullptr( bodyName ) == 1 ) )
        {
            isAccelerationDependentOnBody = 1;
        }

        return isAccelerationDependentOnBody;
    }

    //! Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. a non-translational integrated state
     *  and adding it to the existing partial block.
     *  \param partialMatrix Block of partial derivatives of where current partial is to be added.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which partial is to be computed.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     */
    void wrtNonTranslationalStateOfAdditionalBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType,
            const bool addContribution = true )
    {
        partialOfDirectGravityOnCentralBody_->
                        wrtNonTranslationalStateOfAdditionalBody(
                            partialMatrix, stateReferencePoint, integratedStateType, true );
        partialOfDirectGravityOnBodyUndergoingAcceleration_->
                        wrtNonTranslationalStateOfAdditionalBody(
                            partialMatrix, stateReferencePoint, integratedStateType, false );
    }

    //! Function for determining if the acceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        if( partialOfDirectGravityOnCentralBody_->
                isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
                    stateReferencePoint, integratedStateType ) ||
                partialOfDirectGravityOnBodyUndergoingAcceleration_->
                isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
                    stateReferencePoint, integratedStateType ) )
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
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );

        return orbit_determination::createMergedParameterPartialFunction(
                    partialFunctionFromDirectGravity, partialFunctionFromCentralGravity );
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );

        return orbit_determination::createMergedParameterPartialFunction(
                    partialFunctionFromDirectGravity, partialFunctionFromCentralGravity );
    }

    //! Function to set a dependency of this partial object w.r.t. a given double parameter.
    /*!
     * Function to set a dependency of this partial object w.r.t. a given double parameter. If a dependency exists, the given
     * partial is recomputed on every call of updateParameterPartials.
     * \param parameter Partial w.r.t. which dependency is to be checked and set.
     * \return Size (number of columns) of parameter partial. Zero if no dependency, 1 otherwise.
     */
    int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        // Check parameter dependency of direct acceleration.
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromDirectGravity.second > 0 )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        // Check parameter dependency of indirect acceleration.
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );
        if( partialFunctionFromCentralGravity.second > 0 )
        {
            partialOfDirectGravityOnCentralBody_->setParameterPartialUpdateFunction( parameter );
        }

        // If either dependency exists, create combined partial.
        if( partialFunctionFromCentralGravity.second > 0 || partialFunctionFromDirectGravity.second > 0 )
        {
            parameterDoublePartialFunctions_[ parameter ] =
                    getCombinedCurrentDoubleParameterFunction(
                        partialOfDirectGravityOnBodyUndergoingAcceleration_,
                        partialOfDirectGravityOnCentralBody_,
                        parameter, partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second, 1 );
            isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
            currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, 1 );
        }
        return std::max( partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second );
    }

    //! Function to set a dependency of this partial object w.r.t. a given vector parameter.
    /*!
     * Function to set a dependency of this partial object w.r.t. a given vector parameter. If a dependency exists, the given
     * partial is recomputed on every call of updateParameterPartials.
     * \param parameter Partial w.r.t. which dependency is to be checked and set.
     * \return Size (number of columns) of parameter partial. Zero if no dependency, size of parameter otherwise.
     */
    int setParameterPartialUpdateFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        // Check parameter dependency of direct acceleration.
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromDirectGravity =
                partialOfDirectGravityOnBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromDirectGravity.second > 0 )
        {
            partialOfDirectGravityOnBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        // Check parameter dependency of indirect acceleration.
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromCentralGravity =
                partialOfDirectGravityOnCentralBody_->getParameterPartialFunction( parameter );
        if( partialFunctionFromCentralGravity.second > 0 )
        {
            partialOfDirectGravityOnCentralBody_->setParameterPartialUpdateFunction( parameter );
        }

        // If either dependency exists, create combined partial.
        if( partialFunctionFromCentralGravity.second > 0 || partialFunctionFromDirectGravity.second > 0 )
        {
            parameterVectorPartialFunctions_[ parameter ] =
                    getCombinedCurrentVectorParameterFunction(
                        partialOfDirectGravityOnBodyUndergoingAcceleration_,
                        partialOfDirectGravityOnCentralBody_,
                        parameter, partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second, 1 );
            isCurrentVectorParameterPartialSet_[ parameter ] = 0;
            currentVectorParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, parameter->getParameterSize( ) );

        }
        return std::max( partialFunctionFromDirectGravity.second, partialFunctionFromCentralGravity.second );
    }

    //! Function for updating partials  w.r.t. the bodies' positions
    /*!
     *  Function for updating common blocks of partial to current state. For the third body gravitational acceleration,
     *  the update functions of both constituent accelerations are updated.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->update( currentTime );
        partialOfDirectGravityOnCentralBody_->update( currentTime );

        currentTime_ = currentTime;
    }

    //! Function to get partial derivative object of direct acceleration from centralBodyName on acceleratedBody.
    /*!
     * Function to get the partial derivative object of direct acceleration from centralBodyName on acceleratedBody.
     * \return Partial derivative object of direct acceleration from centralBodyName on acceleratedBody.
     */
    std::shared_ptr< DirectGravityPartial > getPartialOfDirectGravityOnBodyUndergoingAcceleration( )
    {
        return partialOfDirectGravityOnBodyUndergoingAcceleration_;
    }

    //! Function to get partial derivative object of direct acceleration from acceleratingBody on acceleratedBody.
    /*!
     * Function to get the partial derivative object of direct acceleration from acceleratingBody on acceleratedBody.
     * \return Partial derivative object of direct acceleration from acceleratingBody on acceleratedBody.
     */
    std::shared_ptr< DirectGravityPartial > getPartialOfDirectGravityOnCentralBody( )
    {
        return partialOfDirectGravityOnCentralBody_;
    }

    //! Function to get the name of central body w.r.t. which the acceleration is computed.
    /*!
     * Function to get the name of central body w.r.t. which the acceleration is computed.
     * \return Name of central body w.r.t. which the acceleration is computed.
     */
    std::string getCentralBodyName( )
    {
        return centralBodyName_;
    }


protected:

    //! Function to reset the constituent DirectGravityPartial objects to the current time.
    void resetTimeOfMemberObjects( )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->resetTime( currentTime_ );
        partialOfDirectGravityOnCentralBody_->resetTime( currentTime_ );
    }

    //! Function to update the parameter partials of the constituent DirectGravityPartial objects.
    void updateParameterPartialsOfMemberObjects( )
    {
        partialOfDirectGravityOnBodyUndergoingAcceleration_->updateParameterPartials( );
        partialOfDirectGravityOnCentralBody_->updateParameterPartials( );
    }

private:
    //! Partial derivative object of direct acceleration from acceleratingBody on acceleratedBody.
    std::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnBodyUndergoingAcceleration_;

    //! Partial derivative object of direct acceleration from centralBodyName on acceleratedBody.
    std::shared_ptr< DirectGravityPartial > partialOfDirectGravityOnCentralBody_;

    //! Name of central body w.r.t. which the acceleration is computed.
    std::string centralBodyName_;

};

//! Function to retrieve name of central body of third-body acceleration partial
/*!
 *  Function to retrieve name of central body of third-body acceleration partial, from AccelerationPartial base class object input
 *  \param accelerationPartial Acceleration partial model from which central body name is to be retrieved (must be of derived
 *  class ThirdBodyGravityPartial< T >
 *  \return Name of central body of third-body acceleration partial
 */
inline std::string getCentralBodyNameFromThirdBodyAccelerationPartial(
        const std::shared_ptr< AccelerationPartial > accelerationPartial )
{
    std::string centralBody;
    if( !basic_astrodynamics::isAccelerationFromThirdBody( accelerationPartial->getAccelerationType( ) ) )
    {
        throw std::runtime_error( "Error, requested third body from acceleration partial, but input is incompatible." );
    }
    else
    {
        if( accelerationPartial->getAccelerationType( ) == basic_astrodynamics::third_body_point_mass_gravity )
        {
            centralBody = std::dynamic_pointer_cast< ThirdBodyGravityPartial< CentralGravitationPartial > >(
                        accelerationPartial )->getCentralBodyName( );
        }
        else if( accelerationPartial->getAccelerationType( ) == basic_astrodynamics::third_body_spherical_harmonic_gravity )
        {
            centralBody = std::dynamic_pointer_cast< ThirdBodyGravityPartial< SphericalHarmonicsGravityPartial > >(
                        accelerationPartial )->getCentralBodyName( );
        }
        else if( accelerationPartial->getAccelerationType( ) == basic_astrodynamics::third_body_mutual_spherical_harmonic_gravity )
        {
            centralBody = std::dynamic_pointer_cast< ThirdBodyGravityPartial< MutualSphericalHarmonicsGravityPartial > >(
                        accelerationPartial )->getCentralBodyName( );
        }
        else
        {
            throw std::runtime_error( "Error, requested third body from acceleration partial, input third-body type not recognized." );
        }
    }

    return centralBody;
}

} // namespace acceleration_partials

} // namespace tudat


#endif // TUDAT_THIRDBODYGRAVITYPARTIAL_H
