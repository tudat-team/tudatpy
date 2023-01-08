/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALHARMONICACCELERATIONPARTIAL_H
#define TUDAT_SPHERICALHARMONICACCELERATIONPARTIAL_H

#include "tudat/math/basic/coordinateConversions.h"

#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"

#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/sphericalHarmonicAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/tidalLoveNumberPartialInterface.h"
#include "tudat/astro/orbit_determination/observation_partials/rotationMatrixPartial.h"



namespace tudat
{

namespace acceleration_partials
{

//! Class for calculating partial derivatives of a spherical harmonic gravitational acceleration.
/*!
 *  Class for calculating partial derivatives of a spherical harmonic gravitational acceleration, as calculated by the
 *  SphericalHarmonicsGravitationalAccelerationModel class.
 */
class SphericalHarmonicsGravityPartial: public AccelerationPartial
{
public:

    //! Contructor.
    /*!
     *  Constructor, requires input on the acceleration model as of which partials are to be computed.
     *  If any partials of parameters of the rotation model of the body exerting acceleration are to be calculated,
     *  RotationMatrixPartial objects must be pre-constructed and passed here as a map, with one object for each parameter
     *  wrt which a partial is to be taken.
     *  \param acceleratedBody Name of body undergoing acceleration.
     *  \param acceleratingBody Name of body exerting acceleration.
     *  \param accelerationModel Spherical harmonic gravity acceleration model from which acceleration is calculated wrt
     *  which the object being constructed is to calculate partials.
     *  \param rotationMatrixPartials Map of RotationMatrixPartial, one for each paramater representing a property of the
     *  rotation of the body exerting the acceleration wrt which an acceleration partial will be calculated.
     *  \param tidalLoveNumberPartialInterfaces List of objects to compute partials of tidal gravity field variations, one
     *  per corresponding variation object in acceleratedBody.
     */
    SphericalHarmonicsGravityPartial(
            const std::string& acceleratedBody,
            const std::string& acceleratingBody,
            const std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > accelerationModel,
            const observation_partials::RotationMatrixPartialNamedList& rotationMatrixPartials =
            observation_partials::RotationMatrixPartialNamedList( ),
            const std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >&
            tidalLoveNumberPartialInterfaces =
            std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > >( ) );

    //! Destructor
    ~SphericalHarmonicsGravityPartial( ){ }

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
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration..
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body undergoing acceleration
     *  and adding it to the existing partial block
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 3 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
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
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
    }

    //! Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the velocity of body exerting acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian velocity of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtVelocity_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtVelocity_;
        }
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
        bool doesDependencyExist = false;
        if( ( ( stateReferencePoint.first == acceleratingBody_ ||
              ( stateReferencePoint.first == acceleratedBody_  && accelerationUsesMutualAttraction_ ) )
              && integratedStateType == propagators::body_mass_state ) )
        {
            throw std::runtime_error( "Warning, dependency of central gravity on body masses not yet implemented" );
        }
        else if( stateReferencePoint.first == acceleratingBody_ && integratedStateType == propagators::rotational_state )
        {
            doesDependencyExist = true;
        }
        return doesDependencyExist;
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

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
            const bool addContribution = true )
    {
        if( stateReferencePoint.first == acceleratingBody_ && integratedStateType == propagators::rotational_state )
        {
            Eigen::MatrixXd tempMatrix = Eigen::MatrixXd::Zero( 3, 7 );
            wrtRotationModelParameter( tempMatrix, estimatable_parameters::initial_rotational_body_state, "" );
            partialMatrix.block( 0, 0, 3, 7 ) = ( addContribution ? 1.0 : -1.0 ) * tempMatrix;
        }
    }

    //! Function to create a function returning the current partial w.r.t. a gravitational parameter.
    /*!
     * Function to create a function returning the current partial w.r.t. a gravitational parameter.
     * \param parameterId Identifier of parameter for which the partial is to be created.
     * \return Pair with partial function and paramater partial size. The partial function is non-empty only
     * if the parameterId input represents the gravitational parameter of acceleratingBody_ (or acceleratedBody_ if
     * accelerationUsesMutualAttraction_ is true).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getGravitationalParameterPartialFunction(
            const estimatable_parameters::EstimatebleParameterIdentifier& parameterId );

    //! Function for updating the partial object to current state and time.
    /*!
     *  Function for updating the partial object to current state and time. Calculates the variables that are
     *  used for the calculation of multple partials, to prevent multiple calculations of same function.
     *  \param currentTime Time to which object is to be updated (note that most update functions are time-independent,
     *  since the 'current' state of the bodies is typically updated globally by the NBodyStateDerivative class).
     */
    virtual void update( const double currentTime = TUDAT_NAN );

    //! Function to calculate the partial wrt the gravitational parameter.
    /*!
     *  Function to calculate the partial wrt the gravitational parameter of the central body. Note that in the case of
     *  mutual attraction (see SphericalHarmonicsGravitationalAccelerationModel), the partial wrt the gravitational
     *  parameter of the body exerting acceleration is equal to the partial wrt to the body undergoing the acceleration,
     *  which is zero otherwise.
     *  \return Partial wrt the gravitational parameter of the central body.
     */
    void wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& partialMatrix, const int addPartial = 0 )
    {
        if( gravitationalParameterFunction_( ) != 0.0 )
        {
            if( addPartial == 0 )
            {
                partialMatrix = accelerationFunction_( ) / gravitationalParameterFunction_( );
            }
            else if( addPartial == 1 )
            {
                partialMatrix += accelerationFunction_( ) / gravitationalParameterFunction_( );
            }
            else if( addPartial == -1 )
            {
                partialMatrix -= accelerationFunction_( ) / gravitationalParameterFunction_( );
            }
            else
            {
                throw std::runtime_error( "Error when adding partial of spherical harmonic acceleration w.r.t mu, inout is inconsistent" );
            }
        }
        else
        {
            throw std::runtime_error( "Error cannot compute partial of spherical harminic gravity w.r.t mu for zero value" );
        }
    }

    //! Function to retrieve partial of acceleration wrt the position of body undergoing acceleration, in inertial coordinates.
    /*!
     * Function to retrieve the current partial of the acceleration wrt the position of the body undergoing the acceleration,
     * in inertial coordinates
     * \return Current partial of the acceleration wrt the position of the body undergoing the acceleration, in inertial coordinates.
     */
    Eigen::Matrix3d getCurrentPartialWrtPosition( )
    {
        return currentPartialWrtPosition_;
    }

    //! Function to retrieve partial of acceleration wrt the position of body undergoing acceleration, in body-fixed coordinates.
    /*!
     * Function to retrieve the current partial of the acceleration wrt the position of the body undergoing the acceleration,
     * in body-fixed coordinates
     * \return Current partial of the acceleration wrt the position of the body undergoing the acceleration, in body-fixed coordinates.
     */
    Eigen::Matrix3d getCurrentBodyFixedPartialWrtPosition( )
    {
        return currentBodyFixedPartialWrtPosition_;
    }

protected:

    //! Function to reset the relevant member objects to the current time.
    void resetCurrentTimeOfMemberObjects( )
    {
        for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        {
            tidalLoveNumberPartialInterfaces_.at( i )->resetCurrentTime( );
        }
    }

    //! Function to update the parameter partials of the relevant member objects to the current time.
    void updateParameterPartialsOfMemberObjects( )
    {
        for( unsigned int i = 0; i < tidalLoveNumberPartialInterfaces_.size( ); i++ )
        {
            tidalLoveNumberPartialInterfaces_.at( i )->updateParameterPartials( );
        }
    }

    //! Function to calculate the partial of the acceleration wrt a set of cosine coefficients.
    /*!
     *  Function to calculate the partial of the acceleration wrt a set of cosine coefficients.
     *  The set of coefficients wrt which a partial is to be taken is provided as input.
     *  \param blockIndices List of cosine coefficient indices wrt which the partials are to be taken (first and second
     *  are degree and order for each vector entry).
     *  \param partialDerivatives Matrix of acceleration partials that is set by this function (returned by reference),
     *  with each column containg the partial wrt a single coefficient (in same order as blockIndices).
     */
    void wrtCosineCoefficientBlock(
            const std::vector< std::pair< int, int > >& blockIndices,
            Eigen::MatrixXd& partialDerivatives );

    //! Function to calculate the partial of the acceleration wrt a set of sine coefficients.
    /*!
     *  Function to calculate the partial of the acceleration wrt a set of sine coefficients.
     *  The set of coefficients wrt which a partial is to be taken is provided as input.
     *  \param blockIndices List of sine coefficient indices wrt which the partials are to be taken (first and second
     *  are degree and order for each vector entry).
     *  \param partialDerivatives Matrix of acceleration partials that is set by this function (returned by reference),
     *  with each column containg the partial wrt a single coefficient (in same order as blockIndices).
     */
    void wrtSineCoefficientBlock(
            const std::vector< std::pair< int, int > >& blockIndices,
            Eigen::MatrixXd& partialDerivatives );

    //! Function to calculate an acceleration partial wrt a rotational parameter.
    /*!
     *  Function to calculate an acceleration partial wrt a rotational parameter of the rotation model of the body
     *  exerting the acceleration.
     *  \param accelerationPartial Matrix of partials of spherical harmonic acceleration wrt a rotational parameter
     *  that is set by this function (returned by reference)
     *  \param parameterType Type of parameter wrt which a partial is to be calculated.
     *  An entry of the requested type must be present in the rotationMatrixPartials_ map.
     *  \param secondaryIdentifier Identifier required to unambiguously define the parameter (in addition to information in
     *  parameterType.
     */
    void wrtRotationModelParameter(
            Eigen::MatrixXd& accelerationPartial,
            const estimatable_parameters::EstimatebleParametersEnum parameterType,
            const std::string& secondaryIdentifier );

    //! Function to calculate an acceleration partial wrt a tidal parameter.
    /*!
     *  Function to calculate an acceleration partial wrt a tidal parameter of the gravitu field of the body
     *  exerting the acceleration.
     *  \param coefficientPartialFunctions Function returning a list of partial derivatives of C,S spherical harmonic coefficients
     *  with the degree and order determined by corresponding input arguments of this function. Degree is fixed, and this
     *  vector order corresponds to order of 'orders' input parameter (one vector entry per order). Matrix rows denote
     *  C and S coefficient, respectively. Columns denote entries of tidal property (e.g. first entry for real Love number
     *  component, second for complex part)
     * \param degree Degree of property (e.g. Love number) w.r.t. which partials are to be computed
     * \param orders Orders of property (e.g. Love number) w.r.t. which partials are to be computed
     * \param sumOrders True of the contributions of the various orders are to be summed (i.e. assumed to be constant for all
     * orders at given degree), or if they are handled separately
     * \param parameterSize Size of the parameter w.r.t. which the partials are to be computed
     * \param accelerationPartial Matrix of partials of spherical harmonic acceleration wrt a tidal parameter
     * (returned by reference)
     */
    void wrtTidalModelParameter(
            const std::function< std::vector< Eigen::Matrix< double, 2, Eigen::Dynamic > >( ) > coefficientPartialFunctions,
            const int degree,
            const std::vector< int >& orders,
            const bool sumOrders,
            const int parameterSize,
            Eigen::MatrixXd& accelerationPartial );

    //! Function to return the gravitational parameter used for calculating the acceleration.
    std::function< double( ) > gravitationalParameterFunction_;

    //! Function to return the reference radius used for calculating the acceleration.
    std::function< double( ) > bodyReferenceRadius_;

    //! Function to return the current cosine coefficients of the spherical harmonic gravity field.
    std::function< Eigen::MatrixXd( ) > cosineCoefficients_;

    //! Function to return the current sine coefficients of the spherical harmonic gravity field.
    std::function< Eigen::MatrixXd( ) > sineCoefficients_;

    //! Cache object used for storing calculated values at current time and state for spherical harmonic gravity
    //! calculations.
    std::shared_ptr< basic_mathematics::SphericalHarmonicsCache > sphericalHarmonicCache_;

    //! Function returning position of body undergoing acceleration.
    std::function< Eigen::Vector3d( ) > positionFunctionOfAcceleratedBody_;

    //! Function returning position of body exerting acceleration.
    std::function< Eigen::Vector3d( ) > positionFunctionOfAcceleratingBody_;

    //! Function return current rotation from inertial frame to frame fixed to body exerting acceleration.
    std::function< Eigen::Matrix3d( ) > fromBodyFixedToIntegrationFrameRotation_;

    //! Function to retrieve the current spherical harmonic acceleration.
    std::function< Eigen::Matrix< double, 3, 1 >( ) > accelerationFunction_;

    //! Function to update the acceleration to the current state and time.
    /*!
     *  Function to update the acceleration to the current state and time.
     *  Called when updating an object of this class with the update( time ) function,
     *  in case the partial is called before the acceleration model in the current iteration of the numerical integration.
     */
    std::function< void( const double ) > updateFunction_;

    //! Current cosine coefficients of the spherical harmonic gravity field.
    /*!
     *  Current cosine coefficients of the spherical harmonic gravity field, set by update( time ) function.
     */
    Eigen::MatrixXd currentCosineCoefficients_;

    //! Current sine coefficients of the spherical harmonic gravity field.
    /*!
     *  Current sine coefficients of the spherical harmonic gravity field, set by update( time ) function.
     */
    Eigen::MatrixXd currentSineCoefficients_;

    //! Current body-fixed (w.r.t body exerting acceleration) position of body undergoing acceleration
    /*!
     *  Current body-fixed (w.r.t body exerting acceleration) position of body undergoing acceleration,
     *  set by update( time ) function.
     */
    Eigen::Vector3d bodyFixedPosition_;

    //! Current spherical coordinate of body undergoing acceleration
    /*!
     *  Current spherical coordinate of body undergoing acceleration, in reference frame fixed to body exerting acceleration.
     *  Order of components is radial distance (from center of body), latitude, longitude. Note that the the second entry
     *  differs from the direct output of the cartesian -> spherical coordinates, which produces a colatitude.
     */
    Eigen::Vector3d bodyFixedSphericalPosition_;

    //! The current partial of the acceleration wrt the position of the body undergoing the acceleration.
    /*!
     *  The current partial of the acceleration wrt the position of the body undergoing the acceleration.
     *  The partial wrt the position of the body exerting the acceleration is minus this value.
     *  Value is set by the update( time ) function.
     */
    Eigen::Matrix3d currentPartialWrtPosition_;

    //! The current partial of the acceleration wrt the position of the body undergoing the acceleration,
    //! with both acceleration and position in body-fixed frame.
    /*!
     *  The current partial of the acceleration wrt the position of the body undergoing the acceleration.
     *  with both acceleration and position in body-fixed frame. Value is set by the update( time ) function.
     */
    Eigen::Matrix3d currentBodyFixedPartialWrtPosition_;

    //! The current partial of the acceleration wrt the velocity of the body undergoing the acceleration.
    /*!
     *  The current partial of the acceleration wrt the velocity of the body undergoing the acceleration.
     *  The partial wrt the velocity of the body exerting the acceleration is minus this value.
     * Value is set by the update( time ) function.
     */
    Eigen::Matrix3d currentPartialWrtVelocity_;

    //! Maximum degree of spherical harmonic expansion.
    /*!
     *  Maximum degree of spherical harmonic expansion of body exerting acceleration used in the calculation
     *  of the acceleration.
     */
    int maximumDegree_;

    //! Maximum order of spherical harmonic expansion.
    /*!
     *  Maximum order of spherical harmonic expansion of body exerting acceleration used in the calculation
     *  of the acceleration.
     */
    int maximumOrder_;

    //! Map of RotationMatrixPartial, one for each relevant rotation parameter
    /*!
     *  Map of RotationMatrixPartial, one for each parameter representing a property of the rotation of the
     *  body exerting the acceleration wrt which an acceleration partial will be calculated.
     *  Map is pre-created and set through the constructor.
     */
    observation_partials::RotationMatrixPartialNamedList rotationMatrixPartials_;

    //! List of objects to compute partials of tidal gravity field variations
    /*!
     * List of objects to compute partials of tidal gravity field variations, one per corresponding variation object in
     * acceleratedBody.
     */
    std::vector< std::shared_ptr< orbit_determination::TidalLoveNumberPartialInterface > > tidalLoveNumberPartialInterfaces_;

    //! Boolean denoting whether the mutual attraction between the bodies is taken into account
    /*!
     *  Boolean denoting whether the mutual attraction between the bodies is taken into account, as is the case when
     *  integrting in a (non-rotating) frame centered on the body exerting the acceleration, in which case the
     *  gravitational acceleration of the body undergoing the acceleration on that exerting the acceleration must be taken
     *  into account as an inertial effect.
     */
    bool accelerationUsesMutualAttraction_;

};

} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_SPHERICALHARMONICACCELERATIONPARTIAL_H
