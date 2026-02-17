/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_SOLAR_SYSTEM_METRIC_H
#define TUDAT_SOLAR_SYSTEM_METRIC_H

#include "tudat/astro/gravitation/sphericalHarmonicsGravityModel.h"
#include "tudat/math/basic/linearAlgebra.h"

#include "tudat/astro/gravitation/sphericalHarmonicsGravityField.h"
#include "tudat/astro/relativity/relativisticPotentials.h"
#include "tudat/astro/relativity/metric.h"
#include "tudat/simulation/environment_setup/body.h"

#include <functional>
#include <map>
#include <memory>
#include <vector>
#include <iostream>

namespace tudat
{

class SphericalHarmonicWrapper
{
public:

    //! Constructor.
    /*!
     *  \param sphericalHarmonicGravityField Spherical-harmonic gravity field model.
     *  \param centralBodyPositionFunction Function returning central-body position.
     *  \param maximumDegreeAndOrder Maximum degree/order used in spherical-harmonic evaluations.
     *  \param legendreCache Cache for Legendre polynomial evaluations.
     *  \param currentRotationToBodyFixedFrameDerivative Function returning body-frame rotation derivative.
     */
    SphericalHarmonicWrapper(
            const std::shared_ptr< gravitation::SphericalHarmonicsGravityField > sphericalHarmonicGravityField,
            const std::function< Eigen::Vector3d( ) > centralBodyPositionFunction,
            const std::pair< int, int >& maximumDegreeAndOrder,
            basic_mathematics::LegendreCache* legendreCache,
            const std::function< Eigen::Matrix3d( ) > currentRotationToBodyFixedFrameDerivative ):
        sphericalHarmonicGravityField_( sphericalHarmonicGravityField ),
        currentRotationToBodyFixedFrameDerivative_( currentRotationToBodyFixedFrameDerivative ),
        centralBodyPositionFunction_( centralBodyPositionFunction ),
        maximumDegreeAndOrder_( maximumDegreeAndOrder ),
        legendreCache_( legendreCache )
    {
        using namespace std::placeholders;
        
        gravitation::SphericalHarmonicsBlock cosineBlock(
            std::bind( &gravitation::SphericalHarmonicsGravityField::getCosineCoefficientsReference, sphericalHarmonicGravityField_ ),
            maximumDegreeAndOrder_.first + 1,
            maximumDegreeAndOrder_.second + 1 );

        gravitation::SphericalHarmonicsBlock sineBlock(
            std::bind( &gravitation::SphericalHarmonicsGravityField::getSineCoefficientsReference, sphericalHarmonicGravityField_ ),
            maximumDegreeAndOrder_.first + 1,
            maximumDegreeAndOrder_.second + 1 );

        // Wrap the position functions for the acceleration model (expects output-by-reference)
        std::function< void( Eigen::Vector3d& ) > wrappedCentralBodyPositionFunction =
            [=]( Eigen::Vector3d& output ){ output = centralBodyPositionFunction_( ); };
        std::function< void( Eigen::Vector3d& ) > wrappedEvaluationPositionFunction =
            [this]( Eigen::Vector3d& output ){ output = this->getCurrentEvaluationPosition( ); };

        // Create Spherical Harmonics cache
        basic_mathematics::SphericalHarmonicsCache harmonicsCache( maximumDegreeAndOrder_.first + 1,
                                                                    maximumDegreeAndOrder_.second + 1 );

        sphericalHarmonicPotentialGradientModel_ =
            std::make_shared< gravitation::SphericalHarmonicsGravitationalAccelerationModel >(
                wrappedEvaluationPositionFunction,
                std::bind( &gravitation::GravityFieldModel::getGravitationalParameter, sphericalHarmonicGravityField_ ),
                sphericalHarmonicGravityField_->getReferenceRadius( ),
                cosineBlock,
                sineBlock,
                wrappedCentralBodyPositionFunction,
                std::bind( &reference_frames::RotationWrapper::getRotationQuaternion,
                        sphericalHarmonicGravityField_->getRotationToLocalFrameWrapper( ) ),
                false,
                basic_mathematics::SphericalHarmonicsCache( maximumDegreeAndOrder_.first + 1,
                                                            maximumDegreeAndOrder_.second + 1 ) );


    }

    //! Compute current spherical-harmonic potential at evaluation position.
    /*!
     *  \param currentEvaluationPosition Position at which potential is evaluated.
     *  \return Gravitational potential value.
     */
    double getCurrentPotential( const Eigen::Vector3d& currentEvaluationPosition )
    {
        currentEvaluationPosition_ = currentEvaluationPosition;
        return sphericalHarmonicGravityField_->getGravitationalPotentialFromInertialPosition(
                    currentEvaluationPosition_ - centralBodyPositionFunction_( ),
                    maximumDegreeAndOrder_.first, maximumDegreeAndOrder_.second,
                    legendreCache_ );
    }

    //! Compute potential gradient using spherical-harmonic acceleration model.
    /*!
     *  \param currentEvaluationPosition Position at which gradient is evaluated.
     *  \param currentTime Current evaluation time.
     *  \return Potential gradient vector.
     */
    Eigen::Vector3d getSphericalHarmonicPotentialGradient(
            const Eigen::Vector3d& currentEvaluationPosition, const double currentTime = 0.0 )
    {
        currentEvaluationPosition_ = currentEvaluationPosition;
        sphericalHarmonicPotentialGradientModel_->updateMembers( currentTime );
        return sphericalHarmonicPotentialGradientModel_->getAcceleration( );
    }

    //! Get currently stored evaluation position.
    /*!
     *  \return Current evaluation position.
     */
    Eigen::Vector3d getCurrentEvaluationPosition( )
    {
        return currentEvaluationPosition_;
    }

    //! Get spherical-harmonic acceleration model object.
    /*!
     *  \return Shared pointer to spherical-harmonic acceleration model.
     */
    std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > getSphericalHarmonicPotentialGradientModel( )
    {
        return sphericalHarmonicPotentialGradientModel_;
    }

    //! Get Legendre cache.
    /*!
     *  \return Pointer to Legendre cache.
     */
    basic_mathematics::LegendreCache* getLegendreCache( )
    {
        return legendreCache_;
    }

    //! Get rotation from inertial frame to body-fixed frame.
    /*!
     *  \return Rotation quaternion.
     */
    Eigen::Quaterniond getRotationToBodyFixedFrame( )
    {
        return sphericalHarmonicGravityField_->getRotationToLocalFrameWrapper( )->getRotationQuaternion( );
    }

    //! Get time derivative of rotation to body-fixed frame.
    /*!
     *  \return 3x3 rotation-derivative matrix.
     */
    Eigen::Matrix3d getRotationDerivativeToBodyFixedFrame( )
    {
        return currentRotationToBodyFixedFrameDerivative_( );
    }

    //! Get underlying spherical-harmonic gravity field.
    /*!
     *  \return Shared pointer to gravity field.
     */
    std::shared_ptr< gravitation::SphericalHarmonicsGravityField > getSphericalHarmonicGravityField( )
    {
        return sphericalHarmonicGravityField_;
    }

private:

    std::shared_ptr< gravitation::SphericalHarmonicsGravityField > sphericalHarmonicGravityField_;

    std::shared_ptr< gravitation::SphericalHarmonicsGravitationalAccelerationModel > sphericalHarmonicPotentialGradientModel_;

    std::function< Eigen::Matrix3d( ) > currentRotationToBodyFixedFrameDerivative_;

    std::function< Eigen::Vector3d( ) > centralBodyPositionFunction_;

    Eigen::Vector3d currentEvaluationPosition_;

    std::pair< int, int > maximumDegreeAndOrder_;

    basic_mathematics::LegendreCache* legendreCache_;
};

namespace relativity
{

//! Post-Newtonian solar-system metric model with optional higher-order terms.
class SolarSystemMetric: public Metric
{
public:
    //! Constructor.
    /*!
     *  \param bodyList Names of bodies included in metric evaluation.
     *  \param bodyGravitationalParameterFunctions Functions returning body gravitational parameters.
     *  \param bodyStateFunctions Functions returning body translational states.
     *  \param secondOrderBodyList Indices of bodies included in second-order terms.
     *  \param ppnParameterSet PPN parameter set used by the model.
     *  \param bodyAccelerationFunctions Optional body-acceleration functions.
     *  \param bodyAngularMomentumFunctions Optional body angular-momentum functions keyed by body index.
     *  \param bodySphericalHarmonicGravityWrappers Optional spherical-harmonic wrappers keyed by body index.
     *  \param rotationUpdateFunctions Optional rotation-update functions keyed by body index.
     */
    SolarSystemMetric(
            const std::vector< std::string >& bodyList,
            const std::vector< std::function< double( ) > >& bodyGravitationalParameterFunctions,
            const std::vector< std::function< Eigen::Vector6d( ) > >& bodyStateFunctions,
            const std::vector< int >& secondOrderBodyList,
            const std::shared_ptr< PPNParameterSet > ppnParameterSet,
            const std::vector< std::function< Eigen::Vector3d( const double ) > > bodyAccelerationFunctions =
            ( std::vector< std::function< Eigen::Vector3d( const double ) > >( ) ),
            const std::map< int, std::function< double( ) > >& bodyAngularMomentumFunctions =
            ( std::map< int, std::function< double( ) > >( ) ),
            const std::map< int, std::shared_ptr< SphericalHarmonicWrapper > > bodySphericalHarmonicGravityWrappers =
            ( std::map< int, std::shared_ptr< SphericalHarmonicWrapper > >( ) ),
            const std::map< int, std::function< void( const double ) > > rotationUpdateFunctions =
            ( std::map< int, std::function< void( const double ) > >( ) ) ):
        Metric( ),
        bodyList_( bodyList ),
        bodyGravitationalParameterFunctions_( bodyGravitationalParameterFunctions ),
        bodyStateFunctions_( bodyStateFunctions ),
        secondOrderBodyList_( secondOrderBodyList ),
        bodyAccelerationFunctions_( bodyAccelerationFunctions ),
        bodyAngularMomentumFunctions_( bodyAngularMomentumFunctions ),
        bodySphericalHarmonicGravityWrappers_( bodySphericalHarmonicGravityWrappers ),
        rotationUpdateFunctions_( rotationUpdateFunctions ),
        ppnParameterSet_( ppnParameterSet )
    {
        currentBodyGravitationalParameters_.resize( bodyList_.size( ) );
        currentBodyStates_.resize( bodyList_.size( ) );
        currentBodyDistances_.resize( bodyList_.size( ) );
        currentBodyRelativePositions_.resize( bodyList_.size( ) );
        currentScalarPotentialGradients_.resize( bodyList_.size( ) );
        currentBodyAccelerations_.resize( bodyList_.size( ) );

        if( bodyAccelerationFunctions_.size( ) == 0 )
        {
            for( unsigned int i = 0; i < bodyList.size( ); i++ )
            {
                bodyAccelerationFunctions_.push_back( []( const double ){ return Eigen::Vector3d::Zero(); } );
            }
        }
    }

    //! Copy constructor.
    /*!
     *  \param originalMetric Original metric object.
     */
    SolarSystemMetric( const SolarSystemMetric& originalMetric )
    {
        bodyList_= originalMetric.bodyList_;
        bodyGravitationalParameterFunctions_ = originalMetric.bodyGravitationalParameterFunctions_;
        bodyStateFunctions_ = originalMetric.bodyStateFunctions_;
        secondOrderBodyList_ = originalMetric.secondOrderBodyList_;
        ppnParameterSet_ = originalMetric.ppnParameterSet_;
        bodyAngularMomentumFunctions_ = originalMetric.bodyAngularMomentumFunctions_;
        bodySphericalHarmonicGravityWrappers_ = originalMetric.bodySphericalHarmonicGravityWrappers_;
        rotationUpdateFunctions_ = originalMetric.rotationUpdateFunctions_;
        bodyAccelerationFunctions_ = originalMetric.bodyAccelerationFunctions_;
        bodySphericalHarmonicPotentialFunctions_ = originalMetric.bodySphericalHarmonicPotentialFunctions_;

        currentScalarPotentialGradients_ = originalMetric.currentScalarPotentialGradients_;
        currentSphericalHarmonicPotentials_ = originalMetric.currentSphericalHarmonicPotentials_;
        currentExternalScalarPotentials_ = originalMetric.currentExternalScalarPotentials_;
        currentBodyAccelerations_= originalMetric.currentBodyAccelerations_;
        currentBodyGravitationalParameters_.resize( bodyList_.size( ) );
        currentBodyStates_.resize( bodyList_.size( ) );
        currentBodyDistances_.resize( bodyList_.size( ) );
        currentBodyRelativePositions_.resize( bodyList_.size( ) );
    }

    ~SolarSystemMetric( ){ }

    //! Clone metric object.
    /*!
     *  \return Shared pointer to copied metric object.
     */
    std::shared_ptr< Metric > Clone( )
    {
        return std::make_shared< SolarSystemMetric >( *this );
    }

    //! Update metric and/or Christoffel symbols at requested state and time.
    /*!
     *  \param state Current translational state at evaluation point.
     *  \param time Current evaluation time.
     *  \param updateCurrentMetric If true, update covariant metric.
     *  \param updateCurrentChristoffelSymbols If true, update Christoffel symbols.
     */
    void update( const Eigen::Matrix< double, 6, 1 >& state, const double time,
                 const bool updateCurrentMetric, const bool updateCurrentChristoffelSymbols );

    //! Get current Christoffel symbols.
    /*!
     *  \return Vector of Christoffel-symbol matrices.
     */
    std::vector< Eigen::Matrix< double, 4, 4 > > getCurrentChristoffelSymbols( )
    {
        std::cerr<<"Error, solar system metric chirstoffel symbols not yet implemented"<<std::endl;
        return std::vector< Eigen::Matrix< double, 4, 4 > >( );
    }

    //! Get list of bodies used in metric model.
    /*!
     *  \return Body-name list.
     */
    std::vector< std::string > getBodyList( ) { return bodyList_; }

    //! Get PPN parameter set.
    /*!
     *  \return PPN parameter set.
     */
    std::shared_ptr< PPNParameterSet > getPpnParameterSet( ) { return ppnParameterSet_; }

    //! Get spherical-harmonic wrappers used by metric model.
    /*!
     *  \return Map of spherical-harmonic wrappers keyed by body index.
     */
    std::map< int, std::shared_ptr< SphericalHarmonicWrapper > > getBodySphericalHarmonicGravityWrappers( ) { return bodySphericalHarmonicGravityWrappers_; }

    //! Get current state of one metric body.
    /*!
     *  \param bodyIndex Index of body in metric body list.
     *  \return Current body state.
     */
    Eigen::Vector6d getCurrentBodyStates( const int bodyIndex ) { return currentBodyStates_.at( bodyIndex ); }

    //! Get current vector potential.
    /*!
     *  \return Current vector potential.
     */
    Eigen::Vector3d getCurrentVectorPotential( ) { return currentVectorPotential_; }

    //! Get current scalar potential contribution from a single body.
    /*!
     *  \param bodyIndex Index of body in metric body list.
     *  \return Scalar potential contribution.
     */
    double getCurrentSingleBodyScalarPotential( const int bodyIndex ) { return currentExternalScalarPotentials_.at( bodyIndex ); }

    //! Get total current scalar potential.
    /*!
     *  \return Total scalar potential.
     */
    double getCurrentScalarPotential( ) { return currentTotalExternalScalarPotential_; }

    //! Get current squared scalar-potential contribution terms.
    /*!
     *  \return Squared scalar-potential terms.
     */
    double getCurrentScalarPotentialSquaredTerms( ) { return currentTotalExternalScalarPotentialSquaredTerms_; }

    //! Get current gravitational parameter for one body.
    /*!
     *  \param bodyIndex Index of body in metric body list.
     *  \return Gravitational parameter value.
     */
    double getCurrentGravitationalParameter( const int bodyIndex ) { return currentBodyGravitationalParameters_.at( bodyIndex ); }

    //! Get current distance to one body.
    /*!
     *  \param bodyIndex Index of body in metric body list.
     *  \return Distance from evaluation point to body.
     */
    double getCurrentBodyDistance( const int bodyIndex ) { return currentBodyDistances_.at( bodyIndex ); }

    //! Get current relative position to one body.
    /*!
     *  \param bodyIndex Index of body in metric body list.
     *  \return Relative position vector.
     */
    Eigen::Vector3d getCurrentBodyRelativePosition( const int bodyIndex ) { return currentBodyRelativePositions_.at( bodyIndex ); }

    //! Get current scalar-potential gradient contribution from one body.
    /*!
     *  \param bodyIndex Index of body in metric body list.
     *  \return Scalar-potential gradient vector.
     */
    Eigen::Vector3d getCurrentScalarPotentialGradient( const int bodyIndex ) { return currentScalarPotentialGradients_.at( bodyIndex ); }

    //! Get time partial of scalar potential.
    /*!
     *  \return Scalar-potential time partial.
     */
    double getCurrentScalarPotentialTimePartial( ) { return currentScalarPotentialTimePartial_; }

    //! Get position partial of scalar potential.
    /*!
     *  \return Scalar-potential position partial.
     */
    Eigen::Vector3d getCurrentScalarPotentialPositionPartial( ) { return currentScalarPotentialPositionPartial_; }

    //! Get time partial of vector potential.
    /*!
     *  \return Vector-potential time partial.
     */
    Eigen::Vector3d getCurrentVectorPotentialTimePartial( ) { return currentVectorPotentialTimePartial_; }

    //! Get inherent (non-kinematic) scalar-potential time partial.
    /*!
     *  \return Inherent scalar-potential time partial.
     */
    double getCurrentInherentScalarPotentialTimePartial( ) { return currentInherentScalarPotentialTimePartial_; }

    //! Get inherent (non-kinematic) vector-potential time partial.
    /*!
     *  \return Inherent vector-potential time partial.
     */
    Eigen::Vector3d getCurrentInherentVectorPotentialTimePartial( ) { return currentInherentVectorPotentialTimePartial_; }

    //! Get position partial of vector potential.
    /*!
     *  \return Vector-potential position-partial matrix.
     */
    Eigen::Matrix3d getCurrentVectorPotentialPositionPartial( ) { return currentVectorPotentialPositionPartial_; }

    //! Get current acceleration of one metric body.
    /*!
     *  \param bodyIndex Index of body in metric body list.
     *  \return Body acceleration vector.
     */
    Eigen::Vector3d getCurrentBodyAccelerations( const int bodyIndex ) { return currentBodyAccelerations_.at( bodyIndex ); }

protected:

    void updateCurrentState( const bool updateEvaluationPointIndependentVariables );

    void updateMetric( );

    void updateChristoffelSymbols( );

    void updatePotentialPartials( );

    std::vector< double > calculateScalarPotentials( )
    {
        return calculateFirstOrderExternalScalarPotentials(
                    currentBodyGravitationalParameters_, currentBodyDistances_, currentSphericalHarmonicPotentials_ );
    }

    Eigen::Vector3d calculateVectorPotential( )
    {
        return calculateExternalVectorPotential(
                    currentBodyGravitationalParameters_, currentBodyDistances_, currentBodyStates_ );
    }

    Eigen::Vector3d calculateVectorPotential( const std::vector< double >& scalarPotentials )
    {
        return calculateExternalVectorPotential(
                    scalarPotentials, currentBodyStates_ );
    }

    Eigen::Matrix3d calculateAnisotropicMatrixPotential( )
    {
        Eigen::Matrix3d anisotropicPotential = Eigen::Matrix3d::Zero( );
        for( unsigned int i = 0; i < secondOrderBodyList_.size( ); i++ )
        {
            anisotropicPotential += calculateExternalMatrixPotential(
                        currentBodyGravitationalParameters_.at( secondOrderBodyList_.at( i ) ),
                        currentBodyDistances_.at( secondOrderBodyList_.at( i ) ),
                        currentBodyRelativePositions_.at( secondOrderBodyList_.at( i ) ) );
        }
        return anisotropicPotential;
    }

    std::vector< std::string > bodyList_;
    std::vector< std::function< double( ) > > bodyGravitationalParameterFunctions_;
    std::vector< std::function< Eigen::Vector6d( ) > > bodyStateFunctions_;
    std::vector< int > secondOrderBodyList_;
    std::vector< std::function< Eigen::Vector3d( const double ) > > bodyAccelerationFunctions_;
    std::map< int, std::function< double( ) > > bodyAngularMomentumFunctions_;
    std::map< int, std::function< double( const Eigen::Vector3d& ) > > bodySphericalHarmonicPotentialFunctions_;
    std::map< int, std::shared_ptr< SphericalHarmonicWrapper > > bodySphericalHarmonicGravityWrappers_;
    std::map< int, std::function< void( const double ) > > rotationUpdateFunctions_;
    std::shared_ptr< PPNParameterSet > ppnParameterSet_;

    std::vector< double > currentBodyGravitationalParameters_;
    std::vector< Eigen::Vector6d > currentBodyStates_;
    std::vector< Eigen::Vector3d > currentBodyRelativePositions_;
    std::vector< double > currentBodyDistances_;
    std::map< int, double > currentSphericalHarmonicPotentials_;
    std::vector< double > currentExternalScalarPotentials_;
    std::vector< Eigen::Vector3d > currentScalarPotentialGradients_;
    std::vector< Eigen::Vector3d > currentBodyAccelerations_;

    double currentTotalExternalScalarPotential_;
    double currentTotalExternalScalarPotentialSquaredTerms_;
    double currentSecondOrderPotentialSquaredTerms_;
    Eigen::Vector3d currentVectorPotential_;
    Eigen::Matrix3d currentAnisotropicMatrixPotential_;
    double currentScalarPotentialTimePartial_;
    double currentInherentScalarPotentialTimePartial_;
    Eigen::Vector3d currentScalarPotentialPositionPartial_;
    Eigen::Vector3d currentVectorPotentialTimePartial_;
    Eigen::Vector3d currentInherentVectorPotentialTimePartial_;
    Eigen::Matrix3d currentVectorPotentialPositionPartial_;
};

} // namespace relativity

} // namespace tudat

#endif // TUDAT_SOLAR_SYSTEM_METRIC_H
