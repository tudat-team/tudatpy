/*    Copyright (c) 2010-2017, Delft University of Technology
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

#include "Tudat/Astrodynamics/Gravitation/mutualSphericalHarmonicGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityField.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/sphericalHarmonicAccelerationPartial.h"


namespace tudat
{

namespace acceleration_partials
{
//! Class for calculating partial derivatives of a spherical harmonic gravitational acceleration.
/*!
 *  Class for calculating partial derivatives of a spherical harmonic gravitational acceleration, as calculated by the
 *  SphericalHarmonicsGravitationalAccelerationModel class.
 */
class MutualSphericalHarmonicsGravityPartial: public AccelerationPartial
{
public:

    MutualSphericalHarmonicsGravityPartial(
            const boost::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialOfShExpansionOfBodyExertingAcceleration,
            const boost::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialOfShExpansionOfBodyUndergoingAcceleration,
            const std::string& acceleratedBody, const std::string& acceleratingBody, const bool accelerationUsesMutualAttraction ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::mutual_spherical_harmonic_gravity ),
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_( accelerationPartialOfShExpansionOfBodyExertingAcceleration ),
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_( accelerationPartialOfShExpansionOfBodyUndergoingAcceleration ),
        accelerationUsesMutualAttraction_( accelerationUsesMutualAttraction ){ }

    MutualSphericalHarmonicsGravityPartial(
            const boost::shared_ptr< MutualSphericalHarmonicsGravityPartial > originalAccelerationPartial ):
        AccelerationPartial( originalAccelerationPartial->getAcceleratedBody( ), originalAccelerationPartial->getAcceleratingBody( ),
                             basic_astrodynamics::mutual_spherical_harmonic_gravity ),
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_(
            originalAccelerationPartial->getAccelerationPartialOfShExpansionOfBodyExertingAcceleration( ) ),
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_(
            originalAccelerationPartial->getAccelerationPartialOfShExpansionOfBodyUndergoingAcceleration( ) ),
        accelerationUsesMutualAttraction_( originalAccelerationPartial->getAccelerationUsesMutualAttraction( ) ){ }

    virtual ~MutualSphericalHarmonicsGravityPartial( ) { }

    void wrtPositionOfAcceleratedBody(
                Eigen::Block< Eigen::MatrixXd > partialMatrix,
                const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->wrtPositionOfAcceleratedBody(
                    partialMatrix, addContribution, startRow, startColumn );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->wrtPositionOfAcceleratingBody(
                    partialMatrix, !addContribution, startRow, startColumn );
    }

    void wrtPositionOfAcceleratingBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->wrtPositionOfAcceleratingBody(
                    partialMatrix, addContribution, startRow, startColumn );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->wrtPositionOfAcceleratedBody(
                    partialMatrix, !addContribution, startRow, startColumn );
    }


    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter );

    std::pair< boost::function<void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    void wrtGravitationalParameterOfCentralBody( Eigen::MatrixXd& partialMatrix )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->wrtGravitationalParameterOfCentralBody(
                    partialMatrix, 0 );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->wrtGravitationalParameterOfCentralBody(
                    partialMatrix, -1 );
    }


    virtual void update( const double currentTime );


    int setParameterPartialUpdateFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        int partialSize = 0;

        if( parameter->getParameterName( ).first == estimatable_parameters::gravitational_parameter )
        {
            std::pair< boost::function< void( Eigen::MatrixXd& ) >, int >  partialFunction  =
                    getParameterPartialFunction( parameter );
            partialSize = partialFunction.second;

            if( partialFunction.second > 0 )
            {
                parameterDoublePartialFunctions_[ parameter ] = partialFunction.first;
                isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
                currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, 1 );
            }
        }
        else
        {
            std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromExertingExpansion =
                    accelerationPartialOfShExpansionOfBodyExertingAcceleration_->getParameterPartialFunction( parameter );
            if( partialFunctionFromExertingExpansion.second > 0 )
            {
                accelerationPartialOfShExpansionOfBodyExertingAcceleration_->setParameterPartialUpdateFunction( parameter );
            }

            std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromUndergoingGravity =
                    accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
            if( partialFunctionFromUndergoingGravity.second > 0 )
            {
                accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
            }

            if( partialFunctionFromUndergoingGravity.second > 0 || partialFunctionFromExertingExpansion.second > 0 )
            {
                parameterDoublePartialFunctions_[ parameter ] =
                        getCombinedCurrentDoubleParameterFunction(
                            accelerationPartialOfShExpansionOfBodyExertingAcceleration_,
                            accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_,
                            parameter, partialFunctionFromExertingExpansion.second, partialFunctionFromUndergoingGravity.second, 1 );
                isCurrentDoubleParameterPartialSet_[ parameter ] = 0;
                currentDoubleParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, 1 );
            }
            partialSize = std::max( partialFunctionFromExertingExpansion.second, partialFunctionFromUndergoingGravity.second );
        }
        return partialSize;
    }


    int setParameterPartialUpdateFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromExertingExpansion =
                accelerationPartialOfShExpansionOfBodyExertingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromExertingExpansion.second > 0 )
        {
            accelerationPartialOfShExpansionOfBodyExertingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > partialFunctionFromUndergoingGravity =
                accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->getParameterPartialFunction( parameter );
        if( partialFunctionFromUndergoingGravity.second > 0 )
        {
            accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->setParameterPartialUpdateFunction( parameter );
        }

        if( partialFunctionFromUndergoingGravity.second > 0 || partialFunctionFromExertingExpansion.second > 0 )
        {
            parameterVectorPartialFunctions_[ parameter ] =
                    getCombinedCurrentVectorParameterFunction(
                        accelerationPartialOfShExpansionOfBodyExertingAcceleration_,
                        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_,
                        parameter, partialFunctionFromExertingExpansion.second, partialFunctionFromUndergoingGravity.second, 1 );
            isCurrentVectorParameterPartialSet_[ parameter ] = 0;
            currentVectorParameterPartials_[ parameter ] = Eigen::MatrixXd( accelerationSize_, parameter->getParameterSize( ) );

        }
        return std::max( partialFunctionFromExertingExpansion.second, partialFunctionFromUndergoingGravity.second );
    }

    boost::shared_ptr< SphericalHarmonicsGravityPartial > getAccelerationPartialOfShExpansionOfBodyExertingAcceleration( )
    {
        return accelerationPartialOfShExpansionOfBodyExertingAcceleration_;
    }

    boost::shared_ptr< SphericalHarmonicsGravityPartial > getAccelerationPartialOfShExpansionOfBodyUndergoingAcceleration( )
    {
        return accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_;
    }

    bool getAccelerationUsesMutualAttraction( )
    {
        return accelerationUsesMutualAttraction_;
    }

protected:

    void resetTimeOfMemberObjects( )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->resetTime( currentTime_ );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->resetTime( currentTime_ );
    }

    void updateParameterPartialsOfMemberObjects( )
    {
        accelerationPartialOfShExpansionOfBodyExertingAcceleration_->updateParameterPartials( );
        accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_->updateParameterPartials( );
    }

    boost::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialOfShExpansionOfBodyExertingAcceleration_;

    boost::shared_ptr< SphericalHarmonicsGravityPartial > accelerationPartialOfShExpansionOfBodyUndergoingAcceleration_;

    bool accelerationUsesMutualAttraction_;

};

}

}

#endif // TUDAT_MUTUALSPHERICALHARMONICGRAVITYPARTIAL_H
