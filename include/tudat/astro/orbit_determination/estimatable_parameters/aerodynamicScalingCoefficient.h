/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_AERODYNAMICSCALINGCOEFFICIENT_H
#define TUDAT_AERODYNAMICSCALINGCOEFFICIENT_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

class AerodynamicScalingFactor : public EstimatableParameter< double >
{
public:
    AerodynamicScalingFactor( const std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > > aerodynamicAccelerations,
                              const EstimatebleParametersEnum parameterType,
                              const std::string& associatedBody ):
        EstimatableParameter< double >( parameterType, associatedBody ), aerodynamicAccelerations_( aerodynamicAccelerations )
    {
        if( ( parameterType != drag_component_scaling_factor ) && ( parameterType != side_component_scaling_factor ) &&
            ( parameterType != lift_component_scaling_factor ) )
        {
            throw std::runtime_error( "Error when creating aerodynamic scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterType ) );
        }

        switch( parameterType )
        {
            case drag_component_scaling_factor:
                parameterIndex_ = 0;
                break;
            case side_component_scaling_factor:
                parameterIndex_ = 1;
                break;
            case lift_component_scaling_factor:
                parameterIndex_ = 2;
                break;
            default:
                throw std::runtime_error( "Error when creating aerodynamic scaling parameter, type is inconsistent: " +
                                          std::to_string( parameterType ) );
        }
    }

    AerodynamicScalingFactor( const std::shared_ptr< aerodynamics::AerodynamicAcceleration > aerodynamicAcceleration,
                              const EstimatebleParametersEnum parameterType,
                              const std::string& associatedBody ):
        AerodynamicScalingFactor( std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > >( { aerodynamicAcceleration } ),
                                  parameterType,
                                  associatedBody )
    {}

    ~AerodynamicScalingFactor( ) {}

    double getParameterValue( )
    {
        double parameterValue = aerodynamicAccelerations_.at( 0 )->getComponentScaling( parameterIndex_ );
        for( unsigned int i = 1; i < aerodynamicAccelerations_.size( ); i++ )
        {
            if( aerodynamicAccelerations_.at( i )->getComponentScaling( parameterIndex_ ) != parameterValue )
            {
                std::cerr << "Warning when retrieving aerodynamic acceleration scaling factor from list. List entries are not equal, "
                             "returning first entry"
                          << std::endl;
                break;
            }
        }
        return parameterValue;
    }

    void setParameterValue( double parameterValue )
    {
        for( unsigned int i = 0; i < aerodynamicAccelerations_.size( ); i++ )
        {
            aerodynamicAccelerations_.at( i )->setComponentScaling( parameterValue, parameterIndex_ );
        }
    }

    int getParameterSize( )
    {
        return 1;
    }

protected:
    std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > > aerodynamicAccelerations_;

    int parameterIndex_;
};

class ArcWiseAerodynamicScalingFactor : public EstimatableParameter< Eigen::VectorXd >
{
public:
    ArcWiseAerodynamicScalingFactor( const std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > > aerodynamicAccelerations,
                                     const EstimatebleParametersEnum parameterType,
                                     const std::vector< double > timeLimits,
                                     const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( parameterType, associatedBody ), aerodynamicAccelerations_( aerodynamicAccelerations ),
        timeLimits_( timeLimits )
    {
        if( ( parameterType != arc_wise_drag_component_scaling_factor ) && ( parameterType != arc_wise_side_component_scaling_factor ) &&
            ( parameterType != arc_wise_lift_component_scaling_factor ) )
        {
            throw std::runtime_error( "Error when creating aerodynamic scaling parameter, type is inconsistent: " +
                                      std::to_string( parameterType ) );
        }

        switch( parameterType )
        {
            case arc_wise_drag_component_scaling_factor:
                parameterIndex_ = 0;
                break;
            case arc_wise_side_component_scaling_factor:
                parameterIndex_ = 1;
                break;
            case arc_wise_lift_component_scaling_factor:
                parameterIndex_ = 2;
                break;
            default:
                throw std::runtime_error( "Error when creating aerodynamic scaling parameter, type is inconsistent: " +
                                          std::to_string( parameterType ) );
        }

        // setup of arcwise parameter object

        for( unsigned int i = 0; i < aerodynamicAccelerations_.size( ); i++ )
        {
            if( std::isnan( aerodynamicAccelerations_.at( i )->getComponentScaling( parameterIndex_ ) ) )
            {
                throw std::runtime_error( "Error when creating estimated arcwise Cr coefficient for " + associatedBody +
                                          ", current Cr not initialized" );
            }

            if( aerodynamicAccelerations_.at( i )->getComponentScalingFunction( parameterIndex_ ) != nullptr )
            {
                throw std::runtime_error( "Error when creating estimated arcwise Cr coefficient for " + associatedBody +
                                          ", time-variable Cr function defined" );
            }
        }

        double originalComponentScalingCoefficient = getParameterValueFromAccelerationModel( );
        for( unsigned int i = 0; i < timeLimits.size( ); i++ )
        {
            componentScalingCoefficients_.push_back( originalComponentScalingCoefficient );
        }

        timeLimits_.push_back( std::numeric_limits< double >::max( ) );
        fullComponentScalingCoefficients_ = componentScalingCoefficients_;
        fullComponentScalingCoefficients_.push_back( originalComponentScalingCoefficient );

        coefficientInterpolator_ = std::make_shared< interpolators::PiecewiseConstantInterpolator< double, double > >(
                timeLimits_, fullComponentScalingCoefficients_ );

        typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;

        for( unsigned int i = 0; i < aerodynamicAccelerations_.size( ); i++ )
        {
            aerodynamicAccelerations_.at( i )->setComponentScalingFunction(
                    std::bind( static_cast< double ( LocalInterpolator::* )( const double ) >( &LocalInterpolator::interpolate ),
                               coefficientInterpolator_,
                               std::placeholders::_1 ),
                    parameterIndex_ );
        }
    }

    ArcWiseAerodynamicScalingFactor( const std::shared_ptr< aerodynamics::AerodynamicAcceleration > aerodynamicAcceleration,
                                     const EstimatebleParametersEnum parameterType,
                                     const std::vector< double > timeLimits,
                                     const std::string& associatedBody ):
        ArcWiseAerodynamicScalingFactor(
                std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > >( { aerodynamicAcceleration } ),
                parameterType,
                timeLimits,
                associatedBody )
    {}

    ~ArcWiseAerodynamicScalingFactor( ) {}

    //! Function to ask for values from acceleration model (only for initialization)
    double getParameterValueFromAccelerationModel( )
    {
        double parameterValue = aerodynamicAccelerations_.at( 0 )->getComponentScaling( parameterIndex_ );
        for( unsigned int i = 1; i < aerodynamicAccelerations_.size( ); i++ )
        {
            if( aerodynamicAccelerations_.at( i )->getComponentScaling( parameterIndex_ ) != parameterValue )
            {
                std::cerr << "Warning when retrieving aerodynamic acceleration scaling factor from list. List entries are not equal, "
                             "returning first entry"
                          << std::endl;
                break;
            }
        }
        return parameterValue;
    }

    //! Function to reset the arc-wise values of the component scaling factors to be estimated.
    /*!
     * Function to reset the arc-wise values of the component scaling factors to be estimated. Note that this function is not binding a new
     * piecewise interpolator function, but changes the data in the existsing (member) interpolator. \param parameterValue New value of the
     * radiation pressure coefficient that is to be estimated.
     */

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        if( static_cast< int >( componentScalingCoefficients_.size( ) ) != static_cast< int >( parameterValue.rows( ) ) )
        {
            throw std::runtime_error( "Error when resetting arc-wise component scaling factors, sizes are incompatible" );
        }

        componentScalingCoefficients_ = utilities::convertEigenVectorToStlVector( parameterValue );
        for( unsigned int i = 0; i < componentScalingCoefficients_.size( ); i++ )
        {
            fullComponentScalingCoefficients_[ i ] = componentScalingCoefficients_[ i ];
        }
        fullComponentScalingCoefficients_[ componentScalingCoefficients_.size( ) ] =
                componentScalingCoefficients_.at( componentScalingCoefficients_.size( ) - 1 );
        coefficientInterpolator_->resetDependentValues( fullComponentScalingCoefficients_ );
    }

    //! Function to get the current value of the arcwise component scaling coefficient.
    /*!
     * Function to get the current value of the arcwise component scaling coefficient that is to be estimated.
     * \return Current value of the arcwise component scaling coefficient that is to be estimated.
     */

    Eigen::VectorXd getParameterValue( )
    {
        return utilities::convertStlVectorToEigenVector( componentScalingCoefficients_ );
    }

    int getParameterSize( )
    {
        return componentScalingCoefficients_.size( );
    }

    int getNumberOfArcs( )
    {
        return getParameterSize( );
    }

    int getParameterIndex( )
    {
        return parameterIndex_;
    }

    std::shared_ptr< interpolators::LookUpScheme< double > > getArcTimeLookupScheme( )
    {
        return coefficientInterpolator_->getLookUpScheme( );
    }

protected:
    std::vector< std::shared_ptr< aerodynamics::AerodynamicAcceleration > > aerodynamicAccelerations_;

    int parameterIndex_;

    //! Times at which the arcs are to start (including end time at maximum double value).
    std::vector< double > timeLimits_;

    //! Values of aerodynamic component scaling coefficients in each arc.
    std::vector< double > componentScalingCoefficients_;

    //! Values of aerodynamic component scaling coefficients in each arc, with additional value copied at end.
    std::vector< double > fullComponentScalingCoefficients_;

    //! Interpolator that returns the aerodynamic component scaling coefficient as a function of time.
    std::shared_ptr< interpolators::PiecewiseConstantInterpolator< double, double > > coefficientInterpolator_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_AERODYNAMICSCALINGCOEFFICIENT_H