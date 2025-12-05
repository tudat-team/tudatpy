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

#ifndef TUDAT_EXPONENTIALATMOSPHEREPARAMETER_H
#define TUDAT_EXPONENTIALATMOSPHEREPARAMETER_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/basic_astro/empiricalAcceleration.h"
#include "tudat/math/interpolators/piecewiseConstantInterpolator.h"
#include "tudat/astro/aerodynamics/exponentialAtmosphere.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation the base density of a body's exponential atmosphere model
/*!
 * Interface class for estimation the base density of a body's exponential atmosphere model
 */
class ExponentialAtmosphereParameter : public EstimatableParameter< double >
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param exponentialAtmosphereModel Class defining properties of exponentialAtmosphereModel used for body.
     *  \param associatedBody Body associated with atmosphere mdoel.
     */
    ExponentialAtmosphereParameter( const std::shared_ptr< aerodynamics::ExponentialAtmosphere >& exponentialAtmosphereModel,
                                    const estimatable_parameters::EstimatebleParametersEnum parameterType,
                                    const std::string& associatedBody ):
        EstimatableParameter< double >( parameterType, associatedBody ), exponentialAtmosphereModel_( exponentialAtmosphereModel ),
        parameterSize_( 1 )
    {
        if( ( parameterType != exponential_atmosphere_base_density ) && ( parameterType != exponential_atmosphere_scale_height ) )
        {
            throw std::runtime_error( "Error when creating exponential atmosphere parameter, type is inconsistent: " +
                                      std::to_string( parameterType ) );
        }

    }

    //! Destructor
    ~ExponentialAtmosphereParameter( ) {}

    //! Get value of atmosphere parameter
    /*!
     *  Get value of atmosphere parameter
     *  \return Value of atmosphere parameter
     */
    double getParameterValue( )
    {
        double parameterValue;
        switch( parameterName_.first )
        {
            case exponential_atmosphere_base_density: {
                parameterValue = exponentialAtmosphereModel_->getBaseDensity( );
                break;
            }
            case exponential_atmosphere_scale_height: {
                parameterValue = exponentialAtmosphereModel_->getScaleHeight( );
                break;
            }
            default:
                throw std::runtime_error( "Error when creating exponential atmosphere parameter, type is inconsistent: " +
                          std::to_string( parameterName_.first ) );
        }
        return parameterValue;
    }

    //! Reset value of atmosphere parameter
    /*!
     *  Reset value of atmosphere parameter
     *  \param parameterValue New value of atmosphere parameter
     */
    void setParameterValue( double newParameterValue)
    {
        switch( parameterName_.first )
        {
            case exponential_atmosphere_base_density: {
                exponentialAtmosphereModel_->resetBaseDensity( newParameterValue );
                break;
            }
            case exponential_atmosphere_scale_height: {
                exponentialAtmosphereModel_->resetScaleHeight( newParameterValue );
                break;
            }
            default:
                throw std::runtime_error( "Error when creating exponential atmosphere parameter, type is inconsistent: " +
                    std::to_string( parameterName_.first ) );
        }
    }

    //! Function to retrieve the size of the rtg force vector parameter (always 3)

    int getParameterSize( )
    {
        return parameterSize_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription;
        switch( parameterName_.first )
        {
            case exponential_atmosphere_base_density: {
                parameterDescription = "Base density of exponential atmosphere model";
                break;
            }
            case exponential_atmosphere_scale_height: {
                parameterDescription = "Scale height of exponential atmosphere model";
                break;
            }
            default:
                break;
        }
        return parameterDescription;
    }

    //! Function to retrieve list of components in rtg accelerations that are to be estimated it.
    Eigen::Vector3i getIndices( )
    {
        return accelerationIndices_;
    }

protected:
private:
    //! Class defining properties of atmosphere model
    const std::shared_ptr< aerodynamics::ExponentialAtmosphere > exponentialAtmosphereModel_;

    //! Number of rtg acceleration components that are to be estimated.
    int parameterSize_;

    //! List of component indices in rtg accelerations that are to be estimated.
    Eigen::Vector3i accelerationIndices_;
};



//! Interface class for arc-wise estimation the base density of a body's exponential atmosphere model
/*!
 * Interface class for arc-wise estimation the base density of a body's exponential atmosphere model
 */
class ArcWiseExponentialAtmosphereParameter : public EstimatableParameter< Eigen::VectorXd >
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param exponentialAtmosphereModel Class defining properties of exponentialAtmosphereModel used for body.
     *  \param associatedBody Body associated with atmosphere mdoel.
     */
    ArcWiseExponentialAtmosphereParameter( const std::shared_ptr< aerodynamics::ExponentialAtmosphere >& exponentialAtmosphereModel,
                                    const estimatable_parameters::EstimatebleParametersEnum parameterType,
                                    const std::vector< double > timeLimits,
                                    const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( parameterType, associatedBody ), exponentialAtmosphereModel_( exponentialAtmosphereModel ),
        timeLimits_( timeLimits )
    {

        // setup of arcwise parameter object
        switch( parameterType )
        {
            case arc_wise_exponential_atmosphere_base_density:

                if( std::isnan( exponentialAtmosphereModel_->getBaseDensity( ) ) )
                {
                    throw std::runtime_error( "Error when creating arc-wise estimated atmospheric base density parameter for " + associatedBody +
                                              ", current base density parameter not initialized" );
                }
                if( exponentialAtmosphereModel_->getBaseDensityFunction( ) != nullptr )
                {
                    throw std::runtime_error( "Error when creating arc-wise estimated atmospheric base density parameter for " + associatedBody +
                                              ", time-varying scaling coefficient function already defined" );
                }

                originalParameterValue_ = exponentialAtmosphereModel_->getBaseDensity( );
                for( unsigned int i = 0; i < timeLimits.size( ); i++ )
                {
                    parameterValues_.push_back( originalParameterValue_ );
                }

                break;

            case arc_wise_exponential_atmosphere_scale_height:

                if( std::isnan( exponentialAtmosphereModel_->getScaleHeight( ) ) )
                {
                    throw std::runtime_error( "Error when creating arc-wise estimated atmospheric scale height parameter for " + associatedBody +
                                              ", current base density parameter not initialized" );
                }
                if( exponentialAtmosphereModel_->getScaleHeightFunction( ) != nullptr )
                {
                    throw std::runtime_error( "Error when creating arc-wise estimated atmospheric scale height parameter for " + associatedBody +
                                              ", time-varying scaling coefficient function already defined" );
                }

                originalParameterValue_ = exponentialAtmosphereModel_->getScaleHeight( );
                for( unsigned int i = 0; i < timeLimits.size( ); i++ )
                {
                    parameterValues_.push_back( originalParameterValue_ );
                }

                break;



            default:
                throw std::runtime_error( "Error when creating arc-wise exponential atmosphere parameter, type is inconsistent: " +
                          std::to_string( parameterType ) );
        }


        timeLimits_.push_back( std::numeric_limits< double >::max( ) );
        fullParameterValues_ = parameterValues_;
        fullParameterValues_.push_back( originalParameterValue_ );

        parameterInterpolator_ = std::make_shared< interpolators::PiecewiseConstantInterpolator< double, double > >(
                timeLimits_, fullParameterValues_ );

        typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;

        exponentialAtmosphereModel_->setParameterFunction(
                std::bind( static_cast< double ( LocalInterpolator::* )( const double ) >( &LocalInterpolator::interpolate ),
                           parameterInterpolator_,
                           std::placeholders::_1 ),
                parameterType );

    }

    //! Destructor
    ~ArcWiseExponentialAtmosphereParameter( ) {}



    //! Get arc-wise values of atmosphere parameter
    /*!
     *  Get values of atmosphere parameter
     *  \return Arc-wise values of atmosphere parameter
     */
    Eigen::VectorXd getParameterValue( )
    {
        return utilities::convertStlVectorToEigenVector( parameterValues_ );
    }

    int getParameterSize( )
    {
        return parameterValues_.size( );
    }

    int getNumberOfArcs( )
    {
        return getParameterSize( );
    }

    //! Reset value of atmosphere parameter
    /*!
     *  Reset value of atmosphere parameter
     *  \param parameterValue New value of atmosphere parameter
     */
    void setParameterValue( double newParameterValue, const estimatable_parameters::EstimatebleParametersEnum parameterType)
    {
        switch( parameterName_.first )
        {
            case exponential_atmosphere_base_density: {
                exponentialAtmosphereModel_->resetBaseDensity( newParameterValue );
                break;
            }
            case exponential_atmosphere_scale_height: {
                exponentialAtmosphereModel_->resetScaleHeight( newParameterValue );
                break;
            }
            default:
                throw std::runtime_error( "Error when creating exponential atmosphere parameter, type is inconsistent: " +
                    std::to_string( parameterName_.first ) );
        }
    }

    void setParameterValue( Eigen::VectorXd newParameterValues )
    {
        if( static_cast< int >( parameterValues_.size( ) ) != static_cast< int >( newParameterValues.rows( ) ) )
        {
            throw std::runtime_error( "Error when resetting values of arc-wise atmosphere parameters, sizes are incompatible" );
        }

        parameterValues_ = utilities::convertEigenVectorToStlVector( newParameterValues );
        for( unsigned int i = 0; i < parameterValues_.size( ); i++ )
        {
            fullParameterValues_[ i ] = parameterValues_[ i ];
        }
        fullParameterValues_[ parameterValues_.size( ) ] =
                parameterValues_.at( parameterValues_.size( ) - 1 );
        parameterInterpolator_->resetDependentValues( fullParameterValues_ );
    }


    std::string getParameterDescription( )
    {
        std::string parameterDescription;
        switch( parameterName_.first )
        {
            case exponential_atmosphere_base_density: {
                parameterDescription = "Base density of exponential atmosphere model";
                break;
            }
            case exponential_atmosphere_scale_height: {
                parameterDescription = "Scale height of exponential atmosphere model";
                break;
            }
            default:
                break;
        }
        return parameterDescription;
    }

    //! Function to retrieve list of components in rtg accelerations that are to be estimated it.
    Eigen::Vector3i getIndices( )
    {
        return accelerationIndices_;
    }

    std::shared_ptr< interpolators::LookUpScheme< double > > getArcTimeLookupScheme( )
    {
        return parameterInterpolator_->getLookUpScheme( );
    }

protected:
private:
    //! Class defining properties of atmosphere model
    const std::shared_ptr< aerodynamics::ExponentialAtmosphere > exponentialAtmosphereModel_;

    //! Times at which the arcs are to start (including end time at maximum double value).
    std::vector< double > timeLimits_;

    //! List of parameter values
    std::vector< double > parameterValues_;
    std::vector< double > fullParameterValues_;

    //! Original <double> parameter value;
    double originalParameterValue_;

    //! Number of rtg acceleration components that are to be estimated.
    int parameterSize_;

    //! List of component indices in rtg accelerations that are to be estimated.
    Eigen::Vector3i accelerationIndices_;

    //! Interpolator that returns the atmosphere parameter as a function of time.
    std::shared_ptr< interpolators::PiecewiseConstantInterpolator< double, double > > parameterInterpolator_;
};


}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_EXPONENTIALATMOSPHEREPARAMETER_H
