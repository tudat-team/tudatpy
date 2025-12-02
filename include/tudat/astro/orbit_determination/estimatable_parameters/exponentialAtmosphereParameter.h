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
    ExponentialAtmosphereParameter( std::shared_ptr< aerodynamics::ExponentialAtmosphere > & exponentialAtmosphereModel,
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
                break;
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
                break;
        }
    }

    //! Function to retrieve the size of the rtg force vector parameter (always 3)

    int getParameterSize( )
    {
        return parameterSize_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription
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
    std::shared_ptr< aerodynamics::ExponentialAtmosphere > & exponentialAtmosphereModel_;

    //! Number of rtg acceleration components that are to be estimated.
    int parameterSize_;

    //! List of component indices in rtg accelerations that are to be estimated.
    Eigen::Vector3i accelerationIndices_;
};


//
////! Interface class for estimation the scale height of a body's exponential atmosphere model
///*!
// * Interface class for estimation the scale height of a body's exponential atmosphere model
// */
//class ExponentialAtmosphereScaleHeight : public EstimatableParameter< double >
//{
//public:
//    //! Constructor
//    /*!
//     *  Constructor
//     *  \param exponentialAtmosphereModel Class defining properties of exponentialAtmosphereModel used for body.
//     *  \param associatedBody Body associated with atmosphere mdoel.
//     */
//
//    ExponentialAtmosphereScaleHeight( std::shared_ptr< aerodynamics::ExponentialAtmosphere > & exponentialAtmosphereModel,
//                                      const std::string& associatedBody ):
//        EstimatableParameter< double >( exponential_atmosphere_base_density, associatedBody ), exponentialAtmosphereModel_( exponentialAtmosphereModel ),
//        parameterSize_( 1 )
//    {}
//
//    //! Destructor
//    ~ExponentialAtmosphereScaleHeight( ) {}
//
//    //! Get value of atmosphere scale height
//    /*!
//     *  Get value of atmosphere scale height
//     *  \return Value of atmosphere scale height
//     */
//    double getParameterValue( )
//    {
//        double parameterValue = exponentialAtmosphereModel_->getScaleHeight( );
//        return parameterValue;
//    }
//
//    //! Reset value of atmosphere scale height
//    /*!
//     *  Reset value of atmosphere scale height
//     *  \param parameterValue New value of atmosphere scale height
//     */
//    void setParameterValue( const double parameterValue )
//    {
//        exponentialAtmosphereModel_->resetScaleHeight( parameterValue );
//    }
//
//    //! Function to retrieve the size of the rtg force vector parameter (always 3)
//
//    int getParameterSize( )
//    {
//        return parameterSize_;
//    }
//
//    std::string getParameterDescription( )
//    {
//        std::string parameterDescription = "Scale height of exponential atmosphere model";
//        return parameterDescription;
//    }
//
//    //! Function to retrieve list of components in rtg accelerations that are to be estimated it.
//    Eigen::Vector3i getIndices( )
//    {
//        return accelerationIndices_;
//    }
//
//protected:
//private:
//    //! Class defining properties of atmosphere model
//    std::shared_ptr< aerodynamics::ExponentialAtmosphere > & exponentialAtmosphereModel_;
//
//    //! Number of rtg acceleration components that are to be estimated.
//    int parameterSize_;
//
//    //! List of component indices in rtg accelerations that are to be estimated.
//    Eigen::Vector3i accelerationIndices_;
//};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_EXPONENTIALATMOSPHEREPARAMETER_H
