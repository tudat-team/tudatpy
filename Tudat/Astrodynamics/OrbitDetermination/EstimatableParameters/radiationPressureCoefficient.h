/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_RADIATIONPRESSURECOEFFICIENT_H
#define TUDAT_RADIATIONPRESSURECOEFFICIENT_H

#include <boost/make_shared.hpp>
#include <boost/bind.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/utilities.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"
#include "Tudat/Mathematics/Interpolators/piecewiseConstantInterpolator.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a radiation pressure coefficient
class RadiationPressureCoefficient: public EstimatableParameter< double >
{

public:
    //! Constructor.
    /*!
     * Constructor
     * \param radiationPressureInterface Object containing the radiation pressure coefficient to be estimated.
     * \param associatedBody Name of body containing the radiationPressureInterface object
     */
    RadiationPressureCoefficient(
            std::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface,
            std::string& associatedBody ):
        EstimatableParameter< double >( radiation_pressure_coefficient, associatedBody ),
        radiationPressureInterface_( radiationPressureInterface )
    { }

    //! Destructor.
    ~RadiationPressureCoefficient( ) { }

    //! Function to get the current value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to get the current value of the radiation pressure coefficient that is to be estimated.
     * \return Current value of the radiation pressure coefficient that is to be estimated.
     */
    double getParameterValue( )
    {
        return radiationPressureInterface_->getRadiationPressureCoefficient( );
    }

    //! Function to reset the value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to reset the value of the radiation pressure coefficient that is to be estimated.
     * \param parameterValue New value of the radiation pressure coefficient that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        radiationPressureInterface_->resetRadiationPressureCoefficient( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( ){ return 1; }

protected:

private:

    //! Object containing the radiation pressure coefficient to be estimated.
    std::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface_;
};

//! Interface class for the estimation of an arc-wise (piecewise constant) radiation pressure coefficient
class ArcWiseRadiationPressureCoefficient: public EstimatableParameter< Eigen::VectorXd >
{

public:
    //! Constructor.
    /*!
     * Constructor
     * \param radiationPressureInterface Object containing the radiation pressure coefficient to be estimated.
     * \param timeLimits Times at which the arcs are to start.
     * \param associatedBody Name of body containing the radiationPressureInterface object
     */
    ArcWiseRadiationPressureCoefficient(
            const std::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface,
            const std::vector< double > timeLimits,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( arc_wise_radiation_pressure_coefficient, associatedBody ),
        radiationPressureInterface_( radiationPressureInterface ), timeLimits_( timeLimits )
    {
        double radiationPressureCoefficient = radiationPressureInterface->getRadiationPressureCoefficient( );
        for( unsigned int i = 0; i < timeLimits.size( ); i++ )
        {
            radiationPressureCoefficients_.push_back( radiationPressureCoefficient );
        }

        timeLimits_.push_back( std::numeric_limits< double >::max( ) );
        fullRadiationPressureCoefficients_ = radiationPressureCoefficients_;
        fullRadiationPressureCoefficients_.push_back( radiationPressureCoefficient );


        coefficientInterpolator_ = std::make_shared< interpolators::PiecewiseConstantInterpolator< double, double > >(
                    timeLimits_, fullRadiationPressureCoefficients_ );

        typedef interpolators::OneDimensionalInterpolator< double, double > LocalInterpolator;
        radiationPressureInterface->resetRadiationPressureCoefficientFunction(
                    std::bind(
                        static_cast< double( LocalInterpolator::* )( const double ) >
                        ( &LocalInterpolator::interpolate ), coefficientInterpolator_, std::placeholders::_1 ) );
    }

    //! Destructor.
    ~ArcWiseRadiationPressureCoefficient( ) { }

    //! Function to get the current value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to get the current value of the radiation pressure coefficient that is to be estimated.
     * \return Current value of the radiation pressure coefficient that is to be estimated.
     */
    Eigen::VectorXd getParameterValue( )
    {
        return utilities::convertStlVectorToEigenVector( radiationPressureCoefficients_ );
    }

    //! Function to reset the value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to reset the value of the radiation pressure coefficient that is to be estimated.
     * \param parameterValue New value of the radiation pressure coefficient that is to be estimated.
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        if( static_cast< int >( radiationPressureCoefficients_.size( ) ) !=
                static_cast< int >( parameterValue.rows( ) ) )
        {
            throw std::runtime_error( "Error when resetting arc-wise radiation pressure coefficients, sizes are incompatible" );
        }

        radiationPressureCoefficients_ = utilities::convertEigenVectorToStlVector( parameterValue );
        for( unsigned int i = 0; i < radiationPressureCoefficients_.size( ); i++ )
        {
            fullRadiationPressureCoefficients_[ i ] = radiationPressureCoefficients_[ i ];
        }
        fullRadiationPressureCoefficients_[ radiationPressureCoefficients_.size( ) ] =
                    radiationPressureCoefficients_.at( radiationPressureCoefficients_.size( ) - 1 );
        coefficientInterpolator_->resetDependentValues( fullRadiationPressureCoefficients_ );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value
     */
    int getParameterSize( ){ return radiationPressureCoefficients_.size( ); }


    std::shared_ptr< interpolators::LookUpScheme< double > > getArcTimeLookupScheme( )
    {
        return coefficientInterpolator_->getLookUpScheme( );
    }

protected:

private:

    //! Object containing the radiation pressure coefficient to be estimated.
    std::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface_;

    //! Times at which the arcs are to start (including end time at maximum double value).
    std::vector< double > timeLimits_;

    //! Values of radiation pressure coefficients in each arc.
    std::vector< double > radiationPressureCoefficients_;

    //! Values of radiation pressure coefficients in each arc, with additional value copied at end.
    std::vector< double > fullRadiationPressureCoefficients_;

    //! Interpolator that returns the radiation pressure coefficient as a function of time.
    std::shared_ptr< interpolators::PiecewiseConstantInterpolator< double, double > > coefficientInterpolator_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_RADIATIONPRESSURECOEFFICIENT_H
