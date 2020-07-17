/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_POLARMOTIONAMPLITUDE_H
#define TUDAT_POLARMOTIONAMPLITUDE_H


#include "tudat/astro/orbit_determination/EstimatableParameters/estimatableParameter.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/simulation/environment/body.h"
#include "tudat/simulation/environment/createRotationModel.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's polar motion amplitude (for a full planetary rotational model).
/*!
 *  Interface class for estimation of a body's polar motion amplitude (for a full planetary rotational model).
 *  Interfaces the estimation with the polar motion amplitude of a PlanetaryRotationModel
 *  object
 */
class PolarMotionAmplitude: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param rotationModel PlanetaryRotationModel object of which the polar motion amplitude is a property.
     *  \param associatedBody Name of body of which parameter is a property.
     */
    PolarMotionAmplitude(
            const std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel,
            const std::string& associatedBody):
        EstimatableParameter< Eigen::VectorXd >( polar_motion_amplitude, associatedBody ),
        rotationModel_( rotationModel ),
        maxOrderXpolarMotionCoefficients_( rotationModel->getPlanetaryOrientationAngleCalculator()->getXpolarMotionCoefficients().size() ),
        maxOrderYpolarMotionCoefficients_( rotationModel->getPlanetaryOrientationAngleCalculator()->getYpolarMotionCoefficients().size() )
        { }

    //! Destructor
    ~PolarMotionAmplitude( ) { }


    //! Get value of the amplitudes of the polar motion (starting from the lowest order to the highest one, first cosinus
    //! coefficient directly followed by sinus coefficient for the x- polar motion, and then same order
    //! (cosinus followed by sinus coefficient) for the y-axis polar motion, so 4 amplitude coefficients for each order).
    /*!
     *  Get value of the amplitudes of the polar motion (starting from the lowest order to the highest one, first cosinus
     *  coefficient directly followed by sinus coefficient for the x- polar motion, and then same order
     *  (cosinus followed by sinus coefficient) for the y-axis polar motion, so 4 amplitude coefficients for each order).
     *  \return Polar motion amplitudes.
     */
    Eigen::VectorXd getParameterValue( )
    {
        std::vector< double > parameterValues;

        std::map< double, std::pair< double, double > > xPolarMotionCoefficients = rotationModel_->getPlanetaryOrientationAngleCalculator()
                ->getXpolarMotionCoefficients();
        std::map< double, std::pair< double, double > > yPolarMotionCoefficients = rotationModel_->getPlanetaryOrientationAngleCalculator()
                ->getYpolarMotionCoefficients();

        if ( xPolarMotionCoefficients.size() != yPolarMotionCoefficients.size() ){
            throw std::runtime_error( "Error, unconsistent sizes when comparing x and y polar motion coefficients" );
        }
        else {

            for( std::map< double, std::pair< double, double > >::iterator coefficientsIterator = xPolarMotionCoefficients.begin( );
                 coefficientsIterator != xPolarMotionCoefficients.end( ); coefficientsIterator++ )
            {
                // Retrieve x polar motion coefficients
                parameterValues.push_back( coefficientsIterator->second.first );
                parameterValues.push_back( coefficientsIterator->second.second );

                // Retrieve y polar motion coefficients
                parameterValues.push_back( yPolarMotionCoefficients[ coefficientsIterator->first ].first );
                parameterValues.push_back( yPolarMotionCoefficients[ coefficientsIterator->first ].second );
            }

        }

        return ( utilities::convertStlVectorToEigenVector( parameterValues ) );
    }


    //! Reset value of the polar motion amplitude coefficients (starting from the lowest order to the highest one, first cosinus
    //! coefficient directly followed by sinus coefficient for the x- polar motion, and then same order
    //! (cosinus followed by sinus coefficient) for the y-axis polar motion, so 4 amplitude coefficients for each order).
    /*!
     * Reset value of the polar motion amplitude coefficients (starting from the lowest order to the highest one, first cosinus
     * coefficient directly followed by sinus coefficient for the x- polar motion, and then same order
     * (cosinus followed by sinus coefficient) for the y-axis polar motion, so 4 amplitude coefficients for each order).
     *  \param parameterValue New values for polar motion amplitude coefficients.
     */
    void setParameterValue( const Eigen::VectorXd parameterValue )
    {

        // Retrieve current polar motion coefficients
        std::map< double, std::pair< double, double > > oldXpolarMotionAmplitudeCoefficients
                = rotationModel_->getPlanetaryOrientationAngleCalculator()->getXpolarMotionCoefficients(); 
        std::map< double, std::pair< double, double > > oldYpolarMotionAmplitudeCoefficients
                = rotationModel_->getPlanetaryOrientationAngleCalculator()->getYpolarMotionCoefficients();

        std::map< double, std::pair< double, double > > xPolarMotionAmplitudeCoefficients;
        std::map< double, std::pair< double, double > > yPolarMotionAmplitudeCoefficients;

        int currentCoefficientIndex = 0;

        for ( std::map< double, std::pair< double, double > >::iterator oldCoefficientIterator = oldXpolarMotionAmplitudeCoefficients.begin() ;
              oldCoefficientIterator != oldXpolarMotionAmplitudeCoefficients.end() ; oldCoefficientIterator++ ){

            xPolarMotionAmplitudeCoefficients[ oldCoefficientIterator->first ]
                    = std::make_pair( parameterValue[ currentCoefficientIndex ], parameterValue[ currentCoefficientIndex + 1 ] );

            currentCoefficientIndex += 2;

            yPolarMotionAmplitudeCoefficients[ oldCoefficientIterator->first ]
                    = std::make_pair( parameterValue[ currentCoefficientIndex ], parameterValue[ currentCoefficientIndex + 1 ] );

            currentCoefficientIndex += 2;

        }

        // Reset x polar motion amplitude coefficients.
        rotationModel_->getPlanetaryOrientationAngleCalculator()->resetXpolarMotionCoefficients( xPolarMotionAmplitudeCoefficients);

        // Reset y polar motion amplitude coefficients.
        rotationModel_->getPlanetaryOrientationAngleCalculator()->resetYpolarMotionCoefficients( yPolarMotionAmplitudeCoefficients);

    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value, which is 4 times the maximum order of the polar motion amplitude coefficients here.
     */
    int getParameterSize( )
    {
        return 2 * ( maxOrderXpolarMotionCoefficients_ + maxOrderYpolarMotionCoefficients_ );
    }

protected:

private:

    //! PlanetaryRotationModel object of which the polar motion amplitude is a property
    std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel_;
    int maxOrderXpolarMotionCoefficients_;
    int maxOrderYpolarMotionCoefficients_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_POLARMOTIONAMPLITUDE_H
