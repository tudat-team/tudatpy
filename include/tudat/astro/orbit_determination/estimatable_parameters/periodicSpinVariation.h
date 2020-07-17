/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PERIODICSPINVARIATION_H
#define TUDAT_PERIODICSPINVARIATION_H


#include "tudat/astro/orbit_determination/EstimatableParameters/estimatableParameter.h"
#include "tudat/astro/ephemerides/fullPlanetaryRotationModel.h"
#include "tudat/simulation/environment/body.h"
#include "tudat/simulation/environment/createRotationModel.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's periodic spin variation (for a full planetary rotational model).
/*!
 *  Interface class for estimation of a body's periodic spin variation (for a full planetary rotational model).
 *  Interfaces the estimation with the periodic spin variation (rotational rate corrections) of a PlanetaryRotationModel
 *  object
 */
class PeriodicSpinVariation: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param rotationModel PlanetaryRotationModel object of which the periodic spin variation is a property.
     *  \param associatedBody Name of body of which parameter is a property.
     */
    PeriodicSpinVariation(
            const std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel,
            const std::string& associatedBody):
        EstimatableParameter< Eigen::VectorXd >( periodic_spin_variation, associatedBody ),
        rotationModel_( rotationModel ), maxOrder_( rotationModel->getPlanetaryOrientationAngleCalculator()
                                                    ->getRotationRateCorrections().size() ) { }

    //! Destructor
    ~PeriodicSpinVariation( ) { }


    //! Get value of the periodic spin variation coefficients (starting from the lowest order to the highest one, first cosinus
    //! coefficient directly followed by sinus coefficient for each order) .
    /*!
     *  Get value of pole right ascension and declination (starting from the lowest order to the highest one, first cosinus
    //! coefficient directly followed by sinus coefficient for each order) (eg. [cosinus coefficient lowest order,
    //! sinus coefficient lowest order, cosinus coefficient second to lowest order, ...])
     *  \return Periodic spin variation coefficients (starting from the lowest order to the highest one, first cosinus
    //! coefficient directly followed by sinus coefficient for each order).
     */
    Eigen::VectorXd getParameterValue( )
    {
        std::vector< double > parameterValues;
        std::map< double, std::pair< double, double > > rotationRateCorrections =
                rotationModel_->getPlanetaryOrientationAngleCalculator()->getRotationRateCorrections();

        for( std::map< double, std::pair< double, double > >::iterator correctionsIterator = rotationRateCorrections.begin( );
             correctionsIterator != rotationRateCorrections.end( ); correctionsIterator++ )
        {
            parameterValues.push_back( correctionsIterator->second.first );
            parameterValues.push_back( correctionsIterator->second.second );
        }
        return ( utilities::convertStlVectorToEigenVector( parameterValues ) );
    }


    //! Reset value of the periodic spin variation coefficients (starting from the lowest order to the highest one, first cosinus
    //! coefficient directly followed by sinus coefficient for each order).
    /*!
     *  Reset value of the periodic spin variation coefficients (starting from the lowest order to the highest one, first cosinus
    //! coefficient directly followed by sinus coefficient for each order).
     *  \param parameterValue New values for periodic spin variation coefficients.
     */
    void setParameterValue( const Eigen::VectorXd parameterValue )
    {

        std::map< double, std::pair< double, double > > oldRotationRateCorrections
                = rotationModel_->getPlanetaryOrientationAngleCalculator()->getRotationRateCorrections();

        std::map< double, std::pair< double, double > > rotationRateCorrections;

        int currentCorrectionIndex = 0;
        for ( std::map< double, std::pair< double, double > >::iterator oldCorrectionsIterator = oldRotationRateCorrections.begin() ;
              oldCorrectionsIterator != oldRotationRateCorrections.end() ; oldCorrectionsIterator++ ){

            rotationRateCorrections[ oldCorrectionsIterator->first ]
                    = std::make_pair( parameterValue[ currentCorrectionIndex ], parameterValue[ currentCorrectionIndex + 1 ] );

            currentCorrectionIndex += 2;

        }

        rotationModel_->getPlanetaryOrientationAngleCalculator()->resetRotationRateCorrections(
                    rotationRateCorrections);

    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value, which is 2 times the maximum order of the periodic spin variation coefficients
     *  for this parameter.
     */
    int getParameterSize( )
    {
        return 2 * maxOrder_;
    }

protected:

private:

    //! PlanetaryRotationModel object of which periodic spin variation is a property
    std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel_;
    int maxOrder_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_PERIODICSPINVARIATION_H
