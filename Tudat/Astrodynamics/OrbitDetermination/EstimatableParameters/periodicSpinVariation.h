/*    Copyright (c) 2010-2018, Delft University of Technology
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


#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/fullPlanetaryRotationModel.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/SimulationSetup/EnvironmentSetup/createRotationModel.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's full rotational model orientation angles (angles psi, I and phi, as well as their time derivatives).
/*!
 *  Interface class for estimation of a body's full rotational model orientation angles (angles psi, I and phi, as well as their time derivatives).
 *  Interfaces the estimation with the orientation angles defining a PlanetaryRotationModel
 *  object
 */
class PeriodicSpinVariation: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param rotationModel PlanetaryRotationModel object of which the periodic spin variation is a property
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

    //! Get value of pole right ascension and declination (in that order)
    /*!
     *  Get value of pole right ascension and declination (in that order)
     *  \return Right ascension and declination (in that order)
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

    //! Reset value of pole right ascension and declination (in that order)
    /*!
     *  Reset value of pole right ascension and declination (in that order)
     *  \param parameterValue New right ascension and declination (in that order)
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
     *  \return Size of parameter value, 2 for this parameter
     */
    int getParameterSize( )
    {
        return 2 * maxOrder_;
    }

protected:

private:

    //! PlanetaryRotationModel object of which rotation rate parameter is a property
    std::shared_ptr< ephemerides::PlanetaryRotationModel > rotationModel_;
    int maxOrder_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_PERIODICSPINVARIATION_H
