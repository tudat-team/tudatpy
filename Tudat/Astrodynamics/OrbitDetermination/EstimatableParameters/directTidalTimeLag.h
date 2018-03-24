/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DIRECTTIDALTIMELAG_H
#define TUDAT_DIRECTTIDALTIMELAG_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Gravitation/directTidalDissipationAcceleration.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation the tidal time lag in a direct tidal acceleration model
/*!
 * Interface class for estimation the tidal time lag in a direct tidal acceleration modell (e.g. without modification of
 * deformed body's gravity field. Modifies the tidal time lag parameter of a (set of) DirectTidalDissipationAcceleration objects
 */
class DirectTidalTimeLag: public EstimatableParameter< double >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param tidalAccelerationModels List of acceleration models of which the tidal time lag is to be estimated
     * \param deformedBody Name of body being tidally deformed
     * \param bodiesCausingDeformation List of bodies causing tidal deformation (empty if all bodies causing deformation in
     * AccelerationMap are used)
     */
    DirectTidalTimeLag(
            const std::vector< boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > tidalAccelerationModels,
            const std::string& deformedBody,
            const std::vector< std::string > bodiesCausingDeformation = std::vector< std::string >( ) ):
        EstimatableParameter< double >(
            direct_dissipation_tidal_time_lag, deformedBody ),
        tidalAccelerationModels_( tidalAccelerationModels ),
        bodiesCausingDeformation_( bodiesCausingDeformation )
    {
        for( unsigned int i = 1; i < tidalAccelerationModels_.size( ); i++ )
        {
            // Check whetehr input is fully consistent at iteration 0.
            if( tidalAccelerationModels_.at( i )->getTimeLag( ) != tidalAccelerationModels_.at( 0 )->getTimeLag( ) )
            {
                std::cerr << "Warning when making direct tidal time lag parameter. Time lags are different in model upon creation, but will be estimated to the same value" << std::endl;
            }
        }
    }

    //! Function to retrieve the tidal time lag value
    /*!
     * Function to retrieve the tidal time lag value
     * \return Current tidal time lag value
     */
    double getParameterValue( )
    {
        return tidalAccelerationModels_.at( 0 )->getTimeLag( );
    }

    //! Function to reset the tidal time lag value
    /*!
     * Function to reset the tidal time lag value
     * \param parameterValue New tidal time lag value
     */
    void setParameterValue( double parameterValue )
    {
        for( unsigned int i = 0; i < tidalAccelerationModels_.size( ); i++ )
        {
            tidalAccelerationModels_.at( i )->resetTimeLag( parameterValue );
        }
    }


    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value, 1 for this parameter
     */
    int getParameterSize( )
    {
        return 1;
    }

    //! Function to retrieve the list of acceleration models of which the tidal time lag is to be estimated
    /*!
     * Function to retrieve list of acceleration models of which the tidal time lag is to be estimated
     * \return List of acceleration models of which the tidal time lag is to be estimated
     */
    std::vector< boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > getTidalAccelerationModels( )
    {
        return tidalAccelerationModels_;
    }

    //! Function to retrieve the list of bodies causing tidal deformation
    /*!
     * Function to retrieve list of bodies causing tidal deformation
     * \return List ofbodies causing tidal deformation
     */
    std::vector< std::string > getBodiesCausingDeformation( )
    {
        return bodiesCausingDeformation_;
    }


private:

    //! List of acceleration models of which the tidal time lag is to be estimated
    std::vector< boost::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > tidalAccelerationModels_;

    //! List of bodies causing tidal deformation (empty if all bodies causing deformation in AccelerationMap are used)
    std::vector< std::string > bodiesCausingDeformation_;

};

}

}

#endif // TUDAT_DIRECTTIDALTIMELAG_H
