/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_INVERSETIDALQUALITYFACTOR_H
#define TUDAT_INVERSETIDALQUALITYFACTOR_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/gravitation/directTidalDissipationAcceleration.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimating the inverse of the tidal quality factor Q in a direct tidal acceleration model
/*!
* Interface class for estimating the inverse of the tidal quality factor Q in a direct tidal acceleration model (e.g. without modification of
* deformed body's gravity field. Modifies the tidal quality factor of a (set of) DirectTidalDissipationAcceleration objects
*/
class InverseTidalQualityFactor: public EstimatableParameter< double >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param tidalAccelerationModels List of acceleration models of which the tidal quality factor is to be estimated
     * \param deformedBody Name of body being tidally deformed
     * \param bodiesCausingDeformation List of bodies causing tidal deformation (empty if all bodies causing deformation in
     * AccelerationMap are used)
     */
    InverseTidalQualityFactor(
            const std::vector< std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > tidalAccelerationModels,
            const std::string& deformedBody,
            const std::vector< std::string > bodiesCausingDeformation = std::vector< std::string >( ) ):
            EstimatableParameter< double >( inverse_tidal_quality_factor, deformedBody ),
            tidalAccelerationModels_( tidalAccelerationModels ),
            bodiesCausingDeformation_( bodiesCausingDeformation )
    {
        for( unsigned int i = 1; i < tidalAccelerationModels_.size( ); i++ )
        {
            // Check whether input is fully consistent at iteration 0.
            if ( std::isnan( tidalAccelerationModels_.at( i )->getInverseTidalQualityFactor( ) ) || std::isnan( tidalAccelerationModels_.at( i )->getTidalPeriod(  ) ) )
            {
                throw std::runtime_error( "Error when creating inverse tidal quality factor parameter, no value provided for Q and/or tidal period "
                                          " in the acceleration model." );
            }
            if( tidalAccelerationModels_.at( i )->getInverseTidalQualityFactor( ) != tidalAccelerationModels_.at( 0 )->getInverseTidalQualityFactor( ) )
            {
                std::cerr << "Warning when creating inverse tidal quality factor parameter. "
                             "Quality factors are different in model upon creation, but will be estimated to the same value." << std::endl;
            }
        }
    }

    //! Function to retrieve the inverse of the tidal quality factor.
    /*!
     * Function to retrieve the inverse of the tidal quality foator.
     * \return Current value of the inverse of the tidal quality factor.
     */
    double getParameterValue( )
    {
        return tidalAccelerationModels_.at( 0 )->getInverseTidalQualityFactor( );
    }

    //! Function to reset the value for the inverse of the tidal quality factor
    /*!
     * FFunction to reset the value for the inverse of the tidal quality factor
     * \param parameterValue New value for inverse of the tidal quality factor
     */
    void setParameterValue( double parameterValue )
    {
        for( unsigned int i = 0; i < tidalAccelerationModels_.size( ); i++ )
        {
            tidalAccelerationModels_.at( i )->resetInverseTidalQualityFactor( parameterValue );
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
    std::vector< std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > getTidalAccelerationModels( )
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

    std::string getParameterDescription( )
    {
        std::string parameterDescription =
                getParameterTypeString( parameterName_.first ) + "of " + parameterName_.second.first;
        for ( unsigned int i = 0 ; i < bodiesCausingDeformation_.size( ) ; i++ )
        {
            if ( i == 0 )
            {
                parameterDescription += " due to ";
            }
            parameterDescription += bodiesCausingDeformation_[ i ];
            if ( i != bodiesCausingDeformation_.size( ) - 1 )
            {
                parameterDescription += " & ";
            }
        }
        return parameterDescription;
    }


private:

    //! List of acceleration models of which the tidal time lag is to be estimated
    std::vector< std::shared_ptr< gravitation::DirectTidalDissipationAcceleration > > tidalAccelerationModels_;

    //! List of bodies causing tidal deformation (empty if all bodies causing deformation in AccelerationMap are used)
    std::vector< std::string > bodiesCausingDeformation_;

};

}

}

#endif // TUDAT_INVERSETIDALQUALITYFACTOR_H