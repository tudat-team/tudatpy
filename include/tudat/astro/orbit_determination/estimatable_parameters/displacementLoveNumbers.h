/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DISPLACEMENTLOVENUMBERS_H
#define TUDAT_DISPLACEMENTLOVENUMBERS_H

#include <memory>
#include <vector>

#include "tudat/astro/ground_stations/basicTidalBodyDeformation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of either the displacement Love numbers h_n or the Shida numbers l_n of the basic tidal deformation
//! of a body (see basic_astrodynamics::BasicTidalBodyDeformation), for a set of degrees n. The parameter has one entry per degree, and
//! is applied to all basic tidal deformation models of the body.
class DisplacementLoveNumbers : public EstimatableParameter< Eigen::VectorXd >
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param deformationModels Basic tidal deformation models of the body
     *  \param degrees Degrees of the estimated numbers
     *  \param useShidaNumbers Whether the Shida numbers l_n (true) or the Love numbers h_n (false) are estimated
     *  \param associatedBody Deformed body
     */
    DisplacementLoveNumbers( const std::vector< std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > >& deformationModels,
                             const std::vector< int >& degrees,
                             const bool useShidaNumbers,
                             const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( useShidaNumbers ? displacement_shida_number : displacement_love_number, associatedBody ),
        deformationModels_( deformationModels ), degrees_( degrees ), useShidaNumbers_( useShidaNumbers )
    {
        if( deformationModels_.empty( ) )
        {
            throw std::runtime_error( "Error when creating displacement Love or Shida number parameter of " + associatedBody +
                                      ", no basic tidal deformation model found" );
        }
        for( const auto& model : deformationModels_ )
        {
            for( unsigned int i = 0; i < degrees_.size( ); i++ )
            {
                if( model->getDisplacementLoveNumbers( ).count( degrees_.at( i ) ) == 0 )
                {
                    throw std::runtime_error( "Error when creating displacement Love or Shida number parameter of " + associatedBody +
                                              ", degree " + std::to_string( degrees_.at( i ) ) +
                                              " not in a basic tidal deformation model" );
                }
                if( getNumber( model, degrees_.at( i ) ) != getNumber( deformationModels_.front( ), degrees_.at( i ) ) )
                {
                    throw std::runtime_error( "Error when creating displacement Love or Shida number parameter of " + associatedBody +
                                              ", inconsistent numbers in the basic tidal deformation models" );
                }
            }
        }
    }

    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd value( degrees_.size( ) );
        for( unsigned int i = 0; i < degrees_.size( ); i++ )
        {
            value( i ) = getNumber( deformationModels_.front( ), degrees_.at( i ) );
        }
        return value;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        for( const auto& model : deformationModels_ )
        {
            std::map< int, std::pair< double, double > > loveNumbers = model->getDisplacementLoveNumbers( );
            for( unsigned int i = 0; i < degrees_.size( ); i++ )
            {
                ( useShidaNumbers_ ? loveNumbers[ degrees_.at( i ) ].second : loveNumbers[ degrees_.at( i ) ].first ) = parameterValue( i );
            }
            model->setDisplacementLoveNumbers( loveNumbers );
        }
    }

    int getParameterSize( )
    {
        return degrees_.size( );
    }

    std::vector< std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > > getDeformationModels( ) const
    {
        return deformationModels_;
    }

    std::vector< int > getDegrees( ) const
    {
        return degrees_;
    }

    bool getUseShidaNumbers( ) const
    {
        return useShidaNumbers_;
    }

private:
    //! Love number h_n or Shida number l_n of the given degree of a deformation model
    double getNumber( const std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > model, const int degree ) const
    {
        const std::pair< double, double > loveNumbers = model->getDisplacementLoveNumbers( ).at( degree );
        return useShidaNumbers_ ? loveNumbers.second : loveNumbers.first;
    }

    //! Basic tidal deformation models of the body
    std::vector< std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > > deformationModels_;

    //! Degrees of the estimated numbers
    std::vector< int > degrees_;

    //! Whether the Shida numbers l_n (true) or the Love numbers h_n (false) are estimated
    bool useShidaNumbers_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_DISPLACEMENTLOVENUMBERS_H
