/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DISPLACEMENTLOVENUMBER_H
#define TUDAT_DISPLACEMENTLOVENUMBER_H

#include <algorithm>
#include <memory>
#include <vector>

#include "tudat/astro/ground_stations/basicTidalBodyDeformation.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of either the displacement Love number h_n or the Shida number l_n of the basic tidal deformation
//! of a body (see basic_astrodynamics::BasicTidalBodyDeformation), for a set of degrees n. The parameter has one entry per degree, and
//! is applied to all basic tidal deformation models of the body.
class DisplacementLoveNumber : public EstimatableParameter< Eigen::VectorXd >
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param deformationModels Basic tidal deformation models of the body
     *  \param degrees Degree of each estimated number
     *  \param useShidaNumber Whether the Shida number l_n (true) or the Love number h_n (false) is estimated
     *  \param associatedBody Deformed body
     */
    DisplacementLoveNumber( const std::vector< std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > >& deformationModels,
                            const std::vector< int >& degrees,
                            const bool useShidaNumber,
                            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( useShidaNumber ? displacement_shida_number : displacement_love_number, associatedBody ),
        deformationModels_( deformationModels ), degrees_( degrees )
    {
        if( deformationModels_.empty( ) )
        {
            throw std::runtime_error( "Error when creating displacement Love or Shida number parameter of " + associatedBody +
                                      ", no basic tidal deformation model found (IERS 2010 deformation models are not supported)" );
        }
        if( degrees_.empty( ) )
        {
            throw std::runtime_error( "Error when creating displacement Love or Shida number parameter of " + associatedBody +
                                      ", no degrees provided" );
        }
        for( unsigned int i = 0; i < degrees_.size( ); i++ )
        {
            if( std::count( degrees_.begin( ), degrees_.end( ), degrees_.at( i ) ) > 1 )
            {
                throw std::runtime_error( "Error when creating displacement Love or Shida number parameter of " + associatedBody +
                                          ", degree " + std::to_string( degrees_.at( i ) ) + " provided more than once" );
            }
        }
        for( const auto& model : deformationModels_ )
        {
            for( unsigned int i = 0; i < degrees_.size( ); i++ )
            {
                if( model->getDisplacementLoveNumberPerDegree( ).count( degrees_.at( i ) ) == 0 )
                {
                    throw std::runtime_error( "Error when creating displacement Love or Shida number parameter of " + associatedBody +
                                              ", degree " + std::to_string( degrees_.at( i ) ) +
                                              " not in a basic tidal deformation model" );
                }
                if( extractLoveOrShidaNumber( model, degrees_.at( i ) ) !=
                    extractLoveOrShidaNumber( deformationModels_.front( ), degrees_.at( i ) ) )
                {
                    throw std::runtime_error( "Error when creating displacement Love or Shida number parameter of " + associatedBody +
                                              ", inconsistent Love or Shida number in the basic tidal deformation models" );
                }
            }
        }
    }

    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd value( degrees_.size( ) );
        for( unsigned int i = 0; i < degrees_.size( ); i++ )
        {
            value( i ) = extractLoveOrShidaNumber( deformationModels_.front( ), degrees_.at( i ) );
        }
        return value;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        for( const auto& model : deformationModels_ )
        {
            std::map< int, std::pair< double, double > > loveAndShidaNumberPerDegree = model->getDisplacementLoveNumberPerDegree( );
            for( unsigned int i = 0; i < degrees_.size( ); i++ )
            {
                ( getUseShidaNumber( ) ? loveAndShidaNumberPerDegree[ degrees_.at( i ) ].second
                                       : loveAndShidaNumberPerDegree[ degrees_.at( i ) ].first ) = parameterValue( i );
            }
            model->setDisplacementLoveNumberPerDegree( loveAndShidaNumberPerDegree );
        }
    }

    int getParameterSize( )
    {
        return degrees_.size( );
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription =
                getParameterTypeString( parameterName_.first ) + "of " + parameterName_.second.first + " for degrees: ";
        for( unsigned int i = 0; i < degrees_.size( ); i++ )
        {
            parameterDescription += std::to_string( degrees_.at( i ) ) + ( i + 1 < degrees_.size( ) ? ", " : "." );
        }
        return parameterDescription;
    }

    std::vector< std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > > getDeformationModels( ) const
    {
        return deformationModels_;
    }

    std::vector< int > getDegrees( ) const
    {
        return degrees_;
    }

    bool getUseShidaNumber( ) const
    {
        return parameterName_.first == displacement_shida_number;
    }

private:
    //! Love number h_n or Shida number l_n of the given degree of a deformation model
    double extractLoveOrShidaNumber( const std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > model, const int degree ) const
    {
        const std::pair< double, double > loveAndShidaNumber = model->getDisplacementLoveNumberPerDegree( ).at( degree );
        return getUseShidaNumber( ) ? loveAndShidaNumber.second : loveAndShidaNumber.first;
    }

    //! Basic tidal deformation models of the body
    std::vector< std::shared_ptr< basic_astrodynamics::BasicTidalBodyDeformation > > deformationModels_;

    //! Degree of each estimated number
    std::vector< int > degrees_;
};

}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_DISPLACEMENTLOVENUMBER_H
