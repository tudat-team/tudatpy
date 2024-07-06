/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GRAVITYFIELDVARIATIONPARAMETERS_H
#define TUDAT_GRAVITYFIELDVARIATIONPARAMETERS_H

#include <iostream>
#include <vector>
#include <string>
#include <vector>
#include <map>


#include <type_traits>
#include <memory>
#include <Eigen/Geometry>

#include "tudat/astro/gravitation/polynomialGravityFieldVariations.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/gravitation/periodicGravityFieldVariations.h"

namespace tudat
{

namespace estimatable_parameters
{

class PolynomialGravityFieldVariationsParameters: public EstimatableParameter< Eigen::VectorXd >
{

public:
    //! Constructor.
    /*!
     *  Constructor taking parameter name and associated body. All parameters are identified by a these two variables.
     *  Any additional information that may be required for uniquely defining a parameter is to be defined in the derived class.
     *  \param parameterName Enum value defining the type of the parameter.
     *  \param associatedBody Name of body associated with patameters
     *  \param pointOnBodyId Reference point on body associated with parameter (empty by default).
     */
    PolynomialGravityFieldVariationsParameters(
        const std::shared_ptr< gravitation::PolynomialGravityFieldVariations > polynomialVariationModel,
        const std::map< int, std::vector< std::pair< int, int > > >& cosineBlockIndicesPerPower,
        const std::map< int, std::vector< std::pair< int, int > > >& sineBlockIndicesPerPower,
        const std::string& bodyName ):
        EstimatableParameter< Eigen::VectorXd >( polynomial_gravity_field_variation_amplitudes, bodyName ),
        polynomialVariationModel_( polynomialVariationModel )
        {
            std::map< int, Eigen::MatrixXd > cosineVariations = polynomialVariationModel->getCosineAmplitudes( );
            int cosineIndexCounter = 0;
            for( auto it : cosineBlockIndicesPerPower )
            {
                if( cosineVariations.count( it.first ) == 0 )
                {
                    throw std::runtime_error( "Error when estimating gravity field polynomial corrections of body " + bodyName +
                    ", not polynomial term of order " + std::to_string( it.first ) + " found for cosine coefficients" );
                }

                for( unsigned int i = 0; i < it.second.size( ); i++ )
                {
                    if( it.second.at( i ).first > polynomialVariationModel_->getMaximumDegree( ) ||
                        it.second.at( i ).first < polynomialVariationModel_->getMinimumDegree( ) )
                    {
                        throw std::runtime_error( "Error when estimating gravity field polynomial corrections of body " + bodyName +
                                                  " of polynomial term of order " + std::to_string( it.first ) +
                                                  ", no cosine coefficient variation of degree " + std::to_string( it.second.at( i ).first ) + " found." );
                    }

                    if( it.second.at( i ).second > polynomialVariationModel_->getMaximumOrder( ) ||
                        it.second.at( i ).second < polynomialVariationModel_->getMinimumOrder( ) )
                    {
                        throw std::runtime_error( "Error when estimating gravity field polynomial corrections of body " + bodyName +
                                                  " of polynomial term of order " + std::to_string( it.first ) +
                                                  ", no cosine coefficient variation of order " + std::to_string( it.second.at( i ).second ) + " found." );
                    }
                    cosineCorrectionIndices_.push_back( std::make_tuple( it.first, it.second.at( i ).first, it.second.at( i ).second ) );
                    indexAndPowerPerCosineBlockIndex_[ std::make_pair( it.second.at( i ).first, it.second.at( i ).second ) ].push_back(
                        std::make_pair( cosineIndexCounter, it.first ) );
                    cosineIndexCounter++;
                }
            }

            std::map< int, Eigen::MatrixXd > sineVariations = polynomialVariationModel->getSineAmplitudes( );
            int sineIndexCounter = 0;
            for( auto it : sineBlockIndicesPerPower )
            {
                if( sineVariations.count( it.first ) == 0 )
                {
                    throw std::runtime_error( "Error when estimating gravity field polynomial corrections of body " + bodyName +
                                              ", not polynomial term of order " + std::to_string( it.first ) + " found for sine coefficients" );
                }

                for( unsigned int i = 0; i < it.second.size( ); i++ )
                {
                    if( it.second.at( i ).first > polynomialVariationModel_->getMaximumDegree( ) ||
                        it.second.at( i ).first < polynomialVariationModel_->getMinimumDegree( ) )
                    {
                        throw std::runtime_error( "Error when estimating gravity field polynomial corrections of body " + bodyName +
                                                  " of polynomial term of order " + std::to_string( it.first ) +
                                                  ", no sine coefficient variation of degree " + std::to_string( it.second.at( i ).first ) + " found." );
                    }

                    if( it.second.at( i ).second > polynomialVariationModel_->getMaximumOrder( ) ||
                        it.second.at( i ).second < polynomialVariationModel_->getMinimumOrder( ) )
                    {
                        throw std::runtime_error( "Error when estimating gravity field polynomial corrections of body " + bodyName +
                                                  " of polynomial term of order " + std::to_string( it.first ) +
                                                  ", no sine coefficient variation of order " + std::to_string( it.second.at( i ).second ) + " found." );
                    }

                    sineCorrectionIndices_.push_back( std::make_tuple( it.first, it.second.at( i ).first, it.second.at( i ).second ) );
                    indexAndPowerPerSineBlockIndex_[ std::make_pair( it.second.at( i ).first, it.second.at( i ).second ) ].push_back(
                        std::make_pair( sineIndexCounter, it.first ) );
                    sineIndexCounter++;
                }

            }
        }

    //! Virtual destructor.
    ~PolynomialGravityFieldVariationsParameters( ) { }

    //! Pure virtual function to retrieve the value of the parameter
    /*!
     *  Pure virtual function to retrieve the value of the parameter
     *  \return Current value of parameter.
     */
    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd polynomialCorrections = Eigen::VectorXd::Zero( getParameterSize( ) );
        for( unsigned int i = 0; i < cosineCorrectionIndices_.size( ); i++ )
        {
            polynomialCorrections( i ) = polynomialVariationModel_->getCosineAmplitudesReference( ).at( std::get< 0 >( cosineCorrectionIndices_.at( i ) ) )
                ( std::get< 1 >( cosineCorrectionIndices_.at( i ) ) - polynomialVariationModel_->getMinimumDegree( ), std::get< 2 >( cosineCorrectionIndices_.at( i ) )  - polynomialVariationModel_->getMinimumOrder( ) );
        }

        for( unsigned int i = 0; i < sineCorrectionIndices_.size( ); i++ )
        {
            polynomialCorrections( i + cosineCorrectionIndices_.size( ) ) = polynomialVariationModel_->getSineAmplitudesReference( ).at( std::get< 0 >( sineCorrectionIndices_.at( i ) ) )
                ( std::get< 1 >( sineCorrectionIndices_.at( i ) ) - polynomialVariationModel_->getMinimumDegree( ), std::get< 2 >( cosineCorrectionIndices_.at( i ) ) - polynomialVariationModel_->getMinimumOrder( ) );
        }
        return polynomialCorrections;
    }

    //! Pure virtual function to (re)set the value of the parameter.
    /*!
     *  Pure virtual function to (re)set the value of the parameter.
     *  \param parameterValue to which the parameter is to be set.
     */
    void setParameterValue( const Eigen::VectorXd parameterValue )
    {
        if( parameterValue.rows( ) != getParameterSize( ) )
        {
            throw std::runtime_error( "Error when resetting polynomial gravity field variation parameter; sizes are incompatible" );
        }
        std::map< int, Eigen::MatrixXd > cosineVariations = polynomialVariationModel_->getCosineAmplitudes( );
        std::map< int, Eigen::MatrixXd > sineVariations = polynomialVariationModel_->getSineAmplitudes( );

        for( unsigned int i = 0; i < cosineCorrectionIndices_.size( ); i++ )
        {
            cosineVariations[ std::get< 0 >( cosineCorrectionIndices_.at( i ) ) ]( std::get< 1 >( cosineCorrectionIndices_.at( i ) ) - polynomialVariationModel_->getMinimumDegree( ),
                std::get< 2 >( cosineCorrectionIndices_.at( i ) ) - polynomialVariationModel_->getMinimumOrder( ) ) =
                parameterValue( i );
        }

        for( unsigned int i = 0; i < sineCorrectionIndices_.size( ); i++ )
        {
            sineVariations[ std::get< 0 >( sineCorrectionIndices_.at( i ) ) ]( std::get< 1 >( sineCorrectionIndices_.at( i ) ) - polynomialVariationModel_->getMinimumDegree( ),
                std::get< 2 >( sineCorrectionIndices_.at( i ) ) - polynomialVariationModel_->getMinimumOrder( ) ) =
                parameterValue( i + cosineCorrectionIndices_.size( ) );
        }
        polynomialVariationModel_->resetCosineAmplitudes( cosineVariations );
        polynomialVariationModel_->resetSineAmplitudes( sineVariations );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Pure virtual function to retrieve the size of the parameter (i.e. 1 for double parameters)
     *  \return Size of parameter value.
     */
    int getParameterSize( )
    {
        return cosineCorrectionIndices_.size( ) + sineCorrectionIndices_.size( );
    }

    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > getIndexAndPowerPerCosineBlockIndex( )
    {
        return indexAndPowerPerCosineBlockIndex_;
    }

    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > getIndexAndPowerPerSineBlockIndex( )
    {
        return indexAndPowerPerSineBlockIndex_;
    }
    std::shared_ptr< gravitation::PolynomialGravityFieldVariations > getPolynomialVariationModel( )
    {
        return polynomialVariationModel_;
    }


protected:

    std::shared_ptr< gravitation::PolynomialGravityFieldVariations > polynomialVariationModel_;

    std::vector< std::tuple< int, int, int > > cosineCorrectionIndices_;

    std::vector< std::tuple< int, int, int > > sineCorrectionIndices_;

    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerCosineBlockIndex_;

    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerSineBlockIndex_;

};


class PeriodicGravityFieldVariationsParameters: public EstimatableParameter< Eigen::VectorXd >
{

public:
    //! Constructor.
    /*!
     *  Constructor taking parameter name and associated body. All parameters are identified by a these two variables.
     *  Any additional information that may be required for uniquely defining a parameter is to be defined in the derived class.
     *  \param parameterName Enum value defining the type of the parameter.
     *  \param associatedBody Name of body associated with patameters
     *  \param pointOnBodyId Reference point on body associated with parameter (empty by default).
     */
    PeriodicGravityFieldVariationsParameters(
        const std::shared_ptr<gravitation::PeriodicGravityFieldVariations> periodicVariationModel,
        const std::map<int, std::vector<std::pair<int, int> > > &cosineBlockIndicesPerPeriod,
        const std::map<int, std::vector<std::pair<int, int> > > &sineBlockIndicesPerPeriod,
        const std::string &bodyName ) :
        EstimatableParameter<Eigen::VectorXd>( periodic_gravity_field_variation_amplitudes, bodyName ),
        periodicVariationModel_( periodicVariationModel )
    {
        std::vector< Eigen::MatrixXd > cosineVariations = periodicVariationModel->getCosineShAmplitudesCosineTime( );
        int cosineIndexCounter = 0;
        for( auto it : cosineBlockIndicesPerPeriod )
        {
            if( static_cast< int >( cosineVariations.size( ) ) < it.first )
            {
                throw std::runtime_error( "Error when estimating gravity field periodic corrections of body " + bodyName +
                                          ", not periodic term of order " + std::to_string( it.first ) + " found for cosine coefficients" );
            }

            for( unsigned int i = 0; i < it.second.size( ); i++ )
            {
                if( it.second.at( i ).first > periodicVariationModel_->getMaximumDegree( ) ||
                    it.second.at( i ).first < periodicVariationModel_->getMinimumDegree( ) )
                {
                    throw std::runtime_error( "Error when estimating gravity field periodic corrections of body " + bodyName +
                                              " of periodic term of period index " + std::to_string( it.first ) +
                                              ", no cosine coefficient variation of degree " + std::to_string( it.second.at( i ).first ) + " found." );
                }

                if( it.second.at( i ).second > periodicVariationModel_->getMaximumOrder( ) ||
                    it.second.at( i ).second < periodicVariationModel_->getMinimumOrder( ) )
                {
                    throw std::runtime_error( "Error when estimating gravity field periodic corrections of body " + bodyName +
                                              " of periodic term of period index " + std::to_string( it.first ) +
                                              ", no cosine coefficient variation of order " + std::to_string( it.second.at( i ).second ) + " found." );
                }
                cosineCorrectionIndices_.push_back( std::make_tuple( it.first, it.second.at( i ).first, it.second.at( i ).second ) );
                indexAndPowerPerCosineBlockIndex_[ std::make_pair( it.second.at( i ).first, it.second.at( i ).second ) ].push_back(
                    std::make_pair( cosineIndexCounter, it.first ) );
                cosineIndexCounter++;
            }
        }

        std::vector< Eigen::MatrixXd > sineVariations = periodicVariationModel->getSineShAmplitudesCosineTime( );
        int sineIndexCounter = 0;
        for( auto it : sineBlockIndicesPerPeriod )
        {
            if( static_cast< int >( sineVariations.size( ) ) < it.first )
            {
                throw std::runtime_error( "Error when estimating gravity field periodic corrections of body " + bodyName +
                                          ", not periodic term of order " + std::to_string( it.first ) + " found for sine coefficients" );
            }

            for( unsigned int i = 0; i < it.second.size( ); i++ )
            {
                if( it.second.at( i ).first > periodicVariationModel_->getMaximumDegree( ) ||
                    it.second.at( i ).first < periodicVariationModel_->getMinimumDegree( ) )
                {
                    throw std::runtime_error( "Error when estimating gravity field periodic corrections of body " + bodyName +
                                              " of periodic term of period index " + std::to_string( it.first ) +
                                              ", no sine coefficient variation of degree " + std::to_string( it.second.at( i ).first ) + " found." );
                }

                if( it.second.at( i ).second > periodicVariationModel_->getMaximumOrder( ) ||
                    it.second.at( i ).second < periodicVariationModel_->getMinimumOrder( ) )
                {
                    throw std::runtime_error( "Error when estimating gravity field periodic corrections of body " + bodyName +
                                              " of periodic term of period index " + std::to_string( it.first ) +
                                              ", no sine coefficient variation of order " + std::to_string( it.second.at( i ).second ) + " found." );
                }
                sineCorrectionIndices_.push_back( std::make_tuple( it.first, it.second.at( i ).first, it.second.at( i ).second ) );
                indexAndPowerPerSineBlockIndex_[ std::make_pair( it.second.at( i ).first, it.second.at( i ).second ) ].push_back(
                    std::make_pair( sineIndexCounter, it.first ) );
                sineIndexCounter++;
            }
        }
    }

    //! Virtual destructor.
    ~PeriodicGravityFieldVariationsParameters( ) { }

    //! Pure virtual function to retrieve the value of the parameter
    /*!
     *  Pure virtual function to retrieve the value of the parameter
     *  \return Current value of parameter.
     */
    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd periodicCorrections = Eigen::VectorXd::Zero( getParameterSize( ) );
        for( unsigned int i = 0; i < cosineCorrectionIndices_.size( ); i++ )
        {
            int periodIndex = std::get< 0 >( cosineCorrectionIndices_.at( i ) );
            int degree = std::get< 1 >( cosineCorrectionIndices_.at( i ) ) - periodicVariationModel_->getMinimumDegree( );
            int order = std::get< 2 >( cosineCorrectionIndices_.at( i ) ) - periodicVariationModel_->getMinimumOrder( );

            periodicCorrections( 2 * i ) = periodicVariationModel_->getCosineShAmplitudesCosineTime( ).at( periodIndex )( degree, order );
            periodicCorrections( 2 * i + 1 ) = periodicVariationModel_->getCosineShAmplitudesSineTime( ).at( periodIndex )( degree, order );
        }

        for( unsigned int i = 0; i < sineCorrectionIndices_.size( ); i++ )
        {
            int periodIndex = std::get< 0 >( sineCorrectionIndices_.at( i ) );
            int degree = std::get< 1 >( sineCorrectionIndices_.at( i ) ) - periodicVariationModel_->getMinimumDegree( );
            int order = std::get< 2 >( sineCorrectionIndices_.at( i ) ) - periodicVariationModel_->getMinimumOrder( );

            periodicCorrections( 2 * ( i + cosineCorrectionIndices_.size( ) ) ) = periodicVariationModel_->getSineShAmplitudesCosineTime( ).at( periodIndex )( degree, order );
            periodicCorrections( 2 * ( i + cosineCorrectionIndices_.size( ) ) + 1 ) = periodicVariationModel_->getSineShAmplitudesSineTime( ).at( periodIndex )( degree, order );
        }
        return periodicCorrections;
    }

    //! Pure virtual function to (re)set the value of the parameter.
    /*!
     *  Pure virtual function to (re)set the value of the parameter.
     *  \param parameterValue to which the parameter is to be set.
     */
    void setParameterValue( const Eigen::VectorXd parameterValue )
    {
        if( parameterValue.rows( ) != getParameterSize( ) )
        {
            throw std::runtime_error( "Error when resetting periodic gravity field variation parameter; sizes are incompatible" );
        }
        std::vector< Eigen::MatrixXd > cosineVariationsCosineTime = periodicVariationModel_->getCosineShAmplitudesCosineTime( );
        std::vector< Eigen::MatrixXd > cosineVariationsSineTime = periodicVariationModel_->getCosineShAmplitudesSineTime( );
        std::vector< Eigen::MatrixXd > sineVariationsCosineTime = periodicVariationModel_->getSineShAmplitudesCosineTime( );
        std::vector< Eigen::MatrixXd > sineVariationsSineTime = periodicVariationModel_->getSineShAmplitudesSineTime( );


        for( unsigned int i = 0; i < cosineCorrectionIndices_.size( ); i++ )
        {
            int periodIndex = std::get< 0 >( cosineCorrectionIndices_.at( i ) );
            int degree = std::get< 1 >( cosineCorrectionIndices_.at( i ) ) - periodicVariationModel_->getMinimumDegree( );
            int order = std::get< 2 >( cosineCorrectionIndices_.at( i ) ) - periodicVariationModel_->getMinimumOrder( );

            cosineVariationsCosineTime[ periodIndex ]( degree, order ) = parameterValue( 2 * i );
            cosineVariationsSineTime[ periodIndex ]( degree, order ) = parameterValue( 2 * i + 1 );
        }

        for( unsigned int i = 0; i < sineCorrectionIndices_.size( ); i++ )
        {
            int periodIndex = std::get< 0 >( sineCorrectionIndices_.at( i ) );
            int degree = std::get< 1 >( sineCorrectionIndices_.at( i ) ) - periodicVariationModel_->getMinimumDegree( );
            int order = std::get< 2 >( sineCorrectionIndices_.at( i ) ) - periodicVariationModel_->getMinimumOrder( );

            sineVariationsCosineTime[ periodIndex ]( degree, order ) = parameterValue( 2 * ( i + cosineCorrectionIndices_.size( ) ) );
            sineVariationsSineTime[ periodIndex ]( degree, order ) = parameterValue( 2 * ( i + cosineCorrectionIndices_.size( ) ) + 1 );
        }
        periodicVariationModel_->resetCosineShAmplitudesCosineTime( cosineVariationsCosineTime );
        periodicVariationModel_->resetCosineShAmplitudesSineTime( cosineVariationsSineTime );
        periodicVariationModel_->resetSineShAmplitudesCosineTime( sineVariationsCosineTime );
        periodicVariationModel_->resetSineShAmplitudesSineTime( sineVariationsSineTime );

    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Pure virtual function to retrieve the size of the parameter (i.e. 1 for double parameters)
     *  \return Size of parameter value.
     */
    int getParameterSize( )
    {
        return 2 * ( cosineCorrectionIndices_.size( ) + sineCorrectionIndices_.size( ) );
    }

    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > getIndexAndPowerPerCosineBlockIndex( )
    {
        return indexAndPowerPerCosineBlockIndex_;
    }

    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > getIndexAndPowerPerSineBlockIndex( )
    {
        return indexAndPowerPerSineBlockIndex_;
    }
    std::shared_ptr< gravitation::PeriodicGravityFieldVariations > getPeriodicVariationModel( )
    {
        return periodicVariationModel_;
    }


protected:

    std::shared_ptr< gravitation::PeriodicGravityFieldVariations > periodicVariationModel_;

    std::vector< std::tuple< int, int, int > > cosineCorrectionIndices_;

    std::vector< std::tuple< int, int, int > > sineCorrectionIndices_;

    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerCosineBlockIndex_;

    std::map< std::pair< int, int >, std::vector< std::pair< int, int > > > indexAndPowerPerSineBlockIndex_;

};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_GRAVITYFIELDVARIATIONPARAMETERS_H
