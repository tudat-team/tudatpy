/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_GRAVITYDERIVATIVE_H
#define TUDAT_GRAVITYDERIVATIVE_H

#include <vector>
#include <map>
#include <string>

#include <memory>
#include <functional>

#include "tudat/astro/basic_astro/gravityDeformationModel.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"

namespace tudat
{

namespace propagators
{

//! Class for computing the derivative of the gravity field of a set of bodies.
/*!
 *  Class for computing the derivative of the gravity field of a set of bodies.
 */
template< typename StateScalarType = double, typename TimeType = double >
class GravityStateDerivative : public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:
    //! Constructor
    /*!
     * Constructor, sets the gravity deformation models and the bodies that are to be propagated
     * (single gravity deformation model per body).
     * \param massRateModels Map of model per body that is to be used for the gravity deformation computation.
     * \param bodiesToIntegrate List of bodies for which the gravity coefficients are to be propagated.
     */
    GravityStateDerivative(
            const std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > >& gravityDeformationModels,
            const std::vector< std::string >& bodiesToIntegrate ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >( propagators::gravity_deformation_state, bodiesToIntegrate ),
        bodiesToIntegrate_( bodiesToIntegrate )
    {
        for( std::map< std::string, std::shared_ptr< basic_astrodynamics::GravityDeformationModel > >::const_iterator modelIterator =
                     gravityDeformationModels.begin( );
             modelIterator != gravityDeformationModels.end( );
             modelIterator++ )
        {
            gravityDeformationModels_[ modelIterator->first ].push_back( modelIterator->second );
        }

        for( unsigned int i = 0; i < bodiesToIntegrate_.size( ); i++ )
        {
            if( gravityDeformationModels_.count( bodiesToIntegrate_.at( i ) ) == 0 )
            {
                std::cerr << "Warning, propagating gravity coefficients of body " + bodiesToIntegrate_.at( i ) +
                                " without any gravity deformation model."
                          << std::endl;
                gravityDeformationModels_[ bodiesToIntegrate_.at( i ) ] =
                        std::vector< std::shared_ptr< basic_astrodynamics::GravityDeformationModel > >( );
            }
        }

        verifyInput( );
    }

    //! Constructor
    /*!
     * Constructor, sets the gravity deformation models and the bodies that are to be propagated.
     * \param massRateModels Map of models per body that are to be used for the gravity deformation computation.
     * \param bodiesToIntegrate List of bodies for which the gravity coefficients are to be propagated.
     */
    GravityStateDerivative( const std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > >&
                                    gravityDeformationModels,
                            const std::vector< std::string >& bodiesToIntegrate ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >( propagators::gravity_deformation_state ),
        gravityDeformationModels_( gravityDeformationModels ), bodiesToIntegrate_( bodiesToIntegrate )
    {
        for( unsigned int i = 0; i < bodiesToIntegrate_.size( ); i++ )
        {
            if( gravityDeformationModels_.count( bodiesToIntegrate_.at( i ) ) == 0 )
            {
                std::cerr << "Warning, propagating gravity coefficients of body " + bodiesToIntegrate_.at( i ) +
                                " without any gravity deformation model."
                          << std::endl;
                gravityDeformationModels_[ bodiesToIntegrate_.at( i ) ] =
                        std::vector< std::shared_ptr< basic_astrodynamics::GravityDeformationModel > >( );
            }
        }

        verifyInput( );
    }

    //! Destructor
    virtual ~GravityStateDerivative( ) {}

    //! Calculates the state derivative of the system of equations for the gravity field coefficients
    /*!
     * Calculates the state derivative of the system of equations for the gravity field coefficients
     * The environment and gravity deformation models (updateStateDerivativeModel) must be
     * updated before calling this function.
     * \param time Time at which the state derivative is to be calculated.
     * \param stateOfSystemToBeIntegrated Current gravity coefficients of the bodies that are propagated
     * \param stateDerivative Gravity deformation of the bodies for which the gravity coefficients are
     * propagated, in the same order as bodiesToIntegrate_
     */
    void calculateSystemStateDerivative( const TimeType time,
                                         const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
                                         Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative.setZero( );

        // Iterate over all gravity deformation models, retrieve value and put into corresponding entry.
        for( gravityDeformationModelIterator_ = gravityDeformationModels_.begin( );
             gravityDeformationModelIterator_ != gravityDeformationModels_.end( );
             gravityDeformationModelIterator_++ )
        {
            int currentIndex = -1;
            auto it = std::find( bodiesToIntegrate_.begin( ), bodiesToIntegrate_.end( ), gravityDeformationModelIterator_->first );
            if( it != bodiesToIntegrate_.end( ) ) currentIndex = static_cast< int >( it - bodiesToIntegrate_.begin( ) );

            stateDerivative( currentIndex, 0 ) = 0.0;
            for( unsigned int i = 0; i < gravityDeformationModelIterator_->second.size( ); i++ )
            {
                stateDerivative.block( currentIndex * 5, 0, 5, 1 ) +=
                        gravityDeformationModelIterator_->second.at( i )->getDeformation( ).template cast< StateScalarType >( );
            }
        }
    }

    //! Function to clear reference/cached values of gravity state derivative model
    /*!
     * Function to clear reference/cached values of gravity state derivative model. All gravity deformation models' current times
     * are reset to ensure that they are all recalculated.
     */
    void clearStateDerivativeModel( )
    {
        // Reset all gravity deformation times (to allow multiple evaluations at same time, e.g. stage 2  and 3 in RK4 integrator)
        for( gravityDeformationModelIterator_ = gravityDeformationModels_.begin( );
             gravityDeformationModelIterator_ != gravityDeformationModels_.end( );
             gravityDeformationModelIterator_++ )
        {
            for( unsigned int i = 0; i < gravityDeformationModelIterator_->second.size( ); i++ )
            {
                gravityDeformationModelIterator_->second.at( i )->resetCurrentTime( );
            }
        }
    }

    //! Function to update the gravity state derivative model to the current time.
    /*!
     * Function to update the gravity state derivative model to the current time.
     * Note that this function only updates the state derivative model itself, the
     * environment models must be updated before calling this function
     * \param currentTime Time to which the gravity state derivative is to be updated.
     */
    void updateStateDerivativeModel( const TimeType currentTime )
    {
        // Update local variables of gravity deformation model objects.
        for( gravityDeformationModelIterator_ = gravityDeformationModels_.begin( );
             gravityDeformationModelIterator_ != gravityDeformationModels_.end( );
             gravityDeformationModelIterator_++ )
        {
            for( unsigned int i = 0; i < gravityDeformationModelIterator_->second.size( ); i++ )
            {
                gravityDeformationModelIterator_->second.at( i )->updateMembers( static_cast< double >( currentTime ) );
            }
        }
    }

    //! Function included for compatibility purposes with base class, local and global representation is equal for gravity
    //! deformation model. Function returns (by reference) input internalSolution.
    void convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution,
            const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        currentCartesianLocalSoluton = internalSolution;
    }

    //! Function included for compatibility purposes with base class, input and output representation is equal for gravity
    //! deformation rate model. Function returns input outputSolution.
    Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution,
            const TimeType& time )
    {
        return outputSolution;
    }

    //! Function included for compatibility purposes with base class, input and output representation is equal for gravity
    //! deformation model. Function returns  (by reference) input internalSolution.
    void convertToOutputSolution( const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution,
                                  const TimeType& time,
                                  Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSolution )
    {
        currentCartesianLocalSolution = internalSolution;
    }

    //! Function to get the total size of the state of propagated gravity coefficients.
    /*!
     * Function to get the total size of the state of propagated gravity coeficients.
     * Equal to number of bodies for which the gravity coefficients are propagated.
     * \return Size of propagated gravity coefficients.
     */
    int getConventionalStateSize( )
    {
        return 5 * bodiesToIntegrate_.size( );
    }

    //! Get list of bodies for which the gravity coefficients are to be propagated.
    /*!
     * Get list of bodies for which the gravity coefficients are to be propagated.
     * \return List of bodies for which the gravity coefficients are to be propagated.
     */
    std::vector< std::string > getBodiesToIntegrate( )
    {
        return bodiesToIntegrate_;
    }

    //! Function to return the current gravity deformation for a single body
    /*!
     *  Function to return the current gravity deformation for a single body
     * \param bodyName Name of body for which the gravity deformation is to be returned.
     * \return Current gravity deformation for requested body
     */
    Eigen::VectorXd getTotalGravityDeformationForBody( const std::string bodyName )
    {
        // Check if body is propagated.
        Eigen::VectorXd totalGravityDeformation = Eigen::VectorXd::Zero( 5 );

        if( std::find( bodiesToIntegrate_.begin( ), bodiesToIntegrate_.end( ), bodyName ) == bodiesToIntegrate_.end( ) )
        {
            std::string errorMessage =
                    "Error when getting total gravity deformation for body " + bodyName + ", no such deformation is found";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( gravityDeformationModels_.count( bodyName ) != 0 )
            {
                for( unsigned int i = 0; i < gravityDeformationModels_.at( bodyName ).size( ); i++ )
                {
                    totalGravityDeformation += gravityDeformationModels_.at( bodyName ).at( i )->getDeformation( );
                }
            }
        }
        return totalGravityDeformation;
    }

    const basic_astrodynamics::GravityDeformationModelMap& getGravityDeformationModels( ) const
    {
        return gravityDeformationModels_;
    }

    bool isStateToBePostProcessed( )
    {
        return true;
    }

private:
    void verifyInput( )
    {
        for( unsigned int i = 0; i < bodiesToIntegrate_.size( ); i++ )
        {
            if( gravityDeformationModels_.count( bodiesToIntegrate_.at( i ) ) == 0 )
            {
                throw std::runtime_error( "Error, requested propagation of gravity deformation of body " + bodiesToIntegrate_.at( i ) +
                                          ", but no gravity deformation model provided" );
            }
        }

        for( auto it : gravityDeformationModels_ )
        {
            if( std::find( bodiesToIntegrate_.begin( ), bodiesToIntegrate_.end( ), it.first ) == bodiesToIntegrate_.end( ) )
            {
                throw std::runtime_error(
                        "Error, provided gravity deformation models for body " + it.first +
                        ", but this body is not included in list of bodies for which the gravity coefficients are to be propagated." );
            }
        }
    }

    //! Map of models per body that are to be used for the gravity field deformation computation.
    std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > > gravityDeformationModels_;

    //! Predefined iterator to save (de-)allocation time.
    std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::GravityDeformationModel > > >::const_iterator
            gravityDeformationModelIterator_;

    //! List of bodies for which the gravity field coefficients are to be propagated.
    /*!
     * List of bodies for which the gravity field coefficients are to be propagated.
     */
    std::vector< std::string > bodiesToIntegrate_;
};

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_GRAVITYDERIVATIVE_H
