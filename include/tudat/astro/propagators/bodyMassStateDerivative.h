/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_BODYMASSSTATEDERIVATIVE_H
#define TUDAT_BODYMASSSTATEDERIVATIVE_H

#include <vector>
#include <map>
#include <string>

#include <memory>
#include <functional>

#include "tudat/astro/basic_astro/massRateModel.h"
#include "tudat/astro/propagators/singleStateTypeDerivative.h"


namespace tudat
{

namespace propagators
{

//! Class for computing the derivative of the mass of a set of bodies.
/*!
 *  Class for computing the derivative of the mass of a set of bodies. Note that the local and global states are equal
 *  for mass propagation, both represent the physical body mass (in kg).
 */
template< typename StateScalarType = double, typename TimeType = double >
class BodyMassStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor, sets the mass rate models and the bodies that are to be propagated (single mass rate model per body).
     * \param massRateModels Map of model per body that is to be used for the mass rate computation.
     * \param bodiesToIntegrate List of bodies for which the mass is to be propagated. Note that this vector have
     * more entries than the massRateModels map, as a body's mass can be 'propagated' with no rate model (i.e. constant
     * mass).
     */
    BodyMassStateDerivative(
            const std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > >& massRateModels,
            const std::vector< std::string >& bodiesToIntegrate ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >(
            propagators::body_mass_state ),
        bodiesToIntegrate_( bodiesToIntegrate )
    {
        for( std::map< std::string, std::shared_ptr< basic_astrodynamics::MassRateModel > >::const_iterator modelIterator
             = massRateModels.begin( ); modelIterator != massRateModels.end( ); modelIterator++ )
        {
           massRateModels_[ modelIterator->first ].push_back( modelIterator->second );
        }

        for( unsigned int i = 0; i < bodiesToIntegrate_.size( ); i++ )
        {
            if( massRateModels_.count( bodiesToIntegrate_.at( i ) ) == 0 )
            {
                massRateModels_[ bodiesToIntegrate_.at( i ) ] =
                        std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > >( );
            }
        }

        verifyInput( );
    }

    //! Constructor
    /*!
     * Constructor, sets the mass rate models and the bodies that are to be propagated.
     * \param massRateModels Map of models per body that are to be used for the mass rate computation.
     * \param bodiesToIntegrate List of bodies for which the mass is to be propagated. Note that this vector have
     * more entries than the massRateModels map, as a body's mass can be 'propagated' with no rate model (i.e. constant
     * mass).
     */
    BodyMassStateDerivative(
            const std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > >&
            massRateModels,
            const std::vector< std::string >& bodiesToIntegrate ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >(
            propagators::body_mass_state ),
        massRateModels_( massRateModels ), bodiesToIntegrate_( bodiesToIntegrate )
    {
        for( unsigned int i = 0; i < bodiesToIntegrate_.size( ); i++ )
        {
            if( massRateModels_.count( bodiesToIntegrate_.at( i ) ) == 0 )
            {
                massRateModels_[ bodiesToIntegrate_.at( i ) ] =
                        std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > >( );
            }
        }

        verifyInput( );
    }


    //! Destructor
    virtual ~BodyMassStateDerivative( ){ }

    //! Calculates the state derivative of the system of equations for the mass dynamics
    /*!
     * Calculates the state derivative of the system of equations for the mass dynamics
     * The environment and acceleration models (updateStateDerivativeModel) must be
     * updated before calling this function.
     * \param time Time at which the state derivative is to be calculated.
     * \param stateOfSystemToBeIntegrated Current masses of the bodies that are propagated
     * \param stateDerivative Mass rates of the bodies for which the mass is propagated, in the same order as
     * bodiesToIntegrate_
     */
    void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative.setZero( );

        // Iterate over all mass rate models, retrieve value and put into corresponding entry.
        int currentIndex = 0;
        for( massRateModelIterator_ = massRateModels_.begin( );
             massRateModelIterator_ != massRateModels_.end( );
             massRateModelIterator_++ )
        {
            stateDerivative( currentIndex, 0 ) = 0.0;
            for( unsigned int i = 0; i < massRateModelIterator_->second.size( ); i++ )
            {
                stateDerivative( currentIndex, 0 ) += static_cast< StateScalarType >(
                            massRateModelIterator_->second.at ( i )->getMassRate( ) );

                currentIndex++;
            }

        }
    }

    //! Function to clear reference/cached values of body mass state derivative model
    /*!
     * Function to clear reference/cached values of body mass state derivative model. All mass rate models' current times
     * are reset to ensure that they are all recalculated.
     */
    void clearStateDerivativeModel( )
    {
        // Reset all mass rate times (to allow multiple evaluations at same time, e.g. stage 2  and 3 in RK4 integrator)
        for( massRateModelIterator_ = massRateModels_.begin( );
             massRateModelIterator_ != massRateModels_.end( );
             massRateModelIterator_++ )
        {
            for( unsigned int i = 0; i < massRateModelIterator_->second.size( ); i++ )
            {
                massRateModelIterator_->second.at ( i )->resetTime( TUDAT_NAN );
            }
        }
    }

    //! Function to update the mass state derivative model to the current time.
    /*!
     * Function to update the mass state derivative model to the urrent time.
     * cNote that this function only updates the state derivative model itself, the
     * environment models must be updated before calling this function
     * \param currentTime Time to which the mass state derivative is to be updated.
     */
    void updateStateDerivativeModel( const TimeType currentTime )
    {


        // Update local variables of mass rate model objects.
        for( massRateModelIterator_ = massRateModels_.begin( );
             massRateModelIterator_ != massRateModels_.end( );
             massRateModelIterator_++ )
        {
            for( unsigned int i = 0; i < massRateModelIterator_->second.size( ); i++ )
            {
                massRateModelIterator_->second.at ( i )->updateMembers( static_cast< double >( currentTime ) );
            }
        }
    }

    //! Function included for compatibility purposes with base class, local and global representation is equal for mass rate
    //! model. Function returns (by reference)  input internalSolution.
    void convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        currentCartesianLocalSoluton = internalSolution;
    }

    //! Function included for compatibility purposes with base class, input and output representation is equal for mass rate
    //! model. Function returns input outputSolution.
    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time )
    {
        return outputSolution;
    }

    //! Function included for compatibility purposes with base class, input and output representation is equal for mass rate
    //! model. Function returns  (by reference) input internalSolution.
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution,
            const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        currentCartesianLocalSoluton = internalSolution;
    }

    //! Function to get the total size of the state of propagated masses.
    /*!
     * Function to get the total size of the state of propagated masses. Equal to number of bodies for which the mass
     * is propagated.
     * \return Size of propagated mass state.
     */
    virtual int getConventionalStateSize( )
    {
        return bodiesToIntegrate_.size( );
    }

    //! Get list of bodies for which the mass is to be propagated.
    /*!
     * Get list of bodies for which the mass is to be propagated.
     * \return List of bodies for which the mass is to be propagated.
     */
    std::vector< std::string > getBodiesToIntegrate( )
    {
        return bodiesToIntegrate_;
    }

    //! Function to return the current mass rate for a single body
    /*!
     *  Function to return the current mass rate for a single body
     * \param bodyName Name of body fo which the mass rate is to be returned.
     * \return Current mass rate for requested body
     */
    double getTotalMassRateForBody( const std::string bodyName )
    {
        // Check if body is propagated.
        double totalMassRate = 0.0;

        if( std::find( bodiesToIntegrate_.begin( ),
                       bodiesToIntegrate_.end( ),
                       bodyName ) == bodiesToIntegrate_.end( ) )
        {
            std::string errorMessage = "Error when getting total mass rate for body " + bodyName +
                    ", no such acceleration is found";
            throw std::runtime_error( errorMessage );
        }
        else
        {
            if( massRateModels_.count( bodyName ) != 0 )
            {
                for( unsigned int i = 0; i < massRateModels_.at( bodyName ).size( ); i++  )
                {
                    totalMassRate += massRateModels_.at( bodyName ).at ( i )->getMassRate( );
                }
            }
        }
        return totalMassRate;
    }

private:

    void verifyInput( )
    {
        for( unsigned int i = 0; i < bodiesToIntegrate_.size( ); i++ )
        {
            if( massRateModels_.count( bodiesToIntegrate_.at( i ) ) == 0 )
            {
                throw std::runtime_error( "Error, requested propagation of mass of body " +
                                          bodiesToIntegrate_.at( i ) +
                                          ", but no mass rate models provided" );
            }
        }

        for( auto it : massRateModels_ )
        {
            if( std::find( bodiesToIntegrate_.begin( ),
                           bodiesToIntegrate_.end( ),
                           it.first ) == bodiesToIntegrate_.end( ) )
            {
                throw std::runtime_error( "Error, provided mass rate models for body " +
                                          it.first +
                                          ", but this body is not included in list of bodies for which mass is to be propagated." );
            }
        }
    }

    //! Map of models per body that are to be used for the mass rate computation.
    std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > > massRateModels_;

    //! Predefined iterator to save (de-)allocation time.
    std::map< std::string, std::vector< std::shared_ptr< basic_astrodynamics::MassRateModel > > >::const_iterator massRateModelIterator_;

    //! List of bodies for which the mass is to be propagated.
    /*!
     * List of bodies for which the mass is to be propagated. Note that this vector have
     * more entries than the massRateModels map, as a body's mass can be 'propagated' with no rate model (i.e. constant
     * mass).
     */
    std::vector< std::string > bodiesToIntegrate_;

};

} // namespace propagators

} // namespace tudat

#endif // TUDAT_BODYMASSSTATEDERIVATIVE_H
