/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PROPAGATIONTERMINATIONDETAILS_H
#define TUDAT_PROPAGATIONTERMINATIONDETAILS_H

#include <iostream>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

#include <cereal/access.hpp>
#include <cereal/cereal.hpp>
#include <cereal/types/base_class.hpp>
#include <cereal/types/vector.hpp>
#include "tudat/io/serialization/file_io_declarations.h"

namespace tudat
{

namespace propagators
{

class HybridPropagationTerminationCondition;

//! Possible events that can trigger the termination of a propagation
enum PropagationTerminationReason {
    propagation_never_run,
    unknown_propagation_termination_reason,
    termination_condition_reached,
    runtime_error_caught_in_propagation,
    nan_or_inf_detected_in_state
};

inline std::string getPropagationTerminationReasonString( const PropagationTerminationReason terminationReason,
                                                          const bool exactFinalCondition = false )
{
    std::string reasonString;
    switch( terminationReason )
    {
        case propagation_never_run:
            reasonString = "Propagation was never run";
            break;
        case unknown_propagation_termination_reason:
            reasonString = "Unknown termination reason";
            break;
        case termination_condition_reached:
            if( !exactFinalCondition )
            {
                reasonString = "Propagation successful; termination condition exceeded";
            }
            else
            {
                reasonString = "Propagation successful; exact termination condition reached by iteration";
            }
            break;
        case runtime_error_caught_in_propagation:
            reasonString = "Exception caught";
            break;
        case nan_or_inf_detected_in_state:
            reasonString = "NaN or Inf detected in state";
            break;
        default:
            throw std::runtime_error( "Error when getting propagation termination reason string; type not recognized" );
    }
    return reasonString;
}

//! Class for storing details on the propagation termination
class PropagationTerminationDetails
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param propagationTerminationReason Reason for termination
     * \param terminationOnExactCondition True if exact termination condition is used, false if not, -1 if neither is relevant
     */
    PropagationTerminationDetails( const PropagationTerminationReason propagationTerminationReason,
                                   const bool terminationOnExactCondition = 0 ):
        propagationTerminationReason_( propagationTerminationReason ), terminationOnExactCondition_( terminationOnExactCondition )
    {}

    //! Destructor
    virtual ~PropagationTerminationDetails( ) {}

    //! Function to retrieve reason for termination
    /*!
     * Function to retrieve reason for termination
     * \return Reason for termination
     */
    PropagationTerminationReason getPropagationTerminationReason( )
    {
        return propagationTerminationReason_;
    }

    //! Function to retrieve boolean to denote whether exact termination conditions are used.
    /*!
     * Function to retrieve boolean to denote whether exact termination conditions are used.
     * \return Boolean to denote whether exact termination conditions are used.
     */
    bool getTerminationOnExactCondition( )
    {
        return terminationOnExactCondition_;
    }

    std::string getTerminationReasonString( )
    {
        std::string baseString = getPropagationTerminationReasonString( propagationTerminationReason_, terminationOnExactCondition_ );
        return baseString;
    }

    bool operator==( const PropagationTerminationDetails& rhs ) const
    {
        return equals( rhs );
    }
    bool operator!=( const PropagationTerminationDetails& rhs ) const
    {
        return !equals( rhs );
    }

protected:
    //! Reason for termination
    PropagationTerminationReason propagationTerminationReason_;

    //! Boolean to denote whether exact termination conditions are used.
    /*!
     *  Boolean to denote whether exact termination conditions are used. True if exact termination condition is used,
     *  false if not, -1 if neither is relevant.
     */
    bool terminationOnExactCondition_;

    //! Protected default constructor for deserialization (also accessible by derived classes)
    PropagationTerminationDetails( ): propagationTerminationReason_( propagation_never_run ), terminationOnExactCondition_( false ) {}

    virtual bool equals( const PropagationTerminationDetails& rhs ) const
    {
        return ( propagationTerminationReason_ == rhs.propagationTerminationReason_ ) &&
                ( terminationOnExactCondition_ == rhs.terminationOnExactCondition_ );
    }

public:
    //! Save termination details to a binary file
    TUDAT_DECLARE_FILE_IO_POLYMORPHIC( PropagationTerminationDetails )

private:
    friend class cereal::access;

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( CEREAL_NVP( propagationTerminationReason_ ) );
        ar( CEREAL_NVP( terminationOnExactCondition_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( CEREAL_NVP( propagationTerminationReason_ ) );
        ar( CEREAL_NVP( terminationOnExactCondition_ ) );
    }
};

//! Class for storing details on the propagation termination when using hybrid termination conditions
class PropagationTerminationDetailsFromHybridCondition : public PropagationTerminationDetails
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param terminationOnExactCondition True if exact termination condition is used, false if not.
     * \param terminationCondition Hybrid termination conditions that were used
     */
    PropagationTerminationDetailsFromHybridCondition( const bool terminationOnExactCondition,
                                                      const std::shared_ptr< HybridPropagationTerminationCondition > terminationCondition );

    //! Destructor
    ~PropagationTerminationDetailsFromHybridCondition( ) {}

    //! Function to retrieve list of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
    /*!
     * Function to retrieve list of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
     * \return List of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
     */
    std::vector< bool > getWasConditionMetWhenStopping( )
    {
        if( terminationOnExactCondition_ )
        {
            std::cerr << "Warning when retrieving list of conditions that were met in hybrid propagation "
                         "termination details. propagation was terminated on exact conditions using root finder, "
                         "list of conditions may not be reliable"
                      << std::endl;
        }
        return isConditionMetWhenStopping_;
    }

protected:
    bool equals( const PropagationTerminationDetails& rhs ) const override
    {
        if( !PropagationTerminationDetails::equals( rhs ) )
        {
            return false;
        }

        const auto* rhsCast = dynamic_cast< const PropagationTerminationDetailsFromHybridCondition* >( &rhs );
        if( rhsCast == nullptr )
        {
            return false;
        }

        return ( isConditionMetWhenStopping_ == rhsCast->isConditionMetWhenStopping_ );
    }

private:
    //! List of booleans, denoting for each of the constituent stopping conditions whether or not is was met.
    std::vector< bool > isConditionMetWhenStopping_;

    friend class cereal::access;

    PropagationTerminationDetailsFromHybridCondition( ): PropagationTerminationDetails( ) {}

    template< class Archive >
    void save( Archive& ar ) const
    {
        ar( cereal::base_class< PropagationTerminationDetails >( this ) );
        ar( CEREAL_NVP( isConditionMetWhenStopping_ ) );
    }

    template< class Archive >
    void load( Archive& ar )
    {
        ar( cereal::base_class< PropagationTerminationDetails >( this ) );
        ar( CEREAL_NVP( isConditionMetWhenStopping_ ) );
    }
};

}  // namespace propagators

}  // namespace tudat

#endif  // TUDAT_PROPAGATIONTERMINATIONDETAILS_H
