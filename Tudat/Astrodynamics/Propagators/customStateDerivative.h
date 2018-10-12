/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_CUSTOMSTATEDERIVATIVE_H
#define TUDAT_CUSTOMSTATEDERIVATIVE_H

#include <functional>

#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"

namespace tudat
{

namespace propagators
{

//! Model to compute the derivative of a custom state (i.e. a state for which the physical significance is not 'known'
//! to the rest of the code
template< typename StateScalarType = double, typename TimeType = double >
class CustomStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:

    typedef Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > StateVectorType;

    //! Constructor
    /*!
     * Constructor
     * \param stateDerivativeModel Function to compute the state derivative, as a function of current time and state.
     * \param stateSize Size of the custom state that is propagated.
     */
    CustomStateDerivative(
            const std::function< StateVectorType( const TimeType, const StateVectorType& )> stateDerivativeModel,
            const int stateSize ):
        SingleStateTypeDerivative< StateScalarType, TimeType >( custom_state ),
        stateDerivativeModel_( stateDerivativeModel ), stateSize_( stateSize ){ }

    //! Calculates the custom state derivative
    /*!
     * Calculates the custom state derivative
     * \param time Time at which the state derivative is to be calculated.
     * \param stateOfSystemToBeIntegrated Current custom states
     * \param stateDerivative State derivative of custom states (returned by reference).
     */
    void calculateSystemStateDerivative(
            const TimeType time,
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative = stateDerivativeModel_( time, stateOfSystemToBeIntegrated );

    }

    //! Function included for consistency, not used in this derived class.
    void clearStateDerivativeModel( )
    {  }

    //! Function included for consistency, not used in this derived class.
    void updateStateDerivativeModel( const TimeType currentTime )
    { }

    //! Function included for compatibility purposes with base class, local and global representation is equal for custom
    //! model. Function returns (by reference)  input internalSolution.
    void convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentLocalSolution )
    {
        currentLocalSolution = internalSolution;
    }

    //! Function included for compatibility purposes with base class, input and output representation is equal for custom
    //! model. Function returns input outputSolution.
    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time )
    {
        return outputSolution;
    }

    //! Function included for compatibility purposes with base class, input and output representation is equal for custom
    //! model. Function returns  (by reference) input internalSolution.
    void convertToOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution,
            const TimeType& time,
            Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentLocalSolution )
    {
        currentLocalSolution = internalSolution;
    }

    //! Function to get the total size of the state of propagated masses.
    /*!
     * Function to get the total size of the state of propagated masses.
     * \return Size of propagated custom state.
     */
    virtual int getConventionalStateSize( )
    {
        return stateSize_;
    }


private:

    //! Function to compute the state derivative, as a function of current time and state
    std::function< StateVectorType( const TimeType, const StateVectorType& )> stateDerivativeModel_;

    //! Size of the custom state that is propagated.
    int stateSize_;
};

} // namespace propagators

} // namespae tudat
#endif // TUDAT_CUSTOMSTATEDERIVATIVE_H
