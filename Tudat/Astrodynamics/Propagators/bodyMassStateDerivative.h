#ifndef BODYMASSSTATEDERIVATIVE_H
#define BODYMASSSTATEDERIVATIVE_H

/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include <vector>
#include <map>
#include <string>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"


namespace tudat
{

namespace propagators
{

template< typename StateScalarType = double, typename TimeType = double >
class BodyMassStateDerivative: public propagators::SingleStateTypeDerivative< StateScalarType, TimeType >
{
public:

    using propagators::SingleStateTypeDerivative< StateScalarType, TimeType >::calculateSystemStateDerivative;

    BodyMassStateDerivative(
            const std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > >& massRateModels,
            const std::vector< std::string >& bodiesToIntegrate ):
        propagators::SingleStateTypeDerivative< StateScalarType, TimeType >(
            propagators::body_mass_state ),
        massRateModels_( massRateModels ), bodiesToIntegrate_( bodiesToIntegrate ){ }

    //! Destructor
    virtual ~BodyMassStateDerivative( ){ }

    void calculateSystemStateDerivative(
                const TimeType time,
                const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& stateOfSystemToBeIntegrated,
                Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > > stateDerivative )
    {
        stateDerivative.setZero( );

        int currentIndex = 0;
        for( massRateModelIterator_ = massRateModels_.begin( );
             massRateModelIterator_ != massRateModels_.end( );
             massRateModelIterator_++ )
        {
            stateDerivative( currentIndex, 0 ) = static_cast< StateScalarType >(
                        massRateModelIterator_->second->getMassRate( ) );
            currentIndex++;
        }
    }

    void updateStateDerivativeModel( const TimeType currentTime )
    {
        // Reser all acceleration times (to allow multiple evaluations at same time, e.g. stage 2
        // and 3 in RK4 integrator)
        for( massRateModelIterator_ = massRateModels_.begin( );
             massRateModelIterator_ != massRateModels_.end( );
             massRateModelIterator_++ )
        {
            massRateModelIterator_->second->resetTime( TUDAT_NAN );
        }

        for( massRateModelIterator_ = massRateModels_.begin( );
             massRateModelIterator_ != massRateModels_.end( );
             massRateModelIterator_++ )
        {
            massRateModelIterator_->second->resetTime( static_cast< double >( currentTime ) );
        }
    }

    void convertCurrentStateToGlobalRepresentation(
                const Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 >& internalSolution, const TimeType& time,
                Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        currentCartesianLocalSoluton = internalSolution;
    }


    virtual Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const TimeType& time )
    {
        return outputSolution;
    }


    void convertToOutputSolution(
                const Eigen::Matrix< StateScalarType, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const TimeType& time,
                Eigen::Block< Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > currentCartesianLocalSoluton )
    {
        currentCartesianLocalSoluton = internalSolution;
    }

    virtual int getStateSize( )
    {
        return bodiesToIntegrate_.size( );
    }

private:

    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > > massRateModels_;

    std::map< std::string, boost::shared_ptr< basic_astrodynamics::MassRateModel > >::const_iterator massRateModelIterator_;

    std::vector< std::string > bodiesToIntegrate_;

};

} // namespace propagators

} // namespace tudat

#endif // BODYMASSSTATEDERIVATIVE_H
