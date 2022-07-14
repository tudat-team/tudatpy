/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#include "tudat/astro/orbit_determination/observation_partials/differencedObservationPartial.h"
namespace tudat
{

namespace observation_partials
{


    void DifferencedObservablePartialScaling::update( const std::vector< Eigen::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd,
                 const Eigen::VectorXd currentObservation )
    {
        try
        {
            firstPartialScaling_->update(
                        utilities::getVectorEntries( linkEndStates, firstIndices_ ), utilities::getVectorEntries( times, firstIndices_ ),
                        fixedLinkEnd, Eigen::VectorXd::Constant( currentObservation.rows( ), TUDAT_NAN ) );
            secondPartialScaling_->update(
                        utilities::getVectorEntries( linkEndStates, secondIndices_ ), utilities::getVectorEntries( times, secondIndices_ ),
                        fixedLinkEnd, Eigen::VectorXd::Constant( currentObservation.rows( ), TUDAT_NAN ) );;
        }
        catch( const std::exception& caughtException )
        {
            std::string exceptionText = std::string( caughtException.what( ) );
            throw std::runtime_error( "Error when computing differenced observation partial scaling, error: " + exceptionText );
        }



    }

}

}


