/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONBIAS_H
#define TUDAT_OBSERVATIONBIAS_H

#include <iostream>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/make_shared.hpp>
#include <boost/function.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace observation_models
{


//! Base class (non-functional) for describing observation biases
/*!
 * Base class (non-functional) for describing observation biases. In this context, an observation bias denotes any deviation
 * from the ideal observable between two reference points.
 */
template< int ObservationSize = 1 >
class ObservationBias
{
public:

    //! Constructor
    ObservationBias( ){ }

    //! Destructor
    virtual ~ObservationBias( ){ }

    //! Pure virtual function to retrieve the observation bias.
    /*!
     * Pure virtual function to retrieve the observation bias as a function of the observation times and states (which are
     * typically computed by the ObservationModel)
     * \param linkEndTimes List of times at each link end during observation.
     * \param linkEndStates List of states at each link end during observation.
     * \return Observation bias at given times and states.
     */
    virtual Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates ) = 0;

    //! Function to return the size of the associated observation
    /*!
     * Function to return the size of the associated observation
     * \return Size of the associated observation
     */
    int getObservationSize( )
    {
        return ObservationSize;
    }
};

//! Class for a constant observation bias of a given size
template< int ObservationSize = 1 >
class ConstantObservationBias: public ObservationBias< ObservationSize >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param observationBias Constant (entry-wise) observation bias.
     */
    ConstantObservationBias( const Eigen::Matrix< double, ObservationSize, 1 > observationBias ):
        observationBias_( observationBias ){ }

    //! Destructor
    ~ConstantObservationBias( ){ }

    //! Function to retrieve the constant observation bias.
    /*!
     * Function to retrieve the constant observation bias.
     * \param linkEndTimes List of times at each link end during observation (unused).
     * \param linkEndStates List of states at each link end during observation (unused).
     * \return Constant observation bias.
     */
    Eigen::Matrix< double, ObservationSize, 1 > getObservationBias(
            const std::vector< double >& linkEndTimes,
            const std::vector< Eigen::Matrix< double, 6, 1 > >& linkEndStates )
    {
        return observationBias_;
    }


private:

    //! Constant (entry-wise) observation bias.
    Eigen::Matrix< double, ObservationSize, 1 > observationBias_;

};

}

}

#endif // TUDAT_OBSERVATIONMODEL_H
