/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EULERANGLEOBSERVABLEPARTIALS_H
#define TUDAT_EULERANGLEOBSERVABLEPARTIALS_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <functional>

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Mathematics/BasicMathematics/rotationRepresentations.h"

namespace tudat
{

namespace observation_partials
{


//! Class to compute the partial derivatives of a 3-1-3 Euler angle observable.
class EulerAngleObervationPartialWrtCurrentRotationalState: public ObservationPartial< 3 >
{

public:

    //! Local typedef for return type (list of partial matrices and associated evaluation times).
    typedef std::vector< std::pair< Eigen::Matrix< double, 3, Eigen::Dynamic >, double > >
    EulerAngleObservationPartialReturnType;

    //! Constructor
    /*!
     * Constructor
     * \param parameterIdentifier Id of parameter for which instance of class computes partial derivatives.
     */
    EulerAngleObervationPartialWrtCurrentRotationalState(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ):
        ObservationPartial< 3 >( parameterIdentifier ){ }

    //! Destructor
    ~EulerAngleObervationPartialWrtCurrentRotationalState( ) { }

    //! Function to calculate the observation partial(s) at required time.
    /*!
     *  Function to calculate the observation partial(s) at required time. Time is
     *  typically obtained from evaluation of observation model.
     *  \param states Link end state. (Not used; vector contains NaN entry)
     *  \param times Link end time (size 1 std::vector).
     *  \param linkEndOfFixedTime Link end that is kept fixed when computing the observable.
     *  \param currentObservation Value of the observation for which the partial is to be computed (default NaN for
     *  compatibility purposes)
     *  \return Vector of pairs containing partial values and associated times.
     */
    EulerAngleObservationPartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime,
            const Eigen::Vector3d& currentObservation = Eigen::Vector3d::Constant( TUDAT_NAN ) )
    {        
        EulerAngleObservationPartialReturnType returnPartial;

        currentTime_ = times.at( 0 );
        currentPartial_.setZero( );
        currentPartial_.block( 0, 0, 3, 4 ) =
                -basic_mathematics::calculateEulerAngle313WrtQuaternionPartialFromEulerAngles(
                                                            currentObservation );
        returnPartial.push_back( std::make_pair( currentPartial_, currentTime_ ) );
        return returnPartial;
    }

protected:

    //! Partial derivative, as computed by last call to calculatePartial function
    Eigen::Matrix< double, 3, 7 > currentPartial_;

    //! Pre-declared time variable to be used in calculatePartial function.
    double currentTime_;
};

}

}


#endif // TUDAT_EULERANGLEOBSERVABLEPARTIALS_H
