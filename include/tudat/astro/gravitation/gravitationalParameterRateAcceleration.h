/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *          Genova, A., Mazarico, E., Goossens, S., Lemoine, F. G., Neumann, G. A., Smith, D. E., & Zuber, M. T. (2018).
 *          Solar system expansion and strong equivalence principle as seen by the NASA MESSENGER mission.
 *          Nature Communications, 9, 289. DOI: https://doi.org/10.1038/s41467-017-02558-1
 */

#ifndef TUDAT_GRAVITATIONALPARAMETERRATEACCELERATION_H
#define TUDAT_GRAVITATIONALPARAMETERRATEACCELERATION_H

#include <functional>

#include "tudat/astro/basic_astro/accelerationModel.h"

namespace tudat
{
namespace gravitation
{

//! Class for the point-mass acceleration due to a linear variation of the gravitational parameter of the exerting body.
/*!
 * Point-mass acceleration due to the change of the gravitational parameter mu of the exerting body since a reference epoch
 * t0, with a constant relative rate rho = (d mu / dt) / mu (Genova et al., 2018, Eq. 11):
 * a = -mu rho (t - t0) r / |r|^3, with r the position of the accelerated body w.r.t. the exerting body. It is to be used next to
 * a point-mass (or EIH) acceleration with the gravitational parameter mu at t0.
 *
 * This class introduces the time variation of mu in this acceleration ONLY. The gravitational parameter of the gravity field of
 * the exerting body remains constant (mu at t0), and that constant value is used by all other models, such as other accelerations
 * (including the point-mass or EIH acceleration next to this one), light-time corrections and the mass of the body.
 */
class GravitationalParameterRateAcceleration : public basic_astrodynamics::AccelerationModel< Eigen::Vector3d >
{
public:
    //! Constructor
    /*!
     * Constructor
     * \param acceleratedBodyPositionFunction Function returning the position of the body undergoing the acceleration.
     * \param exertingBodyPositionFunction Function returning the position of the body exerting the acceleration.
     * \param gravitationalParameterFunction Function returning the gravitational parameter of the exerting body (at the
     * reference epoch).
     * \param relativeGravitationalParameterRate Relative rate rho = (d mu / dt) / mu of the gravitational parameter [1/s].
     * \param referenceEpoch Epoch t0 at which the gravitational parameter is equal to that of gravitationalParameterFunction.
     */
    GravitationalParameterRateAcceleration( const std::function< Eigen::Vector3d( ) >& acceleratedBodyPositionFunction,
                                            const std::function< Eigen::Vector3d( ) >& exertingBodyPositionFunction,
                                            const std::function< double( ) >& gravitationalParameterFunction,
                                            const double relativeGravitationalParameterRate,
                                            const double referenceEpoch ):
        acceleratedBodyPositionFunction_( acceleratedBodyPositionFunction ), exertingBodyPositionFunction_( exertingBodyPositionFunction ),
        gravitationalParameterFunction_( gravitationalParameterFunction ),
        relativeGravitationalParameterRate_( relativeGravitationalParameterRate ), referenceEpoch_( referenceEpoch )
    {}

    void updateMembers( const double currentTime ) override
    {
        if( this->currentTime_ != currentTime )
        {
            currentRelativePosition_ = acceleratedBodyPositionFunction_( ) - exertingBodyPositionFunction_( );
            currentTimeSinceReferenceEpoch_ = currentTime - referenceEpoch_;
            this->currentAcceleration_ = -gravitationalParameterFunction_( ) * relativeGravitationalParameterRate_ *
                    currentTimeSinceReferenceEpoch_ * currentRelativePosition_ / std::pow( currentRelativePosition_.norm( ), 3 );
            this->currentTime_ = currentTime;
        }
    }

    Eigen::Vector3d getCurrentRelativePosition( ) const
    {
        return currentRelativePosition_;
    }

    double getCurrentTimeSinceReferenceEpoch( ) const
    {
        return currentTimeSinceReferenceEpoch_;
    }

    std::function< double( ) > getGravitationalParameterFunction( ) const
    {
        return gravitationalParameterFunction_;
    }

    double getRelativeGravitationalParameterRate( ) const
    {
        return relativeGravitationalParameterRate_;
    }

    void setRelativeGravitationalParameterRate( const double relativeGravitationalParameterRate )
    {
        relativeGravitationalParameterRate_ = relativeGravitationalParameterRate;
        this->currentTime_ = TUDAT_NAN;
    }

private:
    std::function< Eigen::Vector3d( ) > acceleratedBodyPositionFunction_;

    std::function< Eigen::Vector3d( ) > exertingBodyPositionFunction_;

    std::function< double( ) > gravitationalParameterFunction_;

    //! Relative rate (d mu / dt) / mu of the gravitational parameter of the exerting body [1/s]
    double relativeGravitationalParameterRate_;

    double referenceEpoch_;

    Eigen::Vector3d currentRelativePosition_;

    double currentTimeSinceReferenceEpoch_;
};

}  // namespace gravitation
}  // namespace tudat
#endif  // TUDAT_GRAVITATIONALPARAMETERRATEACCELERATION_H
