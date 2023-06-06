/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SOLARCORONACORRECTION_H
#define TUDAT_SOLARCORONACORRECTION_H

#include "tudat/astro/observation_models/observableTypes.h"
#include "tudat/astro/observation_models/corrections/lightTimeCorrection.h"

namespace tudat
{

namespace observation_models
{

class SolarCoronaCorrection: public LightTimeCorrection
{
public:

    SolarCoronaCorrection(
            const LightTimeCorrectionType lightTimeCorrectionType,
            const ObservableType baseObservableType ):
        LightTimeCorrection( lightTimeCorrectionType )
    {
        if ( isRadiometricObservableType( baseObservableType ) )
        {
            if ( isGroupVelocityBasedObservableType( baseObservableType ) )
            {
                sign_ = 1;
            }
            else if ( isPhaseVelocityBasedObservableType( baseObservableType ) )
            {
                sign_ = -1;
            }
            else
            {
                throw std::runtime_error( "Error when creating solar corona correction: radiometric correction not "
                                          "recognized." );
            }
        }
        else
        {
            throw std::runtime_error( "Error when creating solar corona correction: correction is only valid for "
                                      "radiometric types." );
        }
    }

protected:

    double computeMinimumDistanceOfLineOfSight(
            Eigen::Vector6d transmitterState,
            Eigen::Vector6d receiverState,
            Eigen::Vector6d sunState );

    // Sign of the correction (+1 or -1)
    int sign_;

private:

};

class AndersonSolarCoronaCorrection: public SolarCoronaCorrection
{
public:

private:

};


} // namespace observation_models

} // namespace tudat

#endif //TUDAT_SOLARCORONACORRECTION_H
