/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_LIGHTTIMECORRECTION_H
#define TUDAT_LIGHTTIMECORRECTION_H

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{

namespace observation_models
{

//! Enum defining different types of light time corrections.
enum LightTimeCorrectionType
{
    first_order_relativistic
};

class LightTimeCorrection
{
public:
    LightTimeCorrection( const LightTimeCorrectionType lightTimeCorrectionType ):
        lightTimeCorrectionType_( lightTimeCorrectionType ){ }

    virtual ~LightTimeCorrection( ){ }

    virtual double calculateLightTimeCorrection(
            const basic_mathematics::Vector6d& transmitterState,
            const basic_mathematics::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) = 0;

    LightTimeCorrectionType getLightTimeCorrectionType( )
    {
        return lightTimeCorrectionType_;
    }

protected:
    LightTimeCorrectionType lightTimeCorrectionType_;

};

}

}

#endif // TUDAT_LIGHTTIMECORRECTION_H
