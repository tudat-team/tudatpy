/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#include "Tudat/Astrodynamics/EarthOrientation/terrestrialTimeScaleConverter.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"
#include "Tudat/Mathematics/Interpolators/jumpDataLinearInterpolator.h"

namespace tudat
{

namespace earth_orientation
{

//! Function to get current time list at double precision
template< >
CurrentTimes< double >& TerrestrialTimeScaleConverter::getCurrentTimeList< double >( )
{
    return currentTimes_;
}

//! Function to get current time list at Time precision
template< >
CurrentTimes< Time >& TerrestrialTimeScaleConverter::getCurrentTimeList< Time >( )
{
    return currentTimesSplit_;
}

//! Function to get ground station position used on last call to updateTimes function at double precision
template< >
Eigen::Vector3d& TerrestrialTimeScaleConverter::getPreviousGroundStationPosition< double >( )
{
    return previousEarthFixedPosition_;
}

//! Function to get ground station position used on last call to updateTimes function at Time precision
template< >
Eigen::Vector3d& TerrestrialTimeScaleConverter::getPreviousGroundStationPosition< Time >( )
{
    return previousEarthFixedPositionSplit_;
}

//! Function to reset ground station position used on last call to updateTimes function at double precision
template< >
void TerrestrialTimeScaleConverter::setCurrentGroundStation< double >( const Eigen::Vector3d& currentGroundStation )
{
    previousEarthFixedPosition_ = currentGroundStation;
}

//! Function to reset ground station position used on last call to updateTimes function at Time precision
template< >
void TerrestrialTimeScaleConverter::setCurrentGroundStation< Time >( const Eigen::Vector3d& currentGroundStation )
{
    previousEarthFixedPositionSplit_ = currentGroundStation;
}



//! Function to create the default Earth time scales conversion object
boost::shared_ptr< TerrestrialTimeScaleConverter > createDefaultTimeConverter( const boost::shared_ptr< EOPReader > eopReader )
{
    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator =
            getDefaultUT1CorrectionCalculator( );
    boost::shared_ptr< interpolators::JumpDataLinearInterpolator < double, double > > ut1MinusUtcInterpolator =
            boost::make_shared< interpolators::JumpDataLinearInterpolator< double, double > >(
                eopReader->getUt1MinusUtcMapInSecondsSinceJ2000( ), 0.5, 1.0 );
    return boost::make_shared< TerrestrialTimeScaleConverter >
            ( ut1MinusUtcInterpolator, shortPeriodUt1CorrectionCalculator );
}

}

}
