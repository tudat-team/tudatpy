#include "Tudat/Astrodynamics/EarthOrientation/terrestrialTimeScaleConverter.h"
#include "Tudat/External/SofaInterface/sofaTimeConversions.h"
#include "Tudat/Mathematics/Interpolators/jumpDataLinearInterpolator.h"

namespace tudat
{

namespace earth_orientation
{


boost::shared_ptr< TerrestrialTimeScaleConverter > createDefaultTimeConverter( boost::shared_ptr< EOPReader > eopReader )
{
    using namespace tudat::interpolators;
    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > shortPeriodUt1CorrectionCalculator =
            getDefaultUT1CorrectionCalculator( );

    boost::shared_ptr< JumpDataLinearInterpolator < double, double > > ut1MinusUtcInterpolator =
            boost::make_shared< JumpDataLinearInterpolator< double, double > >(
                eopReader->getUt1MinusUtcMapInSecondsSinceJ2000( ), 0.5, 1.0 ); // d(UT1-

    return boost::make_shared< TerrestrialTimeScaleConverter >
            ( ut1MinusUtcInterpolator, shortPeriodUt1CorrectionCalculator );
}

template< >
CurrentTimes< double >& TerrestrialTimeScaleConverter::getCurrentTimeList< double >( )
{
    return currentTimes_;
}


template< >
CurrentTimes< Time >& TerrestrialTimeScaleConverter::getCurrentTimeList< Time >( )
{
    return currentTimesSplit_;
}


}

}
