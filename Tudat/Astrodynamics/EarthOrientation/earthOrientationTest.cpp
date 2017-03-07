#include <iostream>

#include <boost/make_shared.hpp>

#include "Astrodynamics/EarthOrientation/eopReader.h"
#include "Astrodynamics/EarthOrientation/earthSiderealTimeCalculator.h"
#include "Astrodynamics/EarthOrientation/shortPeriodEarthOrientationCorrectionCalculator.h"
#include "Basics/utilities.h"
#include "External/SofaInterface/fundamentalArguments.h"
#include "InputOutput/writeDataToFile.h"
#include "Mathematics/Interpolators/jumpDataLinearInterpolator.h"

int main( )
{

    using namespace tudat::earth_orientation;
    using namespace tudat::sofa_interface;
    using namespace tudat;
    using namespace tudat::interpolators;

    double testTdb = 1.0E6;
    calculateDoodsonFundamentalArguments( testTdb );

    double testTai = convertTTtoTAI( testTdb );
    double testUtc = sofa_interface::convertTAItoUTC( testTai );

    std::cout<<testTdb-testTai<<" "<<testTai-testUtc<<std::endl;

    Eigen::Vector3d stationSphericalPosition;
    stationSphericalPosition << 6.378E6, mathematical_constants::PI/2.0, 0.0;
    Eigen::Vector3d stationCartesianPosition = coordinate_conversions::convertSphericalToCartesian(
                stationSphericalPosition );

    boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > ut1CorrectionCalculator = getDefaultUT1CorrectionCalculator( );

    EOPReader eopReader = EOPReader( );

    std::map< double, double > ut1OffsetMap = eopReader.getUt1MinusUtcMapInSecondsSinceJ2000( );

    boost::shared_ptr< JumpDataLinearInterpolator < double, double > > ut1OffsetInterpolator =
            boost::make_shared< JumpDataLinearInterpolator< double, double > >( ut1OffsetMap, 0.5, 1.0 );

    // Leap second was added at end of this day.
    double testMjd = 51078.0;
    double testSecondsSinceJ2000 = ( testMjd - ( basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD ) ) * 86400.0;
    double testTaiSinceJ2000 = sofa_interface::convertUTCtoTAI( testSecondsSinceJ2000 );

    std::cout<<testSecondsSinceJ2000-testTaiSinceJ2000<<std::endl;

    std::map< double, double > leapInterpolateTest_;
    std::map< double, double > shortPeriodResultMap;

    double timeStep = 10.0 * 60.0;
    int numberOfSteps = 6 * 24 * 300;
    double currentTime = testSecondsSinceJ2000;


    //output::writeMapHistoryToFile( leapInterpolateTest_, "leapSecondTest.dat" );
    //output::writeMapHistoryToFile( shortPeriodResultMap, "shortPeriodUt1.dat" );

    EarthSiderealTimeCalculator siderealTimeCalculator = EarthSiderealTimeCalculator(
                ut1OffsetInterpolator, ut1CorrectionCalculator );

    std::map< double, Eigen::VectorXd > timeMap;

    Eigen::VectorXd currentTimes;
    currentTimes.resize( 5 );

    std::map< double, Eigen::VectorXd > fundamentalArgumentDifference;

    Eigen::VectorXd preciseArguments;
    Eigen::VectorXd approximateArguments;

    std::map< double, Eigen::VectorXd > tdbTtDifferences;

    Eigen::Vector3d tdbTtValues;
    double ut1FractionOfDay;
    for( int i = 0; i < numberOfSteps; i++ )
    {
        currentTimes( 0 ) = currentTime;
        currentTimes( 1 ) = siderealTimeCalculator.getCurrentTime( tt_scale, tai_scale, currentTime, stationCartesianPosition );
        currentTimes( 2 ) = siderealTimeCalculator.getCurrentTime( tt_scale, tdb_scale, currentTime, stationCartesianPosition );
        currentTimes( 3 ) = siderealTimeCalculator.getCurrentTime( tt_scale, utc_scale, currentTime, stationCartesianPosition );
        currentTimes( 4 ) = siderealTimeCalculator.getCurrentTime( tt_scale, ut1_scale, currentTime, stationCartesianPosition );

        timeMap[ currentTime ] = currentTimes;

        preciseArguments = sofa_interface::calculateDoodsonFundamentalArguments( currentTimes( 2 ),
                                                                                 currentTime,
                                                                                 currentTimes( 4 ),
                                                                                 basic_astrodynamics::JULIAN_DAY_ON_J2000 );
        approximateArguments = sofa_interface::calculateDoodsonFundamentalArguments( currentTimes( 0 ) );
        fundamentalArgumentDifference[ currentTime ] = approximateArguments - preciseArguments;

        //leapInterpolateTest_[ currentTime ] = ut1OffsetInterpolator->interpolate( currentTime );
        //shortPeriodResultMap[ currentTime ] = ut1CorrectionCalculator->getShortPeriodUt1Correction(
        //            currentTime );
        //testTaiSinceJ2000 = sofa_interface::convertUTCtoTAI( currentTime, basic_astrodynamics::JULIAN_DAY_ON_J2000 );
        //std::cout<<i<<" "<<currentTime-testTaiSinceJ2000<<" "<<leapInterpolateTest_[ currentTime ]<<std::endl;

        ut1FractionOfDay = ( currentTimes( 4 ) -
                             std::fmod( currentTimes( 4 ), physical_constants::JULIAN_DAY ) ) / physical_constants::JULIAN_DAY;
        tdbTtValues.x( ) = sofa_interface::getTDBminusTT( currentTimes( 0 ),
                                                          ut1FractionOfDay,
                                                          0.0, 0.0, 0.0 );
        tdbTtValues.y( ) = sofa_interface::getTDBminusTT( currentTimes( 0 ),
                                                          ut1FractionOfDay,
                                                          0.0, 6378.0E3, 0.0 );
        tdbTtValues.z( ) = sofa_interface::getTDBminusTT( currentTimes( 0 ),
                                                          ut1FractionOfDay,
                                                          mathematical_constants::PI, 6378.0E3, 0.0 );

        tdbTtDifferences[ currentTime ] = tdbTtValues;

        currentTime += timeStep;
    }

    output::writeVectorHistoryToFile( timeMap, "timeMap.dat" );
    output::writeVectorHistoryToFile( fundamentalArgumentDifference, "argumentDifference.dat" );
    output::writeVectorHistoryToFile( tdbTtDifferences, "tdbTtDifferences.dat" );


}
