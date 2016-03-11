//#include <omp.h>

//#include "Astrodynamics/BasicAstrodynamics/timeConversions.h"
//#include "External/SofaInterface/sofaTimeConversions.h"
//#include "InputOutput/writeDataToFile.h"

//int main( )
//{
//    using namespace tudat;

//    using namespace tudat::basic_astrodynamics;
//    using namespace tudat::basic_mathematics;
//    using namespace tudat::sofa_interface;

//    boost::gregorian::date startDate = boost::gregorian::date( 2009, 02, 01 );

//    double startJulianDaysSinceJ2000 = calculateJulianDaySinceEpoch( startDate, 0.0, 0.0 );

//    std::cout<<std::setprecision(16)<<"Start time: "<<startJulianDaysSinceJ2000<<std::endl;

//    double timeStep = 1.0;

//    std::map< double, double > geocenterDifference;
//    std::map< double, double > greenwichLatitudeAndLongitudeDifference;
//    std::map< double, double > antiGreenwichLatitudeAndLongitudeDifference;
//    std::map< double, double > northPoleDifference;
//    std::map< double, double > northPoleDifference2;

//    double tdbMinusTT;
//    double approximateUt1;
//    double currentDeltaAt;

//    time_t tstart, tend;
//    tstart = time(0);

//    for( unsigned int i = 0; i < 86400; i++ )
//    {
//        double currentSecondsInDay = static_cast< double >( i ) * timeStep;

//        approximateUt1 = convertTTtoUTC( currentSecondsInDay, startJulianDaysSinceJ2000 );

//       // std::cout<<"App UT1 "<<approximateUt1<<" "<<currentSecondsInDay<<std::endl;

//        currentDeltaAt = getDeltaAtFromUtc( approximateUt1 / 86400.0, startJulianDaysSinceJ2000 );

//        tdbMinusTT = iauDtdb( startJulianDaysSinceJ2000, currentSecondsInDay / 86400.0,
//                              approximateUt1 / 86400.0, 0.0, 0.0, 0.0 );
//        geocenterDifference[ currentSecondsInDay ] = tdbMinusTT + TTMTAI + currentDeltaAt;

//        tdbMinusTT = iauDtdb( startJulianDaysSinceJ2000, currentSecondsInDay / 86400.0,
//                              approximateUt1 / 86400.0, 0.0, 6378.0, 0.0 );
//        greenwichLatitudeAndLongitudeDifference[ currentSecondsInDay ] = tdbMinusTT + TTMTAI + currentDeltaAt;
//                                                                                                //TAi-UTC

//        tdbMinusTT = iauDtdb( startJulianDaysSinceJ2000, currentSecondsInDay / 86400.0,
//                              approximateUt1 / 86400.0, mathematical_constants::PI, 6378.0, 0.0 );
//        antiGreenwichLatitudeAndLongitudeDifference[ currentSecondsInDay ] = tdbMinusTT + TTMTAI + currentDeltaAt;

//        tdbMinusTT = iauDtdb( startJulianDaysSinceJ2000, currentSecondsInDay / 86400.0,
//                              approximateUt1 / 86400.0, 0.0, 0.0, 6378.0 );
//        northPoleDifference[ currentSecondsInDay ] = tdbMinusTT + TTMTAI + currentDeltaAt;

//        tdbMinusTT = iauDtdb( startJulianDaysSinceJ2000, currentSecondsInDay / 86400.0,
//                              approximateUt1 / 86400.0, mathematical_constants::PI, 0.0, 6378.0 );
//        northPoleDifference2[ currentSecondsInDay ] = tdbMinusTT + TTMTAI + currentDeltaAt;
//    }

//    tend = time(0);
//    std::cout << "It took "<< difftime(tend, tstart) <<" second(s)."<< std::endl;

//    output::writeDoubleMapToFile( geocenterDifference, "geocenterDifference.dat" );
//    output::writeDoubleMapToFile( greenwichLatitudeAndLongitudeDifference, "greenwichEquatorDifference.dat" );
//    output::writeDoubleMapToFile( antiGreenwichLatitudeAndLongitudeDifference, "antiGreenwichEquatorDifference.dat" );
//    output::writeDoubleMapToFile( northPoleDifference, "northPoleDifference.dat" );
//    output::writeDoubleMapToFile( northPoleDifference2, "northPoleDifference2.dat" );

//}
