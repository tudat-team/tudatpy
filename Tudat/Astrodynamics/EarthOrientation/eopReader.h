#ifndef EOPREADER_H
#define EOPREADER_H

#include <map>
#include <fstream>
#include <sstream>
#include <string>
#include <iostream>
#include <iomanip>

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <boost/lexical_cast.hpp>
#include <boost/format.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Basics/utilities.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"


namespace tudat
{

struct DailyIersEopCorrections
{
    std::map< double, Eigen::Vector2d > cipInItrs; // xPole, yPole
    std::map< double, Eigen::Vector2d > cipInGcrsCorrection; // dX, dY
    std::map< double, double > ut1MinusUtc;
    std::map< double, double > lengthOfDayOffset;
};

class EOPReader
{
public:
    EOPReader( std::string eopFile = tudat::input_output::getDataFilesRootPath( ) + "/EarthOrientation/eopc04_08_IAU2000.62-now",
               std::string format = "C04",
               IAUConventions = iau_2006 );

    std::map< double, double > getUt1MinusUtcMapRaw( )
    {
        return ut1MinusUtc;
    }

    std::map< double, double > getUt1MinusUtcMapInSecondsSinceJ2000( )
    {
        return utilities::linearlyScaleKeyOfMap< double, double >
                ( ut1MinusUtc,
                  basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD,
                  physical_constants::JULIAN_DAY );

    }

    std::map< double, double > getLengthOfDayMapRaw( )
    {
        return lengthOfDayOffset;
    }

    std::map< double, double > getLengthOfDayMapInSecondsSinceJ2000( )
    {
        return utilities::linearlyScaleKeyOfMap< double, double >
                ( lengthOfDayOffset,
                  basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD,
                  physical_constants::JULIAN_DAY );

    }

    std::map< double, Eigen::Vector2d > getCipInItrsMapRaw( )
    {
        return cipInItrs;
    }

    std::map< double, Eigen::Vector2d > getCipInItrsMapInSecondsSinceJ2000( )
    {
        return utilities::linearlyScaleKeyOfMap< double, Eigen::Vector2d >
                ( cipInItrs,
                  basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD,
                  physical_constants::JULIAN_DAY );

    }

    std::map< double, Eigen::Vector2d > getCipInGcrsCorrectionMapRaw( )
    {
        return cipInGcrsCorrection;
    }

    std::map< double, Eigen::Vector2d > getCipInGcrsCorrectionMapInSecondsSinceJ2000( )
    {
        return utilities::linearlyScaleKeyOfMap< double, Eigen::Vector2d >
                ( cipInGcrsCorrection,
                  basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD,
                  physical_constants::JULIAN_DAY );

    }



private:
    void readEopFile( std::string fileName );

    std::map< double, Eigen::Vector2d > cipInItrs; // xPole, yPole

    std::map< double, Eigen::Vector2d > cipInGcrsCorrection; // dX, dY

    std::map< double, double > ut1MinusUtc;

    std::map< double, double > lengthOfDayOffset;

};

}

#endif // EOPREADER_H
