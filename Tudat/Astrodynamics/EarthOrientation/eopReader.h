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

#ifndef TUDAT_EOPREADER_H
#define TUDAT_EOPREADER_H

#include <map>
#include <string>
#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <Eigen/Core>

#include "Tudat/Astrodynamics/BasicAstrodynamics/physicalConstants.h"
#include "Tudat/Mathematics/BasicMathematics/mathematicalConstants.h"

#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Basics/utilities.h"
#include "Tudat/External/SofaInterface/earthOrientation.h"


namespace tudat
{

namespace earth_orientation
{

//! Class used to read Earth Orientation Parameters (EOP) from file
class EOPReader
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param eopFile Name of EOP file that is to be used
     * \param format Identifier for file format that is provied
     * \param nutationTheory Nutation theory w.r.t. which the EOP data is given.
     */
    EOPReader(
            const std::string& eopFile = tudat::input_output::getEarthOrientationDataFilesPath( ) + "eopc04_08_IAU2000.62-now.txt",
            const std::string& format = "C04",
            const basic_astrodynamics::IAUConventions nutationTheory = basic_astrodynamics::iau_2006 );

    //! Function to retrieve the data of UT1-UTC, as provided in the EOP file.
    std::map< double, double > getUt1MinusUtcMapRaw( )
    {
        return ut1MinusUtc;
    }

    //! Function to retrieve the data of UT1-UTC, with map key seconds since J2000
    std::map< double, double > getUt1MinusUtcMapInSecondsSinceJ2000( )
    {
        return utilities::linearlyScaleKeyOfMap< double, double >
                ( ut1MinusUtc,
                  basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD,
                  physical_constants::JULIAN_DAY );

    }

    //! Function to retrieve the data of LOD offset, as provided in the EOP file.
    /*!
     * Function to retrieve the data of LOD offset, as provided in the EOP file.
     * \return Data of LOD offset, as provided in the EOP file.
     */
    std::map< double, double > getLengthOfDayMapRaw( )
    {
        return lengthOfDayOffset;
    }

    //! Function to retrieve the data of LOD offset, with map key seconds since J2000
    /*!
     * Function to retrieve the data of LOD offset, with map key seconds since J2000
     * \return Data of LOD offset, with map key seconds since J2000
     */
    std::map< double, double > getLengthOfDayMapInSecondsSinceJ2000( )
    {
        return utilities::linearlyScaleKeyOfMap< double, double >
                ( lengthOfDayOffset,
                  basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD,
                  physical_constants::JULIAN_DAY );

    }

    //! Function to retrieve the data of CIP in ITRS correction (polar motion), as provided in the EOP file.
    /*!
    *  Function to retrieve the data of CIP in ITRS correction (polar motion), as provided in the EOP file.
    *  \return Data of CIP in ITRS correction (polar motion), as provided in the EOP file.
    */
    std::map< double, Eigen::Vector2d > getCipInItrsMapRaw( )
    {
        return cipInItrs;
    }

    //! Function to retrieve the data of CIP in ITRS correction (polar motion), with map key seconds since J2000
    /*!
    *  Function to retrieve the data of CIP in ITRS correction (polar motion), with map key seconds since J2000
    *  \return Data of CIP in ITRS correction (polar motion), with map key seconds since J2000
    */
    std::map< double, Eigen::Vector2d > getCipInItrsMapInSecondsSinceJ2000( )
    {
        return utilities::linearlyScaleKeyOfMap< double, Eigen::Vector2d >
                ( cipInItrs,
                  basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD,
                  physical_constants::JULIAN_DAY );

    }

    //! Function to retrieve the data of CIP in GCRS correction (nutation), as provided in the EOP file.
    /*!
    *  Function to retrieve the data of CIP in GCRS correction (nutation), as provided in the EOP file.
    *  \return Data of CIP in GCRS correction (nutation), as provided in the EOP file.
    */
    std::map< double, Eigen::Vector2d > getCipInGcrsCorrectionMapRaw( )
    {
        return cipInGcrsCorrection;
    }

    //! Function to retrieve the data of CIP in GCRS correction (nutation), with map key seconds since J2000
    /*!
    *  Function to retrieve the data of CIP in GCRS correction (nutation), with map key seconds since J2000
    *  \return Data of CIP in GCRS correction (nutation), with map key seconds since J2000
    */
    std::map< double, Eigen::Vector2d > getCipInGcrsCorrectionMapInSecondsSinceJ2000( )
    {
        return utilities::linearlyScaleKeyOfMap< double, Eigen::Vector2d >
                ( cipInGcrsCorrection,
                  basic_astrodynamics::JULIAN_DAY_ON_J2000 - basic_astrodynamics::JULIAN_DAY_AT_0_MJD,
                  physical_constants::JULIAN_DAY );

    }

private:

    //! Function to read EOP file
    /*!
     * Function to read EOP file
     * \param fileName EOP file name.
     */
    void readEopFile( const std::string& fileName );

    //! Terrestrial pole position corrections (CIP in ITRS; polar motion), read from file
    std::map< double, Eigen::Vector2d > cipInItrs;

    //! Celestial pole position corrections (CIP in GCRS; nutation), read from file
    std::map< double, Eigen::Vector2d > cipInGcrsCorrection;

    //! Corrections to UT1 - UTC, read from file
    std::map< double, double > ut1MinusUtc;

    //! Corrections to LOD, read from file
    std::map< double, double > lengthOfDayOffset;

};

}

}

#endif // TUDAT_EOPREADER_H
