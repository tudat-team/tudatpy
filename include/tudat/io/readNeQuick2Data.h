/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          CCIR (1967), Atlas of Ionospheric Characteristics, Report 340-6, ITU, Geneva.
 *          S.M. Radicella (2009), The NeQuick model genesis, uses and evolution, Annals of
 *              Geophysics, 52(3/4), doi:10.4401/ag-4597
 */

#ifndef TUDAT_READNEQUICK2DATA_H
#define TUDAT_READNEQUICK2DATA_H

#include <string>
#include <array>
#include <Eigen/Core>

namespace tudat
{
namespace input_output
{

//! CCIR coefficient data for the NeQuick-2 model.
//! Each monthly file contains foF2 coefficients (13 x 76 x 2) and M(3000)F2 coefficients (9 x 49 x 2).
//! The two slices of the third dimension are for low and high solar activity respectively.
struct CcirData {
    //! foF2 coefficients: F2[month][solarActivity](row, col) where month=0..11, solarActivity=0(low)/1(high)
    //! Each matrix is 76 rows x 13 columns
    std::array< std::array< Eigen::MatrixXd, 2 >, 12 > foF2Coefficients;

    //! M(3000)F2 coefficients: FM3[month][solarActivity](row, col) where month=0..11, solarActivity=0(low)/1(high)
    //! Each matrix is 49 rows x 9 columns
    std::array< std::array< Eigen::MatrixXd, 2 >, 12 > m3000F2Coefficients;
};

//! Read all 12 monthly CCIR coefficient files from the given directory.
//! Files are named ccir11.asc (January) through ccir22.asc (December).
CcirData readCcirCoefficients( const std::string& dataPath );

//! Read the MODIP (modified dip latitude) grid file.
//! Returns a matrix with padding for interpolation: (184 x 184) stored as (row, col).
//! The grid covers latitudes -90 to +90 (1 deg step) and longitudes -180 to +180 (2 deg step).
//! Padding rows/columns are added at boundaries for 3rd-order interpolation.
Eigen::MatrixXd readModipGrid( const std::string& dataPath );

}  // namespace input_output
}  // namespace tudat

#endif  // TUDAT_READNEQUICK2DATA_H
