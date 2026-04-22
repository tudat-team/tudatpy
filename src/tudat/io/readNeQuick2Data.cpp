/*    Copyright (c) 2010-2026, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/io/readNeQuick2Data.h"
#include <fstream>
#include <sstream>
#include <iomanip>
#include <stdexcept>

namespace tudat
{
namespace input_output
{

CcirData readCcirCoefficients( const std::string& dataPath )
{
    CcirData data;

    for( int month = 0; month < 12; ++month )
    {
        // File naming: ccir11.asc (January=month 0) through ccir22.asc (December=month 11)
        std::ostringstream oss;
        oss << dataPath << "/ccir" << std::setw( 2 ) << std::setfill( '0' ) << ( month + 11 ) << ".asc";
        std::string filePath = oss.str( );

        std::ifstream file( filePath );
        if( !file )
        {
            throw std::runtime_error( "Cannot open CCIR coefficient file: " + filePath );
        }

        // The Fortran code reads: read(77,'(4e16.8)') F2, FM3
        // F2 is dimensioned (13, 76, 2) and FM3 is dimensioned (9, 49, 2)
        // Fortran column-major storage means the file contains:
        // All F2 values first (13*76*2 = 1976 values), then all FM3 values (9*49*2 = 882 values)
        // In Fortran order: innermost index varies fastest, so F2(j,i,k) reads j=1..13, i=1..76, k=1..2

        // Read all values from the file
        std::vector< double > values;
        double val;
        while( file >> val )
        {
            values.push_back( val );
        }
        file.close( );

        const std::size_t expectedF2 = 13 * 76 * 2;
        const std::size_t expectedFM3 = 9 * 49 * 2;
        if( values.size( ) < expectedF2 + expectedFM3 )
        {
            throw std::runtime_error( "CCIR file " + filePath + " has insufficient data: expected " +
                                      std::to_string( expectedF2 + expectedFM3 ) + " values, got " +
                                      std::to_string( values.size( ) ) );
        }

        // Parse F2 coefficients: F2(j, i, k) with j=1..13 (cols), i=1..76 (rows), k=1..2 (solar activity)
        for( int k = 0; k < 2; ++k )
        {
            data.foF2Coefficients[ month ][ k ] = Eigen::MatrixXd( 76, 13 );
            for( int i = 0; i < 76; ++i )
            {
                for( int j = 0; j < 13; ++j )
                {
                    std::size_t idx = j + 13 * i + 13 * 76 * k;
                    data.foF2Coefficients[ month ][ k ]( i, j ) = values[ idx ];
                }
            }
        }

        // Parse FM3 coefficients: FM3(j, i, k) with j=1..9 (cols), i=1..49 (rows), k=1..2 (solar activity)
        for( int k = 0; k < 2; ++k )
        {
            data.m3000F2Coefficients[ month ][ k ] = Eigen::MatrixXd( 49, 9 );
            for( int i = 0; i < 49; ++i )
            {
                for( int j = 0; j < 9; ++j )
                {
                    std::size_t idx = expectedF2 + j + 9 * i + 9 * 49 * k;
                    data.m3000F2Coefficients[ month ][ k ]( i, j ) = values[ idx ];
                }
            }
        }
    }

    return data;
}

Eigen::MatrixXd readModipGrid( const std::string& dataPath )
{
    std::string filePath = dataPath + "/modip.asc";
    std::ifstream file( filePath );
    if( !file )
    {
        throw std::runtime_error( "Cannot open MODIP grid file: " + filePath );
    }

    // Skip header line
    std::string headerLine;
    std::getline( file, headerLine );

    // The Fortran code uses: dimension pmodip(0:183, -1:182)
    // which is 184 rows x 184 columns
    // The grid data has 181 latitude rows (-90 to +90, 1 deg step)
    // and 181 longitude columns (-180 to +180, 2 deg step)
    // Rows are indexed as i = -90..+90 mapped to array indices 1..181 (with i+91)
    // Columns are indexed as j = -90..+90 mapped to indices 0..180 (with j+90)
    //
    // Fortran modin reads: pmodip(i+lathp+1, j+lnghp) for i=-lathp..lathp, j=-lnghp..lnghp
    // where lathp=90, lnghp=90, so pmodip(1..181, 0..180)
    // Then applies wrapping for interpolation padding.

    const int latPoints = 181;  // -90 to +90
    const int lonPoints = 181;  // -90 to +90 in half-intervals = -180 to +180 in 2-deg steps

    // Allocate with padding: indices 0..183 in lat, -1..182 in lon
    // We use a (184 x 184) matrix with offset mapping
    Eigen::MatrixXd pmodip = Eigen::MatrixXd::Zero( 184, 184 );

    // Read grid data: 181 rows (lat -90 to +90), each with 181 values (lon -180 to +180)
    for( int i = 0; i < latPoints; ++i )
    {
        for( int j = 0; j < lonPoints; ++j )
        {
            double val;
            file >> val;
            // Map to padded array: row index = i+1, col index = j+1 (offset by 1 for padding)
            pmodip( i + 1, j + 1 ) = val;
        }
    }
    file.close( );

    // Apply boundary wrapping for 3rd-order interpolation (same as Fortran modin subroutine)
    // pmodip(0, i) = pmodip(2, mod((i+lnghp), lngp))  — wrap first row
    // pmodip(latp+2, i) = pmodip(latp, mod((i+lnghp), lngp))  — wrap last+1 row
    // pmodip(latp+3, i) = pmodip(latp-1, mod((i+lnghp), lngp))  — wrap last+2 row
    const int latp = 180;
    const int lngp = 180;
    const int lnghp = 90;

    for( int i = 0; i <= lngp; ++i )
    {
        int wrappedCol = ( ( i + lnghp ) % lngp );
        // +1 offset for our 0-based storage
        pmodip( 0, i + 1 ) = pmodip( 2, wrappedCol + 1 );
        pmodip( latp + 2, i + 1 ) = pmodip( latp, wrappedCol + 1 );
        pmodip( latp + 3, i + 1 ) = pmodip( latp - 1, wrappedCol + 1 );
    }

    // Longitude wrapping: pmodip(i, -1) = pmodip(i, lngp-1), pmodip(i, lngp+1) = pmodip(i, 1)
    for( int i = 0; i < 184; ++i )
    {
        pmodip( i, 0 ) = pmodip( i, lngp );   // col -1 -> col lngp-1 (with +1 offset: 0 -> lngp)
        pmodip( i, lngp + 2 ) = pmodip( i, 2 );  // col lngp+1 -> col 1 (with +1 offset: lngp+2 -> 2)
    }

    return pmodip;
}

}  // namespace input_output
}  // namespace tudat
