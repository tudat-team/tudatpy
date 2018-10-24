/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALHARMONICSINECOEFFICIENTS_H
#define TUDAT_SPHERICALHARMONICSINECOEFFICIENTS_H

#include <map>

#include <functional>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of spherical harmonic gravity field sine coefficients
/*!
 * Interface class for estimation of spherical harmonic gravity field sine coefficients. Interfaces the estimation with the
 * member coefficients of an instance of the SphericalHarmonicsGravityField (or nominal coefficients of
 * TimeDependentSphericalHarmonicsGravityField).
 */
class SphericalHarmonicsSineCoefficients: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param getSineCoefficients Function to retrieve the full set of sine coefficients, of which a subset is to be
     * estimated.
     * \param setSineCoefficients Function to reset the full set of sine coefficients, of which a subset is to be
     * estimated.
     * \param blockIndices List of sine coefficient indices which are to be estimated (first and second
     * are degree and order for each vector entry).
     * \param associatedBody Name of body for which sine coefficients are to be estimated.
     */
    SphericalHarmonicsSineCoefficients(
            const std::function< Eigen::MatrixXd( ) > getSineCoefficients,
            const std::function< void( Eigen::MatrixXd ) > setSineCoefficients,
            const std::vector< std::pair< int, int > >& blockIndices,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( spherical_harmonics_sine_coefficient_block, associatedBody ),
        getSineCoefficients_( getSineCoefficients ), setSineCoefficients_( setSineCoefficients ),
        blockIndices_( blockIndices )
    {
        parameterSize_ = blockIndices_.size( );
    }

    //! Destructor
    ~SphericalHarmonicsSineCoefficients( ) { }

    //! Function to retrieve the current values of the sine coefficients that are to be estimated.
    /*!
     * Function to retrieve the current values of the sine coefficients that are to be estimated. The order of the entries
     * in this vector is the same as in the blockIndices_ member vector.
     * \return List of values for sine coefficients that are to be estimated.
     */
    Eigen::VectorXd getParameterValue( );

    //! Function to reset the sine coefficients that are to be estimated.
    /*!
     * Function to reset the sine coefficients that are to be estimated. The order of the entries in this vector
     * is the same as in the blockIndices_ member vector.
     * \param parameterValue List of new values for sine coefficients that are to be estimated.
     */
    void setParameterValue( const Eigen::VectorXd parameterValue );

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter, equal to the length of the blcokIndices_ member vector, i.e. the
     *  number of coefficients that are to be estimated.
     *  \return Size of parameter value.
     */
    int getParameterSize( ) { return parameterSize_; }

    //! Function to retrieve the list of sine coefficient indices which are to be estimated.
    /*!
     * Function to retrieve the list of sine coefficient indices which are to be estimated (first and second are degree
     *  and order for each vector entry).
     * \return List of sine coefficient indices which are to be estimated.
     */
    std::vector< std::pair< int, int > > getBlockIndices( )
    {
        return blockIndices_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription =
                getParameterTypeString( parameterName_.first ) + "of (" + parameterName_.second.first + ")";
        parameterDescription += ", Minimum D/O: (" +
                std::to_string( blockIndices_.at( 0 ).first ) + ", " +
                std::to_string( blockIndices_.at( 0 ).second ) + "), ";

        parameterDescription += "Maximum D/O: (" +
                std::to_string( blockIndices_.at( blockIndices_.size( ) - 1 ).first ) + ", " +
                std::to_string( blockIndices_.at( blockIndices_.size( ) - 1 ).second ) + "). ";
        return parameterDescription;
    }

    //! Function that returns the indices for degree two coefficients (if any)
    /*!
     * Function that returns the indices for degree two coefficients (if any)
     * \param s21Index Index for degree=2, order=1 entry (-1 if none; returned by reference)
     * \param s22Index Index for degree=2, order=2 entry (-1 if none; returned by reference)
     */
    void getDegreeTwoEntries(
            int& s21Index, int& s22Index )
    {
        s21Index = -1;
        s22Index = -1;

        for( unsigned int i = 0; i < blockIndices_.size( ); i++ )
        {

            if( blockIndices_.at( i ).first == 2 && blockIndices_.at( i ).second == 1 )
            {
               s21Index = i;
            }

            if( blockIndices_.at( i ).first == 2 && blockIndices_.at( i ).second == 2 )
            {
               s22Index = i;
            }
        }
    }

protected:

private:

    //! Function to retrieve the full set of sine coefficients, of which a subset is to be estimated.
    std::function< Eigen::MatrixXd( ) > getSineCoefficients_;

    //! Function to reset the full set of sine coefficients, of which a subset is to be estimated.
    std::function< void( Eigen::MatrixXd ) > setSineCoefficients_;

    //! List of sine coefficient indices which are to be estimated
    /*!
      *  List of sine coefficient indices which are to be estimated (first and second are degree and order for
      *  each vector entry).
      */
    std::vector< std::pair< int, int > > blockIndices_;

    //! Number of coefficients that are to be estimated (i.e. length of blockIndices_ vector).
    int parameterSize_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_SPHERICALHARMONICSINECOEFFICIENTS_H
