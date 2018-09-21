/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_SPHERICALHARMONICCOSINECOEFFICIENTS_H
#define TUDAT_SPHERICALHARMONICCOSINECOEFFICIENTS_H

#include <map>

#include <functional>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of spherical harmonic gravity field cosine coefficients
/*!
 * Interface class for estimation of spherical harmonic gravity field cosine coefficients. Interfaces the estimation with the
 * member coefficients of an instance of the SphericalHarmonicsGravityField (or nominal coefficients of
 * TimeDependentSphericalHarmonicsGravityField).
 */
class SphericalHarmonicsCosineCoefficients: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param getCosineCoefficients Function to retrieve the full set of sine coefficients, of which a subset is to be
     * estimated.
     * \param setCosineCoefficients Function to reset the full set of sine coefficients, of which a subset is to be
     * estimated.
     * \param blockIndices List of cosine coefficient indices which are to be estimated (first and second
     * are degree and order for each vector entry).
     * \param associatedBody Name of body for which cosine coefficients are to be estimated.
     */
    SphericalHarmonicsCosineCoefficients(
            const std::function< Eigen::MatrixXd( ) > getCosineCoefficients,
            const std::function< void( Eigen::MatrixXd ) > setCosineCoefficients,
            const std::vector< std::pair< int, int > > blockIndices,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( spherical_harmonics_cosine_coefficient_block, associatedBody ),
        getCosineCoefficients_( getCosineCoefficients ), setCosineCoefficients_( setCosineCoefficients ),
        blockIndices_( blockIndices )
    {
        parameterSize_ = blockIndices_.size( );
    }

    //! Destructor
    ~SphericalHarmonicsCosineCoefficients( ) { }

    //! Function to retrieve the current values of the cosine coefficients that are to be estimated.
    /*!
     * Function to retrieve the current values of the cosine coefficients that are to be estimated. The order of the entries
     * in this vector is the same as in the blockIndices_ member vector.
     * \return List of values for cosine coefficients that are to be estimated.
     */
    Eigen::VectorXd getParameterValue( );

    //! Function to reset the cosine coefficients that are to be estimated.
    /*!
     * Function to reset the cosine coefficients that are to be estimated. The order of the entries in this vector
     * is the same as in the blockIndices_ member vector.
     * \param parameterValue List of new values for cosine coefficients that are to be estimated.
     */
    void setParameterValue( const Eigen::VectorXd parameterValue );

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter, equal to the length of the blcokIndices_ member vector, i.e. the
     *  number of coefficients that are to be estimated.
     *  \return Size of parameter value.
     */
    int getParameterSize( ) { return parameterSize_; }

    //! Function to retrieve the list of cosine coefficient indices which are to be estimated.
    /*!
     * Function to retrieve the list of cosine coefficient indices which are to be estimated (first and second are degree
     *  and order for each vector entry).
     * \return List of cosine coefficient indices which are to be estimated.
     */
    std::vector< std::pair< int, int > > getBlockIndices( ){ return blockIndices_; }


    std::string getParameterDescription( )
    {
        std::string parameterDescription =
                getParameterTypeString( parameterName_.first ) + "of (" + parameterName_.second.first + "), ";
        parameterDescription += "Minimum D/O: (" +
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
     * \param c20Index Index for degree=2, order=0 entry (-1 if none; returned by reference)
     * \param c21Index Index for degree=2, order=1 entry (-1 if none; returned by reference)
     * \param c22Index Index for degree=2, order=2 entry (-1 if none; returned by reference)
     */
    void getDegreeTwoEntries(
            int& c20Index, int& c21Index, int& c22Index )
    {
        c20Index = -1;
        c21Index = -1;
        c22Index = -1;

        for( unsigned int i = 0; i < blockIndices_.size( ); i++ )
        {
            if( blockIndices_.at( i ).first == 2 && blockIndices_.at( i ).second == 0 )
            {
               c20Index = i;
            }

            if( blockIndices_.at( i ).first == 2 && blockIndices_.at( i ).second == 1 )
            {
               c21Index = i;
            }

            if( blockIndices_.at( i ).first == 2 && blockIndices_.at( i ).second == 2 )
            {
               c22Index = i;
            }
        }
    }
protected:

private:


    //! Function to retrieve the full set of sine coefficients, of which a subset is to be estimated.
    std::function< Eigen::MatrixXd( ) > getCosineCoefficients_;

    //! Function to reset the full set of sine coefficients, of which a subset is to be estimated.
    std::function< void( Eigen::MatrixXd ) > setCosineCoefficients_;

    //! List of cosine coefficient indices which are to be estimated
    /*!
      *  List of cosine coefficient indices which are to be estimated (first and second are degree and order for
      *  each vector entry).
      */
    std::vector< std::pair< int, int > > blockIndices_;

    int parameterSize_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_SPHERICALHARMONICCOSINECOEFFICIENTS_H
