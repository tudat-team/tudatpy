/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Gentry, A., Smyth, D., and Oliver, W. The Mark IV Supersonic-Hypersonic
 *        Arbitrary Body Program, Volume II - Program Formulation, Douglas
 *        Aircraft Company, 1973.
 *
 */

#ifndef TUDAT_AERODYNAMIC_COEFFICIENT_GENERATOR_H
#define TUDAT_AERODYNAMIC_COEFFICIENT_GENERATOR_H

#include <vector>

#include <boost/multi_array.hpp>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/io/multiDimensionalArrayWriter.h"
#include "tudat/math/interpolators/multiLinearInterpolator.h"
#include "tudat/basics/basicTypedefs.h"
#include "tudat/basics/utilities.h"

namespace tudat
{

namespace aerodynamics
{

//! Function to print to console which aerodynamic coefficients are being saved.
/*!
 *  Function to print to console which aerodynamic coefficients are being saved.
 *  \param coefficientIndices Indices of coefficients to be saved.
 */
void informUserOnSavedCoefficient( std::vector< unsigned int > coefficientIndices );

//! Base class for aerodynamic coefficient generator.
/*!
 * Abstract base class for aerodynamic analysis method. Stores independent variable values
 * and data points of independent variables. Coefficients are stored in a multi_array of pointers.
 */
template< unsigned int NumberOfIndependentVariables, unsigned int NumberOfCoefficients = 6 >
class AerodynamicCoefficientGenerator: public AerodynamicCoefficientInterface
{
public:

    //! Default base class constructor.
    /*!
     *  Default base class constructor, sets independent variable data points and aerodynamics
     *  reference quantities.
     *  \param dataPointsOfIndependentVariables Vector of vector, where each subvector contains
     *  the data points of each of the independent variables for the coefficient generation.
     *  The number of subvectors must be equal to the number of independent variables. It is
     *  recommended that each of the subvectors is sorted in ascending order.
     *  \param referenceLength Reference length with which aerodynamic moments (about x- and
     *  z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
     *  is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz) (default true).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction (default true).
     */
    AerodynamicCoefficientGenerator(
            const std::vector< std::vector< double > >& dataPointsOfIndependentVariables,
            const double referenceLength,
            const double referenceArea,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames,
            const AerodynamicCoefficientFrames forceCoefficientsFrame = negative_aerodynamic_frame_coefficients,
            const AerodynamicCoefficientFrames momentCoefficientsFrame = body_fixed_frame_coefficients ) :
        AerodynamicCoefficientInterface(
            referenceLength, referenceArea, momentReferencePoint,
            independentVariableNames, forceCoefficientsFrame,
            momentCoefficientsFrame ),
        dataPointsOfIndependentVariables_( dataPointsOfIndependentVariables )
    {
        // Check that the size of dataPointsOfIndependentVariables matches the template parameter.
        if( !( dataPointsOfIndependentVariables_.size( ) == NumberOfIndependentVariables ) )
        {
            throw std::runtime_error( "Error in AerodynamicCoefficientGenerator, input data is inconsistent" );
        }

        boost::array< int, NumberOfIndependentVariables > numberOfPointsPerIndependentVariables;
        for( unsigned int i = 0; i < NumberOfIndependentVariables; i++ )
        {
            numberOfPointsPerIndependentVariables[ i ] = dataPointsOfIndependentVariables_[ i ].
                    size( );
        }

        aerodynamicCoefficients_.resize( numberOfPointsPerIndependentVariables );
    }

    //! Default destructor.
    virtual ~AerodynamicCoefficientGenerator( ) { }

    //! Get the number of points for an independent variable.
    /*!
     * Returns the number of different values of a given independent variable
     * at which coefficients are determined.
     * \param independentVariable Independent variable from which to retrieve
     * number of data points.
     * \return Number of data points for selected independent variable.
     */
    int getNumberOfValuesOfIndependentVariable( const int independentVariable ) const
    {
        return dataPointsOfIndependentVariables_ [ independentVariable ].size( );
    }

    //! Get a value of an independent variable.
    /*!
     * Returns a value of an independent variable at which aerodynamic
     * coefficients are to be determined.
     * \param independentVariable Independent variable index.
     * \param index Index from independent variable list for which to retrieve the value.
     * \return Value of independent variable at index.
     */
    double getIndependentVariablePoint( const int independentVariable, const int index ) const
    {
        return dataPointsOfIndependentVariables_[ independentVariable ][ index ];
    }

    //! Get aerodynamic coefficients.
    /*!
     * Virtual function to return the aerodynamic coefficients at a specified
     * set of independent variables.
     * \param independentVariables Array of values of independent variable indices in
     * dataPointsOfIndependentVariables_ for which to retrieve aerodynamic coefficients.
     * \return vector of coefficients at specified independent variable indices.
     */
    virtual Eigen::Matrix< double, NumberOfCoefficients, 1 > getAerodynamicCoefficientsDataPoint(
            const boost::array< int, NumberOfIndependentVariables > independentVariables ) = 0;

    //! Function to return the complete set of aerodynamic coefficients that have been calculated.
    /*!
     *  Function to return the complete set of aerodynamic coefficients that have been calculated.
     *  \return Complete set of aerodynamic coefficients that have been calculated.
     */
    boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients, 1 >,
    NumberOfIndependentVariables > getAerodynamicCoefficientsTables( )
    {
        return aerodynamicCoefficients_;
    }

    //! Get the data points of the independent variables at which the coefficients are calculated.
    /*!
     *  Get the data points of the independent variables at which the coefficients are calculated.
     *  The aerodynamic coefficients are calculated each set of combinations of the independent
     *  variables.
     *  \return Data points of the independent variables at which the coefficients are calculated.
     */
    std::vector< std::vector< double > > getDataPointsOfIndependentVariables( )
    {
        return dataPointsOfIndependentVariables_;
    }

    //! Save aerodynamic coefficients to file.
    /*!
     *  Save aerodynamic coefficients to file.
     *  \param fileNamesMap Map of paths to files where aerodynamics coefficients are to be saved.
     */
    void saveAerodynamicCoefficientsTables( const std::map< unsigned int, std::string >& fileNamesMap )
    {
        // Inform user on which variable is being saved
        informUserOnSavedCoefficient( utilities::createVectorFromMapKeys( fileNamesMap ) );

        // Write coefficients to file
        input_output::MultiArrayFileWriter< NumberOfIndependentVariables,
                NumberOfCoefficients >::writeMultiArrayAndIndependentVariablesToFiles( fileNamesMap,
                                                                                       dataPointsOfIndependentVariables_,
                                                                                       aerodynamicCoefficients_ );
    }

    //! Save aerodynamic coefficients to a single file.
    /*!
     *  Save aerodynamic coefficients to a single file, in case of one independent variable is used.
     *  \param fileName Paths to file where aerodynamics coefficients are to be saved.
     *  \param coefficientIndices Indices of coefficients to be saved. Default value is all of them.
     */
    void saveAerodynamicCoefficientsTables( const std::string& fileName,
                                            const std::vector< unsigned int > coefficientIndices = { 0, 1, 2, 3, 4, 5 } )
    {
        if ( NumberOfIndependentVariables == 1 )
        {
            // Inform user on which variable is being saved
            informUserOnSavedCoefficient( coefficientIndices );

            // Write coefficients to file
            input_output::MultiArrayFileWriter< 1, NumberOfCoefficients >::
                    writeMultiArrayAndIndependentVariablesToFiles( fileName, coefficientIndices,
                                                                   dataPointsOfIndependentVariables_,
                                                                   aerodynamicCoefficients_ );
        }
        else
        {
            throw std::runtime_error( "Error in aerodynamic coefficient generator. The saveAerodynamicCoefficientsTables with "
                                      "single file path can only be used in case only one independent variable is used. "
                                      "Number of independent variables: " + std::to_string( NumberOfIndependentVariables ) );
        }
    }

    //! Compute the aerodynamic coefficients at current flight condition.
    /*!
     *  Compute the aerodynamic coefficients at current flight conditions (independent variables).
     *  Input is a set of independent variables (doubles) which represent the variables from which
     *  the coefficients are calculated. The physical nature of these variables depends on
     *  the coefficientFunction_ variables. The size of the independent variable vector must be
     *  numberOfIndependentVariables_
     *  \param independentVariables Independent variables of force and moment coefficient
     *  determination implemented by derived class
     *  \param currentTime Time to which coefficients are to be updated (not used in this derived class).
     */
    virtual void updateCurrentCoefficients( const std::vector< double >& independentVariables,
                                            const double currentTime = TUDAT_NAN )
    {
        // Check if the correct number of aerodynamic coefficients is provided.
        if( independentVariables.size( ) != numberOfIndependentVariables_ )
        {
            std::string errorMessage =
                    "Error in AerodynamicCoefficientGenerator, number of input variables is inconsistent " +
                    std::to_string( independentVariables.size( ) ) + ", " +
                    std::to_string( numberOfIndependentVariables_ );
            throw std::runtime_error( errorMessage );
        }

        // Update current coefficients.
        Eigen::Vector6d currentCoefficients = coefficientInterpolator_->interpolate(
                    independentVariables );


        currentForceCoefficients_ = currentCoefficients.segment( 0, 3 );
        currentMomentCoefficients_ = currentCoefficients.segment( 3, 3 );

    }

    void clearBaseData( )
    {
        boost::array< int, NumberOfCoefficients > numberOfPointsPerIndependentVariables;
        for( unsigned int i = 0; i < NumberOfCoefficients; i++ )
        {
            numberOfPointsPerIndependentVariables[ i ] = 0;
        }
        aerodynamicCoefficients_.resize( numberOfPointsPerIndependentVariables );

        for( unsigned int i = 0; i < dataPointsOfIndependentVariables_.size( ); i++ )
        {
            dataPointsOfIndependentVariables_.at( i ).clear( );
        }
        dataPointsOfIndependentVariables_.clear( );

    }

protected:

    //! Generate aerodynamic coefficients.
    /*!
     * Virtual function to generate aerodynamic coefficients for the list of independent variables
     * passed to the constructor.
     */
    virtual void generateCoefficients( ) = 0;

    //! Function to create the coefficient interpolator from the discrete set in
    //! aerodynamicCoefficients_
    void createInterpolator( )
    {
        // Create interpolator for coefficients.
        coefficientInterpolator_ =
                std::make_shared< interpolators::MultiLinearInterpolator< double,
                Eigen::Vector6d, NumberOfIndependentVariables > >
                ( dataPointsOfIndependentVariables_, aerodynamicCoefficients_ );
    }

    //! N-dimensional array containing all computed aerodynamic coefficients.
    /*!
     *  N-dimensional array containing all computer aerodynamic coefficients. The k-th dimension
     *  pertains to coefficients at the k-th independent variable, the data points for which are
     *  defined by dataPointsOfIndependentVariables_ and the physical meaning of which are defined
     *  by independentVariableNames_
     */
    boost::multi_array< Eigen::Matrix< double, NumberOfCoefficients, 1 >,
    NumberOfIndependentVariables > aerodynamicCoefficients_;

    //! Data points of the independent variables at which the coefficients are calculated.
    /*!
     *  Data points of the independent variables at which the coefficients are calculated. The
     *  k-th vector contains the vector of data points to which the k-th independent variables
     *  (defined by independentVariableNames_) are set during the calculation of the aerodynamic
     *  coefficients.
     */
    std::vector< std::vector< double > > dataPointsOfIndependentVariables_;

    //! Interpolator producing continuous aerodynamic coefficients from the discrete calculations
    //! contained in aerodynamicCoefficients_.
    std::shared_ptr< interpolators::Interpolator< double, Eigen::Vector6d > >
    coefficientInterpolator_;
};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_AERODYNAMIC_COEFFICIENT_GENERATOR_H
