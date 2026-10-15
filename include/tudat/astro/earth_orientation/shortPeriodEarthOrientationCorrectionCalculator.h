/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_SHORTPERIODEARTHORIENTATIONCORRECTIONCALCULATOR_H
#define TUDAT_SHORTPERIODEARTHORIENTATIONCORRECTIONCALCULATOR_H

#include <cmath>
#include <stdexcept>
#include <string>

#include <functional>

#include <Eigen/Core>

#include "tudat/basics/basicTypedefs.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/astro/basic_astro/unitConversions.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/earth_orientation/readAmplitudeAndArgumentMultipliers.h"

#include "tudat/interface/sofa/fundamentalArguments.h"
#include "tudat/io/basicInputOutput.h"
#include "tudat/math/basic/numericalDerivative.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/basics/timeType.h"

namespace tudat
{

namespace earth_orientation
{

//! Model for the diurnal and semidiurnal ocean-tide variations of polar motion and UT1. The conventional
//! libration terms are included for either model.
enum class OceanTideEopModel { iers_2010, desai_sibois_2016 };

//! Object to calculate the short period variations in Earth orientaion parameters
/*!
 *  Object to calculate the short period  variations in Earth orientaion parameters, e.g. taking into account
 *  variations due to both libration and ocean tides.
 */
template< typename OutputType >
class ShortPeriodEarthOrientationCorrectionCalculator
{
public:
    //! Constructor, requires files defining the amplitudes and fundamental argument multipliers of the variations.
    /*!
     *  Constructor, requires files defining the sine and cosine amplitudes and fundamental argument multipliers numbers
     *  of the variations for both ocean tide and libration variations. A cut-off amplitude (taken as
     *  the RSS of the amplitudes) below which amplitudes in the files are ignored can be provided.
     *  \param conversionFactor Conversion factor to be used for amplitudes, used to multiply input values, typically for unit
     *  conversion purposes.
     *  \param minimumAmplitude Minimum amplitude that is read from files and considered in calculations.
     *  \param amplitudesFiles List of files with amplitudes for corrections
     *  \param argumentMultipliersFile Fundamental argument multiplier for corrections
     *  \param argumentFunction Fundamental argument functions associated with multipliers, default Delaunay arguments with GMST
     */
    ShortPeriodEarthOrientationCorrectionCalculator(
            const double conversionFactor,
            const double minimumAmplitude,
            const std::vector< std::string >& amplitudesFiles,
            const std::vector< std::string >& argumentMultipliersFile,
            const std::function< Eigen::Vector6d( const double ) > argumentFunction =
                    std::bind( &sofa_interface::calculateApproximateDelaunayFundamentalArgumentsWithGmst, std::placeholders::_1 ),
            const std::shared_ptr< interpolators::InterpolatorGenerationSettings< double > > shortTermInterpolatorSettings = nullptr ):
        argumentFunction_( argumentFunction )
    {
        if( amplitudesFiles.size( ) != argumentMultipliersFile.size( ) )
        {
            throw std::runtime_error( "Error when calling ShortPeriodEarthOrientationCorrectionCalculator, input size is inconsistent" );
        }

        // Read data from files
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > dataFromFile;
        for( unsigned int i = 0; i < amplitudesFiles.size( ); i++ )
        {
            dataFromFile = readAmplitudesAndFundamentalArgumentMultipliers(
                    amplitudesFiles.at( i ), argumentMultipliersFile.at( i ), minimumAmplitude );
            argumentAmplitudes_.push_back( conversionFactor * dataFromFile.first );
            argumentMultipliers_.push_back( dataFromFile.second );
        }

        if( shortTermInterpolatorSettings != nullptr )
        {
            std::function< OutputType( const double ) > correctionFunction = std::bind(
                    &ShortPeriodEarthOrientationCorrectionCalculator< OutputType >::getCorrections, this, std::placeholders::_1 );
            correctionInterpolator_ = interpolators::createOneDimensionalInterpolator< double, OutputType >(
                    correctionFunction, shortTermInterpolatorSettings );
        }
        else
        {
            correctionInterpolator_ = nullptr;
        }
    }

    //! Function to obtain short period corrections.
    /*!
     *  Function to obtain short period corrections, using time as input. Fundamental arguments are calculated internally.
     *  \param ephemerisTime Time (TDB seconds since J2000) at which corretions are to be determined
     *  \return Short period corrections
     */
    OutputType getCorrections( const double& ephemerisTime )
    {
        if( correctionInterpolator_ == nullptr )
        {
            return sumCorrectionTerms( argumentFunction_( ephemerisTime ) );
        }
        else
        {
            try
            {
                return correctionInterpolator_->interpolate( ephemerisTime );
            }
            catch( std::runtime_error& caughtException )
            {
                throw std::runtime_error( "Error in short period EOP calculator.\nOriginal error: " +
                                          std::string( caughtException.what( ) ) );
            }
        }
    }

    //! Function to obtain short period corrections.
    /*!
     *  Function to obtain short period corrections, using fundamental arguments as input.
     *  \param fundamentalArguments Fundamental arguments from which corretions are to be determined
     *  \return Short period corrections
     */
    OutputType getCorrectionsFromFundamentalArgument( const Eigen::Vector6d& fundamentalArguments )
    {
        return sumCorrectionTerms( fundamentalArguments );
    }

    //! Function to obtain the time derivative of the short period corrections.
    /*!
     *  Function to obtain the time derivative of the short period corrections, per TT second. The corrections are
     *  differentiated numerically, with a fourth-order central difference in TT. At each of the nearby TT epochs, UT1 is
     *  shifted from the given UT1 by dUt1DTt times the TT shift, so that the daily Earth orientation data does not need to be
     *  evaluated at these epochs.
     *  \param tt TT (seconds since J2000) at which the rate is to be determined
     *  \param ut1 UT1 (seconds since J2000) at the same instant
     *  \param dUt1DTt Derivative of UT1 with respect to TT at this instant
     *  \param stepSeconds Step size (in seconds of TT) of the central difference
     *  eturn Time derivative of the short period corrections (radians per second for polar motion, seconds per second for UT1)
     */
    OutputType getCorrectionRate( const Time& tt, const Time& ut1, const double dUt1DTt, const double stepSeconds = 32.0 )
    {
        if( !std::isfinite( stepSeconds ) || stepSeconds <= 0.0 || !std::isfinite( dUt1DTt ) )
        {
            throw std::invalid_argument( "Invalid input for the rate of the short period Earth orientation corrections" );
        }
        const std::function< OutputType( const Time& ) > correctionFunction =
                std::bind( &ShortPeriodEarthOrientationCorrectionCalculator< OutputType >::getCorrectionsAtShiftedEpoch,
                           this,
                           tt,
                           ut1,
                           dUt1DTt,
                           std::placeholders::_1 );
        return numerical_derivatives::computeFourthOrderCentralDifference( correctionFunction, tt, stepSeconds );
    }

private:
    //! Function to obtain the short period corrections at a TT epoch close to a reference epoch.
    /*!
     *  Function to obtain the short period corrections at a TT epoch close to a reference epoch, with UT1 shifted from its
     *  value at the reference epoch by dUt1DTt times the difference between the two TT epochs.
     *  \param tt Reference TT (seconds since J2000)
     *  \param ut1 UT1 (seconds since J2000) at the reference epoch
     *  \param dUt1DTt Derivative of UT1 with respect to TT at the reference epoch
     *  \param shiftedTt TT (seconds since J2000) at which the corrections are to be determined
     *  eturn Short period corrections
     */
    OutputType getCorrectionsAtShiftedEpoch( const Time& tt, const Time& ut1, const double dUt1DTt, const Time& shiftedTt )
    {
        const Time shiftedUt1 = ut1 + static_cast< double >( shiftedTt - tt ) * dUt1DTt;
        const Eigen::Vector6d arguments = sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst(
                static_cast< double >( shiftedTt ), static_cast< double >( shiftedTt ), static_cast< double >( shiftedUt1 ) );
        return sumCorrectionTerms( arguments );
    }

    //! Function to sum all the corrcetion terms.
    /*!
     *  Function to sum all the corrcetion terms.
     * \param arguments Values of fundamental arguments
     * \return Total correction at current fundamental arguments
     */
    OutputType sumCorrectionTerms( const Eigen::Vector6d& arguments );

    //! Amplitudes of libration-induced variations
    std::vector< Eigen::MatrixXd > argumentAmplitudes_;

    //! Fundamental argument multipliers of libration-induced variations.
    std::vector< Eigen::MatrixXd > argumentMultipliers_;

    //! Fundamental argument functions associated with multipliers.
    std::function< Eigen::Vector6d( const double ) > argumentFunction_;

    std::shared_ptr< interpolators::OneDimensionalInterpolator< double, OutputType > > correctionInterpolator_;
};

//! Function to retrieve the default UT1 short-period correction calculator
/*!
 * Function to retrieve the default UT1 short-period correction calculator, from Tables  5.1, 8.2 and 8.2. of IERS 2010
 * Conventions. An amplitude cutoff may be provided.
 * \param minimumAmplitude Variation amplitude below which corrections are not taken into account.
 * \return Default UT1 short-period correction calculator
 */
std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > getDefaultUT1CorrectionCalculator(
        const double minimumAmplitude = 0.0 );

//! Function to retrieve the default polar motion short-period correction calculator
/*!
 * Function to retrieve the default polar motion short-period correction calculator, from Tables  5.1, 8.2 and 8.2. of IERS 2010
 * Conventions. An amplitude cutoff may be provided.
 * \param minimumAmplitude Variation amplitude below which corrections are not taken into account.
 * \return Default polar motion short-period correction calculator
 */
std::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > getDefaultPolarMotionCorrectionCalculator(
        const double minimumAmplitude = 0.0 );

}  // namespace earth_orientation

}  // namespace tudat

#endif  // TUDAT_SHORTPERIODEARTHORIENTATIONCORRECTIONCALCULATOR_H
