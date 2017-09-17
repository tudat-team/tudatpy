#ifndef SHORTPERIODEARTHORIENTATIONCORRECTIONCALCULATOR_H
#define SHORTPERIODEARTHORIENTATIONCORRECTIONCALCULATOR_H


#include <string>

#include <boost/function.hpp>
#include <boost/bind.hpp>
#include <boost/make_shared.hpp>

#include <Eigen/Core>

#include "Tudat/Basics/basicTypedefs.h"
#include "Tudat/InputOutput/basicInputOutput.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/unitConversions.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/EarthOrientation/readAmplitudeAndDoodsonNumber.h"
#include "Tudat/External/SofaInterface/fundamentalArguments.h"
#include "Tudat/InputOutput/basicInputOutput.h"


//! Object to calculate the short period (2 days and lower) polar motion.
/*!
 *  Object to calculate the short period (2 days and lower) polar motion, taking into account
 *  variations due to both libration and ocean tides.
 */
namespace tudat
{

namespace earth_orientation
{

template< typename OutputType >
class ShortPeriodEarthOrientationCorrectionCalculator
{
public:
    //! Constructor, requires files defining the amplitudes and doodson numbers of the variations.
    /*!
     *  Constructor, requires files defining the sine and cosine amplitudes and doodson numbers o
     *  of the variations for both ocean tide and libration variations. A cut-off amplitude (taken as
     *  the RSS of the amplitudes) below which amplitudes in the files are ignored can be provided.
     *  \param minimumAmplitude Minimum amplitude that is read from files and considered in calculations.
     *  Default is zero, i.e. all corrections are accepted.
     *  \param librationAmplitudesFile Amplitudes of libration-induced variations, defaults from
     *  2010 IERS Conventions, Tables 5.1a.
     *  \param librationDoodsonMultipliersFile Doodson number of libration-induced variations, defaults from
     *  2010 IERS Conventions, Tables 5.1a
     *  \param oceanTideAmplitudesFile Amplitudes of ocean tide-induced variations, defaults from
     *  2010 IERS Conventions, Tables 8.2a and 8.2b.
     *  \param oceanTideDoodsonMultipliersFile Doodson number of ocean tide-induced variations, defaults from
     *  2010 IERS Conventions, Tables 8.2a and 8.2b.
     */
    ShortPeriodEarthOrientationCorrectionCalculator(
            const double conversionFactor,
            const double minimumAmplitude = 0.0,
            const std::vector< std::string >& amplitudesFiles =
    { tudat::input_output::getEarthOrientationDataFilesPath( ) + "polarMotionLibrationAmplitudesQuasiDiurnalOnly.txt",
            tudat::input_output::getEarthOrientationDataFilesPath( ) + "polarMotionOceanTidesAmplitudes.txt", },
            const std::vector< std::string >& argumentMultipliersFile =
    { tudat::input_output::getEarthOrientationDataFilesPath( ) + "polarMotionLibrationDoodsonMultipliersQuasiDiurnalOnly.txt",
            tudat::input_output::getEarthOrientationDataFilesPath( ) + "polarMotionOceanTidesDoodsonMultipliers.txt" },
            const boost::function< Eigen::Vector6d( const double ) > argumentFunction =
            boost::bind( &sofa_interface::calculateDelaunayFundamentalArgumentsWithGmst, _1 ) ):
        argumentFunction_( argumentFunction )
    {
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd > dataFromFile;

        for( unsigned int i = 0; i < amplitudesFiles.size( ); i++ )
        {
            dataFromFile = readAmplitudesAndDoodsonMultipliers(
                        amplitudesFiles.at( i ), argumentMultipliersFile.at( i ), minimumAmplitude );
            argumentAmplitudes_.push_back( conversionFactor * dataFromFile.first );
            argumentMultipliers_.push_back( dataFromFile.second );
        }
    }

    //! Function to obtain short period polar motion corrections.
    /*!
     *  Function to obtain short period polar motion corrections, using time as input. Doodson arguments are calculated internally.
     *  \param ephemerisTime Time (TDB seconds since J2000) at which corretions are to be determined
     *  \return Short period corrections to x_{p} and y_{p}
     */
    OutputType getCorrections( const double& ephemerisTime )
    {
        return sumCorrectionTerms( argumentFunction_( ephemerisTime ) );
    }

    //! Function to obtain short period polar motion corrections.
    /*!
     *  Function to obtain short period polar motion corrections, using doodson arguments as input.
     *  \param doodsonArguments Doodson arguments from which corretions are to be determined
     *  \return Short period corrections to x_{p} and y_{p}
     */
    OutputType getCorrections( const Eigen::Vector6d& arguments )
    {
        return sumCorrectionTerms( arguments );
    }

private:

    //! Function to sum all the corretion terms to the polar motion.
    /*!
     *  Function to sum all the corretion terms to the polar motion.
     */
    OutputType sumCorrectionTerms( const Eigen::Vector6d& arguments );

    //! Amplitudes of libration-induced variations
    /*!
     *  Amplitudes of libration-induced variations
     */
    std::vector< Eigen::MatrixXd > argumentAmplitudes_;

    //! Doodsons number of libration-induced variations.
    /*!
     *  Doodsons number of libration-induced variations.
     */
    std::vector< Eigen::MatrixXd > argumentMultipliers_;

    boost::function< Eigen::Vector6d( const double ) > argumentFunction_;


};


boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< double > > getDefaultUT1CorrectionCalculator(
        const double minimumAmplitude = 0.0 );

boost::shared_ptr< ShortPeriodEarthOrientationCorrectionCalculator< Eigen::Vector2d > > getDefaultPolarMotionCorrectionCalculator(
        const double minimumAmplitude = 0.0 );

}

}


#endif // SHORTPERIODEARTHORIENTATIONCORRECTIONCALCULATOR_H
