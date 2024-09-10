#ifndef TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H
#define TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H

#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>

//#include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/interface/sofa/fundamentalArguments.h"

namespace tudat
{

namespace basic_astrodynamics
{


//! Function to read BLQ file.
/*!
 *  Function to read BLQ file containing three-dimensional local (up, west, south) site displacement
 *  amplitudes and phases at a set of sites.
 *  \param fileName Location and name of file that is to be read.
 *  \return Vector containig two maps of string vs. Matrix. First map contains site ID vs. amplitude,
 *  seconds map contains site ID vs. phase.
 */
std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > > readBlqFile( const std::string& fileName );

Eigen::Vector3d convertLocalOceanTideDisplacementToPlentocentricDisplacement(
        const Eigen::Vector3d& localDisplacement,
        const std::vector< Eigen::Vector3d >& earthNorthradialUnitVectors );

//! Class to calculate site displacement to to ocean tide deformation.
/*!
 *  Class to calculate site displacement to to ocean tide deformation. Only the 11 main tidal constituents are include,
 *  i.e. Eqs. 5-85 to 5-87 of Moyer(2000), not the full interpolation as adviced by the IERS 2010 Conventions, section 7.1.2.
 *  Object requires data typically obtained from BLQ files generated at http://froste.oso.chalmers.se/loading/.
 */
class OceanTideEarthDeformation: public BodyDeformationModel
{
public:
    //! Constructor, reads and initializes tide displacement amplitudes and phases.
    /*!
     *  Constructor, reads and initializes tide displacement amplitudes and phases at number of sites provided
     *  by BLQ file.
     *  \param blqFile BLQ file containing ocean tide displacement amplitudes and phases.
     */
    OceanTideEarthDeformation( const std::vector< std::string >& blqFiles,
                               const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction =
                                   [](const double time ){ return sofa_interface::calculateDoodsonFundamentalArguments( time ); } );

    OceanTideEarthDeformation( const std::string& blqFile,
                               const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction =
                                   [](const double time ){ return sofa_interface::calculateDoodsonFundamentalArguments( time ); } ):
                                   OceanTideEarthDeformation( std::vector< std::string >( { blqFile } ), doodsonArgumentFunction ){ }

    void addBlqFile( const std::string& blqFile );

    Eigen::Vector3d calculateDisplacementInEnuFrame(
        const double time,
        const std::string& siteIdentifier );

    Eigen::Vector3d calculateDisplacement(
        const double time,
        const Eigen::Vector3d& bodyFixedPosition )
    {
        throw std::runtime_error( "Error, cannot compute ocean tide displacements from position alone." );
    }

    //! Function to calculate site displacement at given time.
    /*!
     * Function to calculate site displacement at given time at given site. Note that the site identifier needs to be
     * identical to that provided in the BLQ file. If doodson arguments are available, it is recommended to use the second,
     * overloaded version of this function for performance reasons.
     * \param ephemerisTime Time in seconds since J2000. Due to the rate of change and amplitude of the effect, the time
     * scale in which this is provided is unlikely to be consequential, to measurement precision, on the resulting displacement.
     * \return Displacement of requested site, at requested time.
     */
    Eigen::Vector3d calculateDisplacement(
        const double time,
        const std::shared_ptr< ground_stations::GroundStationState > stationState );

    //! Function to retrieve the map containing site displacement amplitudes per tide.
    /*!
     *  Function to retrieve the map containing site displacement amplitudes per tide.
     */
    std::map< std::string, Eigen::Matrix< double, 3, 11 > > getAmplitudesMap( )
    { return siteOceanTideAmplitudes_; }

    //! Function to retrieve the map containing site displacement phase per tide.
    /*!
     *  Function to retrieve the map containing site displacement phase per tide.
     */
    std::map< std::string, Eigen::Matrix< double, 3, 11 > > getPhasesMap( )
    { return siteOceanTidePhases_; }

    //! Function to retrieve matrix containing the tide schwiderski factors.
    /*!
     *  Function to retrieve matrix containing the tide schwiderski factors.
     */
    Eigen::Matrix< double, 11, 1 > getSchwiderskiFactors( )
    { return schwiderskiFactors_; }

    //! Function to retrieve the matrix containing the tide doodson multipliers.
    /*!
     *  Function to retrieve the matrix containing the tide doodson multipliers.
     */
    Eigen::Matrix< double, 11, 4 > getDoodsonMultipliers( )
    { return doodsonMultipliers_; }

private:
    //! Map containing site displacement amplitudes per tide.
    /*!
     *  Map containing site displacement amplitudes per tide.
     */
    std::map< std::string, Eigen::Matrix< double, 3, 11 > > siteOceanTideAmplitudes_;

    //! Map containing site displacement phase per tide.
    /*!
     *  Mmap containing site displacement phase per tide.
     */
    std::map< std::string, Eigen::Matrix< double, 3, 11 > > siteOceanTidePhases_;

    //! Matrix containing the tide schwiderski factors.
    /*!
     *  Matrix containing the tide schwiderski factors.
     */
    Eigen::Matrix< double, 11, 1 > schwiderskiFactors_;

    //! Matrix containing the tide doodson multipliers.
    /*!
     *  Matrix containing the tide doodson multipliers.
     */
    Eigen::Matrix< double, 11, 4 > doodsonMultipliers_;

    std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction_;

};

}

}

#endif // TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H
