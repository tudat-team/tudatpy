#ifndef TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H
#define TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H

#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>

// #include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

#include "tudat/io/readOceanLoadingData.h"
#include "tudat/astro/basic_astro/timeConversions.h"
#include "tudat/astro/ground_stations/bodyDeformationModel.h"
#include "tudat/interface/sofa/fundamentalArguments.h"

namespace tudat
{

namespace basic_astrodynamics
{

//! Ocean tidal loading synthesis convention.
enum class OceanTideSynthesis { iersHardisp, legacyElevenConstituents };

//! Function to read BLQ file.
/*!
 *  Function to read BLQ file containing three-dimensional local (up, west, south) site displacement
 *  amplitudes and phases at a set of sites.
 *  \param fileName Location and name of file that is to be read.
 *  \return Vector containig two maps of string vs. Matrix. First map contains site ID vs. amplitude,
 *  seconds map contains site ID vs. phase.
 */
std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > > readBlqFile( const std::string& fileName );

//! Convert local up/west/south to Cartesian using east/north/radial unit vectors.
Eigen::Vector3d convertLocalOceanTideDisplacementToPlentocentricDisplacement(
        const Eigen::Vector3d& localDisplacement,
        const std::vector< Eigen::Vector3d >& earthNorthradialUnitVectors );

//! Class to calculate site displacement to to ocean tide deformation.
/*!
 *  The default IERS HARDISP convention expands the eleven provider coefficients to
 *  342 constituents using complex admittance interpolation (IERS section 7.1.2).
 *  The explicit legacy mode retains the former eleven-constituent Moyer synthesis.
 *  BLQ data are loaded before evaluation and retain their provider frame-origin metadata.
 *  Native phase rates reject UTC offset steps/ambiguity in the current scalar-UTC
 *  representation. This is not a physical tide singularity. Explicit custom argument
 *  callbacks retain their own time basis and interval policy.
 */
class OceanTideEarthDeformation : public BodyDeformationModel
{
public:
    //! Filename adapters; all BLQ input occurs during construction or addBlqFile.
    OceanTideEarthDeformation( const std::vector< std::string >& blqFiles,
                               const OceanTideSynthesis synthesis = OceanTideSynthesis::iersHardisp );

    OceanTideEarthDeformation( const std::string& blqFile, const OceanTideSynthesis synthesis = OceanTideSynthesis::iersHardisp ):
        OceanTideEarthDeformation( std::vector< std::string >{ blqFile }, synthesis )
    {}

    //! Construct from validated, preloaded provider coefficients with no file access.
    OceanTideEarthDeformation( const input_output::OceanLoadingData& data,
                               const OceanTideSynthesis synthesis = OceanTideSynthesis::iersHardisp );

    //! Explicit custom-argument compatibility path, always using the legacy eleven constituents.
    OceanTideEarthDeformation( const std::vector< std::string >& blqFiles,
                               const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction );

    OceanTideEarthDeformation( const std::string& blqFile, const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction ):
        OceanTideEarthDeformation( std::vector< std::string >{ blqFile }, doodsonArgumentFunction )
    {}

    OceanTideEarthDeformation( const input_output::OceanLoadingData& data,
                               const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction );

    void addBlqFile( const std::string& blqFile );

    //! Add preloaded station data; duplicate station identifiers are rejected atomically.
    void addOceanLoadingData( const input_output::OceanLoadingData& data );

    const input_output::OceanLoadingData& getOceanLoadingData( ) const
    {
        return oceanLoadingData_;
    }

    //! Spherical east/north/radial displacement in metres; time is TT seconds since J2000.
    Eigen::Vector3d calculateDisplacementInEnuFrame( const double time, const std::string& siteIdentifier );

    //! True station topocentric ENU (ellipsoid-normal up for oblate station shapes), in metres.
    Eigen::Vector3d calculateDisplacementInEnuFrame( const double time,
                                                     const std::shared_ptr< ground_stations::GroundStationState > stationState );

    OceanTideSynthesis getSynthesis( ) const
    {
        return synthesis_;
    }

    Eigen::Vector3d calculateDisplacement( const double time, const Eigen::Vector3d& bodyFixedPosition )
    {
        throw std::runtime_error( "Error, cannot compute ocean tide displacements from position alone." );
    }

    //! Function to calculate site displacement at given time.
    /*!
     * Function to calculate site displacement at given time at given site. Note that the site identifier needs to be
     * identical to that provided in the BLQ file. If doodson arguments are available, it is recommended to use the second,
     * overloaded version of this function for performance reasons.
     * \param time TT seconds since J2000. HARDISP uses the corresponding UTC time of day.
     * Standalone calls require TT. The existing environment passes its TDB time as TT,
     * retaining the millisecond TDB-TT approximation until explicit time-scale routing is implemented.
     * The Cartesian result uses the geocentric spherical basis, independently of the station shape.
     * \return Displacement of requested site, at requested time.
     */
    Eigen::Vector3d calculateDisplacement( const double time, const std::shared_ptr< ground_stations::GroundStationState > stationState );

    //! Function to retrieve the map containing site displacement amplitudes per tide.
    /*!
     *  Function to retrieve the map containing site displacement amplitudes per tide.
     */
    std::map< std::string, Eigen::Matrix< double, 3, 11 > > getAmplitudesMap( )
    {
        return siteOceanTideAmplitudes_;
    }

    //! Function to retrieve the map containing site displacement phase per tide.
    /*!
     *  Function to retrieve the map containing site displacement phase per tide.
     */
    std::map< std::string, Eigen::Matrix< double, 3, 11 > > getPhasesMap( )
    {
        return siteOceanTidePhases_;
    }

    //! Function to retrieve matrix containing the tide schwiderski factors.
    /*!
     *  Function to retrieve matrix containing the tide schwiderski factors.
     */
    Eigen::Matrix< double, 11, 1 > getSchwiderskiFactors( )
    {
        return schwiderskiFactors_;
    }

    //! Function to retrieve the matrix containing the tide doodson multipliers.
    /*!
     *  Function to retrieve the matrix containing the tide doodson multipliers.
     */
    Eigen::Matrix< double, 11, 4 > getDoodsonMultipliers( )
    {
        return doodsonMultipliers_;
    }

protected:
    void validateDisplacementDerivativeInterval( const double time, const double lower, const double upper ) const override;

private:
    void initializeLegacyConstituents( );
    Eigen::Vector3d calculateGeocentricEnuDisplacement( const double time, const std::string& siteIdentifier );
    OceanTideSynthesis synthesis_;
    input_output::OceanLoadingData oceanLoadingData_;

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
    bool usesNativeUtcPhase_ = true;
};

}  // namespace basic_astrodynamics

}  // namespace tudat

#endif  // TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H
