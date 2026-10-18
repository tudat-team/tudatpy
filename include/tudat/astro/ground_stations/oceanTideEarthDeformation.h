#ifndef TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H
#define TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H

#include <map>
#include <vector>
#include <string>

#include <Eigen/Core>

// #include "Tudat/Mathematics/BasicMathematics/basicMathematicsFunctions.h"

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
 *  BLQ data are loaded at construction (or by addBlqFile).
 *  Native phase rates reject UTC offset steps/ambiguity in the current scalar-UTC
 *  representation. This is not a physical tide singularity. Explicit custom argument
 *  callbacks retain their own time basis and interval policy.
 */
class OceanTideEarthDeformation : public BodyDeformationModel
{
public:
    //! Constructor from BLQ files.
    OceanTideEarthDeformation( const std::vector< std::string >& blqFiles,
                               const OceanTideSynthesis synthesis = OceanTideSynthesis::iersHardisp,
                               const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter = nullptr );

    OceanTideEarthDeformation( const std::string& blqFile,
                               const OceanTideSynthesis synthesis = OceanTideSynthesis::iersHardisp,
                               const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter = nullptr ):
        OceanTideEarthDeformation( std::vector< std::string >{ blqFile }, synthesis, converter )
    {}

    //! Constructor from BLQ data, as returned by readBlqFile (amplitudes in metres and phase lags in radians per station).
    OceanTideEarthDeformation( const std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > >& blqData,
                               const OceanTideSynthesis synthesis = OceanTideSynthesis::iersHardisp,
                               const std::shared_ptr< earth_orientation::TerrestrialTimeScaleConverter >& converter = nullptr );

    //! Explicit custom-argument compatibility path, always using the legacy eleven constituents.
    OceanTideEarthDeformation( const std::vector< std::string >& blqFiles,
                               const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction );

    OceanTideEarthDeformation( const std::string& blqFile, const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction ):
        OceanTideEarthDeformation( std::vector< std::string >{ blqFile }, doodsonArgumentFunction )
    {}

    //! Constructor from BLQ data with a custom argument function, using the legacy eleven constituents.
    OceanTideEarthDeformation( const std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > >& blqData,
                               const std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction ):
        OceanTideEarthDeformation( std::vector< std::string >( ), doodsonArgumentFunction )
    {
        addBlqData( blqData );
    }

    void addBlqFile( const std::string& blqFile );

    //! Add BLQ data, as returned by readBlqFile; duplicate stations are rejected.
    void addBlqData( const std::vector< std::map< std::string, Eigen::Matrix< double, 3, 11 > > >& blqData );

    OceanTideSynthesis getSynthesis( ) const
    {
        return synthesis_;
    }

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
    Eigen::Vector3d evaluateDisplacement( const Time& tdb, const std::shared_ptr< ground_stations::GroundStationState > station ) override;
    Eigen::Vector6d evaluateDisplacementState( const Time& tdb,
                                               const std::shared_ptr< ground_stations::GroundStationState > station ) override;

private:
    Eigen::Vector3d displacementAtArgumentTime( double time, const std::shared_ptr< ground_stations::GroundStationState > station );
    void initializeLegacyConstituents( );
    Eigen::Vector3d calculateGeocentricEnuDisplacement( const double time, const std::string& siteIdentifier );

    //! Method used to compute the displacement from the coefficients.
    OceanTideSynthesis synthesis_;

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

    //! Function returning the Doodson arguments at a given epoch.
    std::function< Eigen::Vector6d( const double ) > doodsonArgumentFunction_;

    //! Whether the tidal arguments are computed from UTC by this class (false if the user supplied the argument function).
    bool usesNativeUtcPhase_ = true;
};

}  // namespace basic_astrodynamics

}  // namespace tudat

#endif  // TUDAT_TUDAT_OCEANTIDEEARTHDEFORMATION_H
