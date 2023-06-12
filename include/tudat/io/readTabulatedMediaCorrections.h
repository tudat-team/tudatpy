/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References:
 *          820-013 TRK-2-23, Media Calibration Interface, Revision C (2008), DSN/JPL
 */

#ifndef TUDAT_READTABULATEDMEDIACORRECTIONS_H
#define TUDAT_READTABULATEDMEDIACORRECTIONS_H

#include "tudat/astro/observation_models/corrections/tabulatedMediaCorrection.h"
#include "tudat/astro/observation_models/observableTypes.h"

namespace tudat
{

namespace input_output
{

// Base class defining a Control Statement Processor (CSP) command
class CspCommand
{
public:
    CspCommand ( )
    { }

    virtual ~CspCommand ( )
    { }

private:

};

// Representation of a CSP command with information about atmospheric corrections, according to TRK-2-23
class AtmosphericCorrectionCspCommand: public CspCommand
{
public:

    /*!
     * Constructor. Parses a CSP command, saving its information.
     * @param cspCommand Vector with CSP command elements and respective values.
     */
    AtmosphericCorrectionCspCommand( std::vector< std::string >& cspCommand );

    /*!
     * Constructor.
     */
    AtmosphericCorrectionCspCommand( ):
        sourceSpecifier_( "" ),
        sourceId_( 0 )
    { }

    // CSP identifier of the model. Can take the values DRY NUPART (dry part of troposphere), WET NUPART (wet part of troposphere)
    // or CHPART (ionosphere).
    std::string modelIdentifier_;

    // CSP identifier of the data types to which the command applies.
    std::string dataTypesIdentifier_;

    // CSP identifier of the type of model used to compute the correction.
    std::string computationSpecifier_;
    // Coefficients of the model.
    std::vector< double > computationCoefficients_;

    // Identification of the DSN station or complex
    std::string groundStationsId_;

    // Identification of the source type: SCID for spacecraft or QUASAR for quasar. Only applicable for ionosphere correctins
    std::string sourceSpecifier_;
    // Id of the source
    int sourceId_;

    // Start time of command validity
    double startTime_;
    // End time of command validity
    double endTime_;

    /*!
     * Converts the provided time into seconds since J2000. Doesn't convert the UTC time to TDB as the difference between
     * the two is negligible for computing the atmospheric corrections.
     *
     * @param yearMonthDay Data in the format YYMMDD, UTC
     * @param hoursMinutesSeconds Time in the format HH:MM:SSSSS, UTC
     * @return Time (seconds) since J2000
     */
    double convertTime( std::string yearMonthDay, std::string hoursMinutesSeconds );

private:

};

// Control Statement Processor (CSP) file.
class CspRawFile
{
public:

    /*!
     * Constructor. Reads file and parses the CSP commands it contains.
     *
     * @param cspFile File name.
     */
    CspRawFile( const std::string& cspFile ):
        fileName_( cspFile )
    {
        // Read CSP commands. Each command saved as one string
        std::vector< std::string > cspCommandsVector = readCspCommandsFile( fileName_ );

        // Parse commands
        parseCspCommands( cspCommandsVector );
    }

    /*!
     * Constructor. Parses the provided CSP commands.
     *
     * @param cspCommandsVector Vector of CSP commands.
     */
    CspRawFile( const std::vector< std::string >& cspCommandsVector ):
        fileName_( "" )
    {
        // Parse commands
        parseCspCommands( cspCommandsVector );
    }

    /*!
     * Constructor. Saves the provided CSP commands.
     *
     * @param cspCommands Vector of parsed CSP commands.
     */
    CspRawFile( const std::vector< std::shared_ptr< CspCommand > >& cspCommands ):
        fileName_( "" ),
        cspCommands_( cspCommands )
    { }

    // Returns the name of the file.
    std::string getFileName( )
    {
        return fileName_;
    }

    // Returns the parsed CSP commands that were contained in the file.
    std::vector< std::shared_ptr< CspCommand > > getCspCommands( )
    {
        return cspCommands_;
    }

private:

    /*!
     * Reads CSP commands file, placing each command into a string.
     *
     * @param file File name.
     * @return Vector of CSP commands.
     */
    std::vector< std::string > readCspCommandsFile( const std::string& file );

    /*!
     * Parses a vector of CSP commands and saves them.
     *
     * @param cspCommandsVector Vector of CSP commands.
     */
    void parseCspCommands( const std::vector< std::string >& cspCommandsVector );

    // Name of the file
    std::string fileName_;

    // Vector of CSP commands
    std::vector< std::shared_ptr< CspCommand > > cspCommands_;
};

/*!
 * Creates the object to compute the media correction, given the information in a CSP command.
 *
 * @param startTime Start time of CSP command.
 * @param endTime End time of CSP command.
 * @param coefficients Coefficients of computation model.
 * @param computationSpecifier Computation model specifier.
 * @return Object to compute DSN tabulated media correction.
 */
std::shared_ptr< observation_models::TabulatedMediaReferenceCorrection > createReferenceCorrection(
        const double startTime,
        const double endTime,
        const std::vector< double >& coefficients,
        const std::string& computationSpecifier );

/*!
 * Returns the names of the ground stations associated with a given CSP ground station identifier. The identifier may
 * correspond to a single ground station or to a complex.
 *
 * @param groundStationIdentifier CSP ground station identifier.
 * @return Vector with the names of the ground stations.
 */
std::vector< std::string > getGroundStationsNames( const std::string& groundStationIdentifier );

/*!
 * Returns the vector of base observable types associated with each CSP observable identifier. Each CSP observable identifier
 * can correspond to one or more tudat observable types.
 *
 * @param observableTypeIdentifier CSP observable identifier.
 * @return Vector of observable types.
 */
std::vector< observation_models::ObservableType > getBaseObservableTypes( const std::string& observableTypeIdentifier );

/*!
 * Checks which file starts first. Used to sort CSP files. Returns true if rawCspData1 starts first, false otherwise.
 *
 * @param rawCspData1 CSP file.
 * @param rawCspData2 CSP file.
 * @return
 */
bool compareAtmosphericCspFileStartDate( std::shared_ptr< CspRawFile > rawCspData1,
                                         std::shared_ptr< CspRawFile > rawCspData2 );

/*!
 * Returns the name of the data source.
 * If the source ID is defined in one of the provided maps, the function returns the name in the map.
 * If the source is a spacecraft and its ID isn't in the spacecraftNamePerSpacecraftId map, an error is thrown.
 * If the source is a quasar and its ID isn't in the quasarNamePerQuasarId map, an error is thrown.
 *
 * @param sourceSpecifier CSP source specifier (SCID or QUASAR)
 * @param sourceId Source ID.
 * @param spacecraftNamePerSpacecraftId Map with the name associated with each spacecraft ID
 * @param quasarNamePerQuasarId Map with the name associated with each quasar ID
 * @return Name of the source
 */
std::string getSourceName(
        const std::string& sourceSpecifier,
        const int sourceId,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId = std::map< int, std::string >( ),
        const std::map< int, std::string >& quasarNamePerQuasarId = std::map< int, std::string >( ) );

/*!
 * Creates the objects to compute the specified atmospheric corrections using the data in the provided CSP files.
 * spacecraftNamePerSpacecraftId and quasarNamePerQuasarId are only necessary for ionospheric corrections (the case in which
 * the data is provided for a specific spacecraft/quasar).
 * Function should not be called directly. Instead, one should call the function associated with the desired correction
 * type (dry troposphere, wet troposhere, or ionosphere).
 *
 * @param rawCspFiles Vector of CSP files.
 * @param modelIdentifier CSP model identifier (DRY NUPART, WET NUPART, or CHPART)
 * @param spacecraftNamePerSpacecraftId Map with the name of the body associated with each spacecraft ID in the CSP file.
 * @param quasarNamePerQuasarId Map with the name of the body associated with each quasar ID in the CSP file.
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type.
 */
observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractAtmosphericCorrection(
        std::vector< std::shared_ptr< CspRawFile > > rawCspFiles,
        const std::string& modelIdentifier,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId = std::map< int, std::string >( ),
        const std::map< int, std::string >& quasarNamePerQuasarId = std::map< int, std::string >( ) );

/*!
 * Creates the objets to compute the dry tropospheric correction adjustments, based on the data in the provided CSP files.
 *
 * @param rawCspFiles Vector of CSP files.
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type. The
 *      corrections don't depend on the data source, hence the spacecraft is always set to empty string ("").
 */
observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractTroposphericDryCorrectionAdjustment(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles );

/*!
 * Creates the objets to compute the wet tropospheric correction adjustments, based on the data in the provided CSP files.
 *
 * @param rawCspFiles Vector of CSP files.
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type. The
 *      corrections don't depend on the data source, hence the spacecraft is always set to empty string ("").
 */
observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractTroposphericWetCorrectionAdjustment(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles );

/*!
 * Creates the objets to compute the ionospheric corrections, based on the data in the provided CSP files.
 * If the source is a spacecraft and its ID isn't in the spacecraftNamePerSpacecraftId map, an error is thrown.
 * If the source is a quasar and its ID isn't in the quasarNamePerQuasarId map, an error is thrown.
 *
 * @param rawCspFiles Vector of CSP files.
 * @param spacecraftNamePerSpacecraftId Map with the name of the body associated with each spacecraft ID in the CSP file.
 * @param quasarNamePerQuasarId Map with the name of the body associated with each quasar ID in the CSP file.
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type.
 */
observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractIonosphericCorrection(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId = std::map< int, std::string >( ),
        const std::map< int, std::string >& quasarNamePerQuasarId = std::map< int, std::string >( ) );

/*!
 * Returns a CSP commands file corresponding to the DSN default tropospheric seasonal model, according to Estefan and
 * Sovers (1994), Fig. 3b.
 * @return
 */
std::shared_ptr< CspRawFile > getDsnDefaultTroposphericSeasonalModelCspFile( );

/*!
 * Creates the objects to compute the default tropospheric dry corrections, according to the DSN default seasonal model.
 *
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type. The
 *      corrections don't depend on the data source, hence the spacecraft is always set to empty string ("").
 */
observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractDefaultTroposphericDryCorrection( );

/*!
 * Creates the objects to compute the default tropospheric wet corrections, according to the DSN default seasonal model.
 *
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type. The
 *      corrections don't depend on the data source, hence the spacecraft is always set to empty string ("").
 */
observation_models::AtmosphericCorrectionPerStationAndSpacecraftType extractDefaultTroposphericWetCorrection( );

} // namespace input_output

} // namespace tudat

#endif //TUDAT_READTABULATEDMEDIACORRECTIONS_H
