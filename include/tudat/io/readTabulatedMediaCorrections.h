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

// TODO: code to parse CSP commands is a bit shitty... would be nice if it wasn't

class CspCommand
{
public:
    CspCommand ( )
    { }

    virtual ~CspCommand ( )
    { }

private:

};

class AtmosphericCorrectionCspCommand: public CspCommand
{
public:

    AtmosphericCorrectionCspCommand( std::vector< std::string >& cspCommand );

    std::string modelIdentifier_;

    std::string dataTypesIdentifier_;

    std::string computationSpecifier_;
    std::vector< double > computationCoefficients_;

    std::string groundStationsId_;

    std::string sourceSpecifier_;
    int sourceId_;

    double startTime_;
    double endTime_;

private:

    double convertTime( std::string yearMonthDay, std::string hoursMinutesSeconds );
};

class CspRawFile
{
public:

    CspRawFile( const std::string& cspFile );

    std::string getFileName( )
    {
        return fileName_;
    }

    std::vector< std::shared_ptr< CspCommand > > getCspCommands( )
    {
        return cspCommands_;
    }

private:

    std::vector< std::string > readCspCommandsFile( std::string file );

    std::string fileName_;

    std::vector< std::shared_ptr< CspCommand > > cspCommands_;
};

std::shared_ptr< observation_models::TabulatedMediaReferenceCorrection > createReferenceCorrection(
        double startTime,
        double endTime,
        std::vector< double > coefficients,
        std::string computationSpecifier );

std::vector< std::string > getGroundStationsNames( std::string groundStationIdentifier );

std::vector< observation_models::ObservableType > getBaseObservableTypes( std::string observableTypeIdentifier );

bool compareAtmosphericCspFileStartDate( std::shared_ptr< CspRawFile > rawCspData1,
                                         std::shared_ptr< CspRawFile > rawCspData2 );

std::string getSourceName(
        std::string& sourceSpecifier,
        int sourceId,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId,
        const std::map< int, std::string >& quasarNamePerQuasarId );

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType createTroposphericCorrection(
        std::vector< std::shared_ptr< CspRawFile > > rawCspFiles,
        const std::string& modelIdentifier,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId = std::map< int, std::string >( ),
        const std::map< int, std::string >& quasarNamePerQuasarId = std::map< int, std::string >( ) );

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType createTroposphericDryCorrection(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles );

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType createTroposphericWetCorrection(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles );

observation_models::AtmosphericCorrectionPerStationAndSpacecraftType createIonosphericCorrection(
        const std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId = std::map< int, std::string >( ),
        const std::map< int, std::string >& quasarNamePerQuasarId = std::map< int, std::string >( ) );

} // namespace input_output

} // namespace tudat

#endif //TUDAT_READTABULATEDMEDIACORRECTIONS_H
