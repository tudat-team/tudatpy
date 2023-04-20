/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_READTABULATEDMEDIACORRECTIONS_H
#define TUDAT_READTABULATEDMEDIACORRECTIONS_H

#include "tudat/astro/observation_models/corrections/tabulatedMediaCorrection.h"
#include "tudat/astro/observation_models/observableTypes.h"

namespace tudat
{

namespace input_output
{

// TODO: code to parse CSP commands is super shitty... would be nice if it wasn't

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

    double startTime_;
    double endTime_;

private:

    double convertTime( std::string yearMonthDay, std::string hoursMinutesSeconds );
};

class IonosphericCorrectionCspCommand: public AtmosphericCorrectionCspCommand
{
public:

    IonosphericCorrectionCspCommand( std::vector< std::string >& cspCommand );

    std::string sourceSpecifier_;

    int sourceId_;

private:

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

observation_models::AtmosphericCorrectionPerStationType createTroposphericCorrection(
        std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles,
        const std::string& modelIdentifier );

observation_models::AtmosphericCorrectionPerStationType createTroposphericDryCorrection(
        std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles );

observation_models::AtmosphericCorrectionPerStationType createTroposphericWetCorrection(
        std::vector< std::shared_ptr< CspRawFile > >& rawCspFiles );

} // namespace input_output

} // namespace tudat

#endif //TUDAT_READTABULATEDMEDIACORRECTIONS_H
