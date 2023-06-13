/*    Copyright (c) 2010-2023, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEATMOSPHERICLIGHTTIMECORRECTION_H
#define TUDAT_CREATEATMOSPHERICLIGHTTIMECORRECTION_H

#include "tudat/io/readTabulatedMediaCorrections.h"
#include "tudat/astro/observation_models/corrections/tabulatedMediaCorrection.h"

namespace tudat
{

namespace observation_models
{

typedef std::map< std::pair< std::string, std::string >, std::map< observation_models::ObservableType,
    std::shared_ptr< TabulatedMediaReferenceCorrectionManager > > > AtmosphericCorrectionPerStationAndSpacecraftType;

/*!
 * Creates the object to compute the media correction, given the information in a CSP command.
 *
 * @param startTime Start time of CSP command.
 * @param endTime End time of CSP command.
 * @param coefficients Coefficients of computation model.
 * @param computationSpecifier Computation model specifier.
 * @return Object to compute DSN tabulated media correction.
 */
std::shared_ptr< TabulatedMediaReferenceCorrection > createReferenceCorrection(
        const double startTime,
        const double endTime,
        const std::vector< double >& coefficients,
        const std::string& computationSpecifier );

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
AtmosphericCorrectionPerStationAndSpacecraftType extractAtmosphericCorrection(
        std::vector< std::shared_ptr< input_output::CspRawFile > > rawCspFiles,
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
AtmosphericCorrectionPerStationAndSpacecraftType extractTroposphericDryCorrectionAdjustment(
        const std::vector< std::shared_ptr< input_output::CspRawFile > >& rawCspFiles );

/*!
 * Creates the objets to compute the wet tropospheric correction adjustments, based on the data in the provided CSP files.
 *
 * @param rawCspFiles Vector of CSP files.
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type. The
 *      corrections don't depend on the data source, hence the spacecraft is always set to empty string ("").
 */
AtmosphericCorrectionPerStationAndSpacecraftType extractTroposphericWetCorrectionAdjustment(
        const std::vector< std::shared_ptr< input_output::CspRawFile > >& rawCspFiles );

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
AtmosphericCorrectionPerStationAndSpacecraftType extractIonosphericCorrection(
        const std::vector< std::shared_ptr< input_output::CspRawFile > >& rawCspFiles,
        const std::map< int, std::string >& spacecraftNamePerSpacecraftId = std::map< int, std::string >( ),
        const std::map< int, std::string >& quasarNamePerQuasarId = std::map< int, std::string >( ) );

/*!
 * Creates the objects to compute the default tropospheric dry corrections, according to the DSN default seasonal model.
 *
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type. The
 *      corrections don't depend on the data source, hence the spacecraft is always set to empty string ("").
 */
AtmosphericCorrectionPerStationAndSpacecraftType extractDefaultTroposphericDryCorrection( );

/*!
 * Creates the objects to compute the default tropospheric wet corrections, according to the DSN default seasonal model.
 *
 * @return Atmospheric corrections mapped by a pair of (ground station, spacecraft) and by the observable type. The
 *      corrections don't depend on the data source, hence the spacecraft is always set to empty string ("").
 */
AtmosphericCorrectionPerStationAndSpacecraftType extractDefaultTroposphericWetCorrection( );


} // namespace observation_models

} // namespace tudat

#endif //TUDAT_CREATEATMOSPHERICLIGHTTIMECORRECTION_H
