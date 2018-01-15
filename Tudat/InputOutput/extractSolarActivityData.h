/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      Data file:
 *                        http://celestrak.com/SpaceData/sw19571001.txt
 *                        http://celestrak.com/SpaceData/sw20110101.txt
 *      Data format explanation:
 *                        http://celestrak.com/SpaceData/SpaceWx-format.asp
 *
 */

#ifndef TUDAT_EXTRACTSOLARACTIVITY_H
#define TUDAT_EXTRACTSOLARACTIVITY_H

#include "Tudat/InputOutput/solarActivityData.h"
#include "Tudat/InputOutput/extractor.h"

#include <boost/shared_ptr.hpp>

namespace tudat
{
namespace input_output
{
namespace solar_activity
{

//! Solar activity extractor class.
/*!
 * This class extracts the numeric information from a ParsedDataLineMapPtr containing parsed
 * solar activity data designed and placec them in a container of class SolarActivityData
 */
class ExtractSolarActivityData : public tudat::input_output::Extractor<
            tudat::input_output::solar_activity::SolarActivityData >
{

public:

    //! Extracts the solar activity data to a SolarActivityData container.
    /*!
     * Extracts the solar activity data from a "ParsedDataLineMap" object and saves it in a
     * "SolarActivityData" contatiner.
     */
    boost::shared_ptr< tudat::input_output::solar_activity::SolarActivityData > extract(
                tudat::input_output::parsed_data_vector_utilities::ParsedDataLineMapPtr data );

protected:

private:

};

}   // namespace solar_activity
}   // namespace input_output
}   // namespace tudat

#endif  // TUDAT_EXTRACTSOLARACTIVITY_H
