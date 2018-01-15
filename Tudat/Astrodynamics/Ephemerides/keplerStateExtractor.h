/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_KEPLER_STATE_EXTRACTOR_H
#define TUDAT_KEPLER_STATE_EXTRACTOR_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/InputOutput/extractor.h"

#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace ephemerides
{

//! Keplerian State extractor class.
/*!
 * This class contains the functionality of extracting the Keplerian elements from a parsed data
 * line map.
 */
class KeplerStateExtractor : public input_output::Extractor< Eigen::Vector6d >
{
public:

    //! Extract the Keplerian Elements.
    /*!
     * Returns a KeplerianElements object containing the orbital parameters found in the input data
     * line map. If one of the elements is not present, the function throws an exception.
     * If no true anomaly field type is found, the function searches for the mean
     * anomaly field type and performs the necessary conversions. The if-statements can be easily
     * extended to deal with other non-standard field types.
     *
     * \param dataLineMap Data map corresponding to a parsed line.
     * \return A KeplerianElements object containing the orbital parameters found in the data map.
     */
    boost::shared_ptr< Eigen::Vector6d > extract( ParsedDataLineMapPtr dataLineMap );

protected:

private:
};

//! Typedef for shared-pointer to KeplerStateExtractor object.
typedef boost::shared_ptr< KeplerStateExtractor > KeplerStateExtractorPointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_KEPLER_STATE_EXTRACTOR_H
