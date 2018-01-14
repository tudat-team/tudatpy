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

#ifndef TUDAT_CARTESIAN_STATE_EXTRACTOR_H
#define TUDAT_CARTESIAN_STATE_EXTRACTOR_H

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/InputOutput/extractor.h"
#include "Tudat/Basics/basicTypedefs.h"

namespace tudat
{
namespace ephemerides
{

//! Cartesian State extractor class.
/*!
 * This class contains the functionality of extracting the Cartesian elements from a parsed data
 * line map.
 */
class CartesianStateExtractor : public input_output::Extractor< Eigen::Vector6d >
{
public:

    //! Extract the Cartesian elements.
    /*!
     * Returns a CartesianElements object containing the cartesian elements found in the input data
     * line map. If one of the elements is not present, the function throws an exception.
     *
     * \param dataLineMap Data map corresponding to a parsed line.
     * \return A CartesianElements object containing the orbital parameters found in the data map.
     */
    boost::shared_ptr< Eigen::Vector6d > extract( ParsedDataLineMapPtr dataLineMap );

protected:

private:
};

//! Typedef for shared-pointer to CartesianStateExtractor object.
typedef boost::shared_ptr< CartesianStateExtractor > CartesianStateExtractorPointer;

} // namespace ephemerides
} // namespace tudat

#endif // TUDAT_CARTESIAN_STATE_EXTRACTOR_H
