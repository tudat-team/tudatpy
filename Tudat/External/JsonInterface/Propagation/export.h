/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_EXPORT_H
#define TUDAT_JSONINTERFACE_EXPORT_H

#include "variable.h"

#include "Tudat/External/JsonInterface/Support/valueAccess.h"
#include "Tudat/External/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace simulation_setup
{

class ExportSettings
{
public:
    //! Constructor.
    ExportSettings( const boost::filesystem::path& outputFile,
                    const std::vector< boost::shared_ptr< propagators::VariableSettings > >& variables ) :
        outputFile( outputFile ), variables( variables ) { }

    //! Destructor.
    virtual ~ExportSettings( ) { }


    //! Path of the output file the results will be written to.
    boost::filesystem::path outputFile;

    //! Variables to export.
    //! The variables will be exported to a table in which each row corresponds to an epoch,
    //! and the columns contain the values of the variables specified in this vector in the provided order.
    std::vector< boost::shared_ptr< propagators::VariableSettings > > variables;

    //! Whether to include the epochs in the first column of the results table.
    bool epochsInFirstColumn = true;

    //! Number of significant digits for the exported results.
    unsigned int numericalPrecision = 10;

    //! Whether to print only the values corresponding to the initial integration step.
    bool onlyInitialStep = false;

    //! Whether to print only the values corresponding to the final integration step.
    bool onlyFinalStep = false;
};

//! Create a `json` object from a shared pointer to a `ExportSettings` object.
void to_json( json& jsonObject, const boost::shared_ptr< ExportSettings >& saveSettings );

//! Create a shared pointer to a `ExportSettings` object from a `json` object.
void from_json( const json& jsonObject, boost::shared_ptr< ExportSettings >& saveSettings );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_EXPORT_H
