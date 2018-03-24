/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_SPICE_H
#define TUDAT_JSONINTERFACE_SPICE_H

#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/JsonInterface/Support/valueAccess.h"
#include "Tudat/JsonInterface/Support/valueConversions.h"

namespace tudat
{

namespace json_interface
{

// SpiceSettings

//! Class containing the settings for Spice used in a simulation.
/*!
 * Class containing the settings for Spice used in a simulation.
 */
class SpiceSettings
{
public:

    //! Empty constructor.
    SpiceSettings( ) { }

    //! Destructor.
    virtual ~SpiceSettings( ) { }


    //! Whether the standard kernel should be preloaded.
    bool useStandardKernels_ = true;

    //! Vector containing the paths to the additional spice kernel files to be loaded if using standard kernels.
    std::vector< boost::filesystem::path > alternativeKernels_;

    //! Vector containing the paths to the spice kernel files to be loaded if using standard kernels.
    std::vector< boost::filesystem::path > kernels_;

    //! Whether all the data from the Spice kernels should be preloaded before the simulation for the interval start
    //! epoch to end epoch (true), or whether the data from Spice should be accessed on request at every step (false).
    /*!
     * Whether all the data from the Spice kernels should be preloaded before the simulation for the interval start
     * epoch to end epoch (true), or whether the data from Spice should be accessed on request at every step (false).
     * <br/>
     * Preloading Spice data generally results in faster propagations, unless:
     * <br/>
     * <ul>
     *  <li>The simulation ends much earlier than the specified maximum simulation end epoch.</li>
     *  <li>The integrator step-size is very large (in the order of several hours or days).</li>
     * </ul>
     */
    bool preloadEphemeris_ = true;

    //! Offsets for the interval for which the spice kernels are to be preloaded.
    /*!
     * Offsets for the interval for which the spice kernels are to be preloaded.
     * <br/>
     * The kernels will be interpolated for the interval:
     * `[ initialEpoch - interpolationOffsets_.first, finalEpoch + interpolationOffsets_.second ]`
     * \remark Ignored if SpiceSettings::preloadEphemeris_ is set to `false`.
     * \remark If not specified, the used values are 10 * interpolationStep_.
     */
    std::pair< double, double > interpolationOffsets_ = { TUDAT_NAN, TUDAT_NAN };

    //! Step-size for the interpolated Spice ephemeris.
    /*!
     * Step-size for the interpolated Spice ephemeris. Ignored if preloadEphemeris_ set to false.
     */
    double interpolationStep_ = 300.0;

    //! Get initial offset for the interpolated Spice ephemeris.
    /*!
     * @copybrief getInitialOffset
     * \remark If not defined by the user (i.e. is NaN), returns `10 * interpolationStep_`.
     * \return Initial offset for the interpolated Spice ephemeris.
     */
    double getInitialOffset( )
    {
        if ( isNaN( interpolationOffsets_.first ) )
        {
            return 10.0 * interpolationStep_;
        }
        else
        {
            return interpolationOffsets_.first;
        }
    }

    //! Get final offset for the interpolated Spice ephemeris.
    /*!
     * @copybrief getFinalOffset
     * \remark If not defined by the user (i.e. is NaN), returns `10 * interpolationStep_`.
     * \return Final offset for the interpolated Spice ephemeris.
     */
    double getFinalOffset( )
    {
        if ( isNaN( interpolationOffsets_.second ) )
        {
            return 10.0 * interpolationStep_;
        }
        else
        {
            return interpolationOffsets_.second;
        }
    }
};

//! Create a `json` object from a shared pointer to a `SpiceSettings` object.
void to_json( nlohmann::json& jsonObject, const boost::shared_ptr< SpiceSettings >& spiceSettings );

//! Create a shared pointer to a `SpiceSettings` object from a `json` object.
void from_json( const nlohmann::json& jsonObject, boost::shared_ptr< SpiceSettings >& spiceSettings );


//! Load in Tudat the Spice kernels specified in \p spiceSettings.
/*!
 * @copybrief loadSpiceKernels
 * \remark Clears any Spice kernel loaded previously.
 * \remark If \p spiceSettings is `NULL`, no kernels are loaded.
 * \param spiceSettings The Spice settings containing the paths to the kernels to be loaded.
 */
void loadSpiceKernels( const boost::shared_ptr< SpiceSettings >& spiceSettings );

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_SPICE_H
