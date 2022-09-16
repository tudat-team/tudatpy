/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEMASSRATEMODELS_H
#define TUDAT_CREATEMASSRATEMODELS_H


#include <vector>
#include <string>

#include "tudat/astro/basic_astro/accelerationModelTypes.h"
#include "tudat/astro/basic_astro/accelerationModel.h"
#include "tudat/astro/basic_astro/massRateModel.h"
#include "tudat/astro/propulsion/massRateFromThrust.h"
#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{

namespace simulation_setup
{

//! Class for providing settings for a mass rate model.
/*!
 *  Class for providing settings for a mass rate model, that defines the models to be used to numerically propagate the
 *  mass of a body during a simulation. If any additional information (besides the type of the mass rate model) is required,
 *  these must be implemented in a derived class (dedicated for each mass rate model type).
 */
class MassRateModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param massRateType Type of mass rate model that is to be used.
     */
    MassRateModelSettings(
            const basic_astrodynamics::AvailableMassRateModels massRateType ):
        massRateType_( massRateType ){ }

    //! Destructor.
    virtual ~MassRateModelSettings( ){ }

    //! Type of mass rate model that is to be used.
    basic_astrodynamics::AvailableMassRateModels massRateType_;

};

//! Class defining the settings for a custom (i.e. predefined function of time) mass rate model.
class CustomMassRateSettings: public MassRateModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param massRateFunction Function returning the mass rate as a function of time.
     */
    CustomMassRateSettings(
            const std::function< double( const double ) > massRateFunction ):
        MassRateModelSettings( basic_astrodynamics::custom_mass_rate_model ),
    massRateFunction_( massRateFunction ){ }

    //! Destructor.
    ~CustomMassRateSettings( ){ }

    //! Function returning the mass rate as a function of time.
    std::function< double( const double ) > massRateFunction_;

};


//! Class defining the settings of a thrust model where the thrust is directly retrieved from (a) model(s) of an engine.
class FromThrustMassRateSettings: public MassRateModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param useAllThrustModels Boolean denoting whether all engines of the associated body are to be combined into a
     * single thrust model.
     * \param associatedThrustSource Name of engine model from which thrust is to be derived (must be empty if
     * useAllThrustModels is set to true)
     */
    FromThrustMassRateSettings(
            const bool useAllThrustModels = 1,
            const std::string& associatedThrustSource = "" ):
        MassRateModelSettings( basic_astrodynamics::from_thrust_mass_rate_model ),
        associatedThrustSource_( { associatedThrustSource } ), useAllThrustModels_( useAllThrustModels ){ }

    FromThrustMassRateSettings(
            const std::vector< std::string > associatedThrustSources ):
        MassRateModelSettings( basic_astrodynamics::from_thrust_mass_rate_model ),
    associatedThrustSource_( associatedThrustSources ), useAllThrustModels_( false ){ }

    //! Destructor
    ~FromThrustMassRateSettings( ){ }

    //! Name of engine model from which thrust is to be derived
    std::vector< std::string > associatedThrustSource_;

    //! Boolean denoting whether all engines of the associated body are to be combined into a single thrust model.
    bool useAllThrustModels_;

};


typedef std::map< std::string, std::vector< std::shared_ptr< MassRateModelSettings > > > SelectedMassRateModelMap;

inline std::shared_ptr< MassRateModelSettings > customMassRate(
        const std::function< double( const double ) > massRateFunction )
{
    return std::make_shared< CustomMassRateSettings >(massRateFunction );
}

inline std::shared_ptr< MassRateModelSettings > fromThrustMassRate(
        const bool useAllThrustModels = 1,
        const std::string& associatedThrustSource = "" )
{
    return std::make_shared< FromThrustMassRateSettings >(useAllThrustModels, associatedThrustSource );
}

//! Function to create a mass rate model
/*!
 * Function to create a mass rate model, from specific settings and the full set of environment models.
 * \param bodyWithMassRate Name of body for which a mass rate model is to be created.
 * \param massRateModelSettings Settings for the mass rate model that is to be created.
 * \param bodies List of pointers to body objects; defines the full simulation environment.
 * \param accelerationModels List of acceleration models that are used during numerical propagation (empty by default).
 * \return Mass rate model that is to be used during numerical propagation.
 */
std::shared_ptr< basic_astrodynamics::MassRateModel > createMassRateModel(
        const std::string& bodyWithMassRate,
        const std::shared_ptr< MassRateModelSettings > massRateModelSettings,
        const SystemOfBodies& bodies,
        const basic_astrodynamics::AccelerationMap& accelerationModels = basic_astrodynamics::AccelerationMap( ) );


//! Function to create a list of mass rate models for a list of bodies.
/*!
 * Function to create a list of mass rate models for a list of bodies, from specific settings and the full set of
 * environment models.
 * \param bodies List of pointers to body objects; defines the full simulation environment.
 * \param massRateModelSettings Settings for the mass rate models that are to be created (key is body id).
 * \param accelerationModels List of acceleration models that are used during numerical propagation (empty by default).
 * \return Mass rate models that are to be used during numerical propagation (key is body id)..
 */
basic_astrodynamics::MassRateModelMap createMassRateModelsMap(
        const SystemOfBodies& bodies,
        const SelectedMassRateModelMap& massRateModelSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModels = basic_astrodynamics::AccelerationMap( ) );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEMASSRATEMODELS_H
