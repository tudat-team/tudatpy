/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h"
#include "Tudat/Astrodynamics/Propulsion/massRateFromThrust.h"
#include "Tudat/SimulationSetup/body.h"

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
class CustomMassRateModelSettings: public MassRateModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param massRateFunction Function returning the mass rate as a function of time.
     */
    CustomMassRateModelSettings(
            const boost::function< double( const double ) > massRateFunction ):
        MassRateModelSettings( basic_astrodynamics::custom_mass_rate_model ),
    massRateFunction_( massRateFunction ){ }

    //! Destructor.
    ~CustomMassRateModelSettings( ){ }

    //! Function returning the mass rate as a function of time.
    boost::function< double( const double ) > massRateFunction_;

};


//! Class defining the settings of a thrust model where the thrust is directly retrieved from (a) model(s) of an engine.
class FromThrustMassModelSettings: public MassRateModelSettings
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param useAllThrustModels Boolean denoting whether all engines of the associated body are to be combined into a
     * single thrust model.
     * \param associatedThroustSource Name of engine model from which thrust is to be derived (must be empty if
     * useAllThrustModels is set to true)
     */
    FromThrustMassModelSettings(
            const bool useAllThrustModels = 1,
            const std::string& associatedThroustSource = "" ):
        MassRateModelSettings( basic_astrodynamics::from_thrust_mass_rate_model ),
    associatedThroustSource_( associatedThroustSource ), useAllThrustModels_( useAllThrustModels ){ }

    //! Destructor
    ~FromThrustMassModelSettings( ){ }

    //! Name of engine model from which thrust is to be derived
    std::string associatedThroustSource_;

    //! Boolean denoting whether all engines of the associated body are to be combined into a single thrust model.
    bool useAllThrustModels_;

};

//! Function to create a mass rate model
/*!
 * Function to create a mass rate model, from specific settings and the full set of environment models.
 * \param bodyWithMassRate Name of body for which a mass rate model is to be created.
 * \param massRateModelSettings Settings for the mass rate model that is to be created.
 * \param bodyMap List of pointers to body objects; defines the full simulation environment.
 * \param accelerationModels List of acceleration models that are used during numerical propagation (empty by default).
 * \return Mass rate model that is to be used during numerical propagation.
 */
boost::shared_ptr< basic_astrodynamics::MassRateModel > createMassRateModel(
        const std::string& bodyWithMassRate,
        const boost::shared_ptr< MassRateModelSettings > massRateModelSettings,
        const NamedBodyMap& bodyMap,
        const basic_astrodynamics::AccelerationMap& accelerationModels = basic_astrodynamics::AccelerationMap( ) );


//! Function to create a list of mass rate models for a list of bodies.
/*!
 * Function to create a list of mass rate models for a list of bodies, from specific settings and the full set of
 * environment models.
 * \param bodyMap List of pointers to body objects; defines the full simulation environment.
 * \param massRateModelSettings Settings for the mass rate models that are to be created (key is body id).
 * \param bodyMap List of pointers to body objects; defines the full simulation environment.
 * \param accelerationModels List of acceleration models that are used during numerical propagation (empty by default).
 * \return Mass rate models that are to be used during numerical propagation (key is body id)..
 */
std::map< std::string, std::vector< boost::shared_ptr< basic_astrodynamics::MassRateModel > > > createMassRateModelsMap(
        const NamedBodyMap& bodyMap,
        const std::map< std::string, std::vector< boost::shared_ptr< MassRateModelSettings > > >& massRateModelSettings,
        const basic_astrodynamics::AccelerationMap& accelerationModels = basic_astrodynamics::AccelerationMap( ) );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEMASSRATEMODELS_H
