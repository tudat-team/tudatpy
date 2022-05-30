/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_VEHICLESYSTEMS_H
#define TUDAT_VEHICLESYSTEMS_H

#include <map>
#include <iostream>

#include <memory>

#include "tudat/astro/system_models/engineModel.h"

namespace tudat
{

namespace system_models
{

//! Wrapper class that contains the relevant hardware systems of a vehicle.
/*!
 *  Wrapper class that contains the relevant hardware systems of a vehicle. Not all member objects need to be set; nullptr
 *  member objects denote that the vehicle does not contain the associated hardware.
 */
class VehicleSystems
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dryMass Total dry mass of the vehicle (not defined; NaN by default).
     */
    VehicleSystems( const double dryMass = TUDAT_NAN ):
        dryMass_( dryMass ){ }

    //! Destructor
    ~VehicleSystems( ){ }

    //! Function to retrieve the engine models
    /*!
     * Function to retrieve the engine models
     * \return Named list of engine models in the vehicle
     */
    std::map< std::string, std::shared_ptr< EngineModel > > getEngineModels( )
    {
        return engineModels_;
    }

    //! Function to set a single engine in the vehicle
    /*!
     * Function to set a single engine in the vehicle. Each engine can be identified by a string. If only a single
     * engine is set, the default (empty string) can be used.
     * \param engineModel Model of engine that is to be set
     * \param engineName Reference id of the engine that is to be set.
     */
    void setEngineModel(
            const std::shared_ptr< EngineModel > engineModel )
    {
        // Check if engine with this name already exists.
        if( engineModels_.count( engineModel->getEngineName( ) ) )
        {
            std::cerr << "Warning, engine model of name " << engineModel->getEngineName( ) << " already exists, overriding old model" << std::endl;
        }

        engineModels_[ engineModel->getEngineName( ) ] = engineModel;
    }

    //! Function to retrieve the total dry mass of the vehicle
    /*!
     * Function to retrieve the  total dry mass of the vehicle
     * \return Total dry mass of the vehicle
     */
    double getDryMass( )
    {
        return dryMass_;
    }

    //! Function to set the current deflection of a single control surface
    /*!
     * Function to set the current deflection of a single control surface
     * \param controlSurfaceId Name of control surface for which deflection is to be set
     * \param deflectionAngle Current deflection of control surface that is to be set.
     */
    void setCurrentControlSurfaceDeflection(
            const std::string& controlSurfaceId, const double deflectionAngle )
    {
        currentControlSurfaceDeflections_[ controlSurfaceId ] =  deflectionAngle;
    }

    //! Function to retrieve the current deflection of a single control surface
    /*!
     * Function to retrieve the current deflection of a single control surface
     * \param controlSurfaceId Name of control surface for which deflection is to be set
     * \return Current deflection of control surface that is requested.
     */
    double getCurrentControlSurfaceDeflection(
            const std::string& controlSurfaceId )
    {
        return currentControlSurfaceDeflections_.at( controlSurfaceId );
    }

    //! Function to (re)set the vehicle nose radius
    /*!
     * Function to (re)set the vehicle nose radius
     * \param noseRadius The  vehicle nose radius that is to be set
     */
    void setNoseRadius( const double noseRadius )
    {
        noseRadius_ = noseRadius;
    }

    //! Function to retrieve the vehicle nose radius
    /*!
     * Function to retrieve the vehicle nose radius
     * \return The vehicle nose radius
     */
    double getNoseRadius( )
    {
        return noseRadius_;
    }

    //! Function to (re)set the vehicle wall emissivity
    /*!
     * Function to (re)set the vehicle wall emissivity
     * \param wallEmissivity The vehicle wall emissivity that is to be set
     */
    void setWallEmissivity( const double wallEmissivity )
    {
        wallEmissivity_ = wallEmissivity;
    }

    //! Function to retrieve the vehicle wall emissivity
    /*!
     * Function to retrieve the vehicle wall emissivity
     * \return The vehicle wall emissivity
     */
    double getWallEmissivity( )
    {
        return wallEmissivity_;
    }

private:

    //! Named list of engine models in the vehicle
    std::map< std::string, std::shared_ptr< EngineModel > > engineModels_;

    //! Total dry mass of the vehicle
    double dryMass_;

    //! List if current control surface deflections (with key the control surface id).
    std::map< std::string, double > currentControlSurfaceDeflections_;

    //! Nose radius of the vehicle (used for heating computations)
    double noseRadius_;

    //! Wall emissivity of the vehicle (used for heating computations)
    double wallEmissivity_;

};

} // namespace system_models

} // namespace tudat

#endif // TUDAT_VEHICLESYSTEMS_H
