/*    Copyright (c) 2010-2016, Delft University of Technology
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

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/SystemModels/engineModel.h"

namespace tudat
{

namespace system_models
{

//! Wrapper class that contains the relevant hardware systems of a vehicle.
/*!
 *  Wrapper class that contains the relevant hardware systems of a vehicle. Not all member objects need to be set; NULL
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

    //! Function to retrieve the engine models
    /*!
     * Function to retrieve the engine models
     * \return Named list of engine models in the vehicle
     */
    std::map< std::string, boost::shared_ptr< EngineModel > > getEngineModels( )
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
            const boost::shared_ptr< EngineModel > engineModel, const std::string engineName = "" )
    {
        // Check if engine with this name already exists.
        if( engineModels_.count( engineName ) )
        {
            std::cerr<<"Warning, engine model of name "<<engineModel<<" already exists, overriding old model"<<std::endl;
        }

        engineModels_[ engineName ] = engineModel;
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

private:

    //! Named list of engine models in the vehicle
    std::map< std::string, boost::shared_ptr< EngineModel > > engineModels_;

    //! Total dry mass of the vehicle
    double dryMass_;

};


} // namespace system_models

} // namespace tudat

#endif // TUDAT_VEHICLESYSTEMS_H
