/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */


#ifndef TUDAT_GROUNDSTATIONPOSITION_H
#define TUDAT_GROUNDSTATIONPOSITION_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/GroundStations/groundStationState.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"


namespace tudat
{

namespace estimatable_parameters
{


//! Class for estimating constant position of a ground station
/*!
 *  Class for estimating constant position of a ground station. Specific parameter is Cartesian x,y and z component of station
 *  in body-fixed frame.
 */
class GroundStationPosition: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param groundStationState Object that represents state of ground station
     *  \param associatedBody Body on which station is located
     *  \param associatedStation Name of station
     */
    GroundStationPosition( const boost::shared_ptr<ground_stations:: GroundStationState > groundStationState,
                           const std::string& associatedBody,
                           const std::string& associatedStation ):
        EstimatableParameter< Eigen::VectorXd  >( ground_station_position, associatedBody, associatedStation ),
        groundStationState_( groundStationState ){ }

    //! Destructor    
    ~GroundStationPosition( ) { }

    //! Get values (Cartesian x,y,z position) of ground station position
    /*!
     *  Get values (Cartesian x,y,z position) of ground station position
     *  \return Values (Cartesian x,y,z position) of ground station position
     */
    Eigen::VectorXd  getParameterValue( )
    {
        return groundStationState_->getNominalCartesianPosition( );
    }

    //! Reset values (Cartesian x,y,z position) of ground station position
    /*!
     *  Reset values (Cartesian x,y,z position) of ground station position
     *  \param parameterValue New values (Cartesian x,y,z position) of ground station position
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        groundStationState_->resetGroundStationPositionAtEpoch( parameterValue );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value, 3 for this parameter
     */
    int getParameterSize( )
    {
        return 3;
    }

protected:

private:

    //! Object that represents state of ground stations
    boost::shared_ptr< ground_stations::GroundStationState > groundStationState_;
};

}

}

#endif // TUDAT_GROUNDSTATIONPOSITION_H
