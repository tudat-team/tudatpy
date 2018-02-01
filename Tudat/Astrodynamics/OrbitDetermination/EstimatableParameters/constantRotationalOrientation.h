/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONSTANTROTATIONALORIENTATION_H
#define TUDAT_CONSTANTROTATIONALORIENTATION_H


#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's constant pole position (right ascension and declination of north pole).
/*!
 *  Interface class for estimation of a body's constant pole position (right ascension and declination of north pole).
 *  Interfaces the estimation with the right ascension and declination Euler angle members of a SimpleRotationalEphemeris
 *  object
 */
class ConstantRotationalOrientation: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param rotationModel SimpleRotationalEphemeris object of which pole position is a property
     *  \param associatedBody Name of body of which parameter is a property.
     */
    ConstantRotationalOrientation(
            const boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( rotation_pole_position, associatedBody ),
        rotationModel_( rotationModel ) { }

    //! Destructor
    ~ConstantRotationalOrientation( ) { }

    //! Get value of pole right ascension and declination (in that order)
    /*!
     *  Get value of pole right ascension and declination (in that order)
     *  \return Right ascension and declination (in that order)
     */
    Eigen::VectorXd getParameterValue( )
    {
        return rotationModel_->getInitialEulerAngles( ).segment( 0, 2 );
    }

    //! Reset value of pole right ascension and declination (in that order)
    /*!
     *  Reset value of pole right ascension and declination (in that order)
     *  \param parameterValue New right ascension and declination (in that order)
     */
    void setParameterValue( const Eigen::VectorXd parameterValue )
    {
        rotationModel_->resetInitialPoleRightAscensionAndDeclination(
                    parameterValue.x( ), parameterValue.y( ) );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value, 2 for this parameter
     */
    int getParameterSize( )
    {
        return 2;
    }

protected:

private:

    //! SimpleRotationalEphemeris object of which rotation rate parameter is a property
    boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel_;
};

} // namespace estimatable_parameters

} // namespace tudat

#endif // TUDAT_CONSTANTROTATIONALORIENTATION_H
