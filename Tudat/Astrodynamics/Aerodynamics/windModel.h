/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_WIND_MODEL_H
#define TUDAT_WIND_MODEL_H

#include <Eigen/Core>

#include <boost/shared_ptr.hpp>
#include <boost/function.hpp>

namespace tudat
{
namespace aerodynamics
{

class WindModel
{
public:
    WindModel( ){ }

    ~WindModel( ){ }

    virtual Eigen::Vector3d getCurrentWindVelocity(
            const double currentAltitude ,
            const double currentLongitude ,
            const double currentLatitude,
            const double currentTime ) = 0;
};

class CustomWindModel: public WindModel
{
public:
    CustomWindModel(
            const boost::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction ):
    windFunction_( windFunction ){ }

    ~CustomWindModel( ){ }

    Eigen::Vector3d getCurrentWindVelocity(
            const double currentAltitude ,
            const double currentLongitude ,
            const double currentLatitude,
            const double currentTime )
    {
        return windFunction_( currentAltitude, currentLongitude, currentLatitude, currentTime );
    }

private:
    boost::function< Eigen::Vector3d( const double, const double, const double, const double ) > windFunction_;
};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_WIND_MODEL_H
