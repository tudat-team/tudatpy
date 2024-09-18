/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MASSRATEMODEL_H
#define TUDAT_MASSRATEMODEL_H

#include <map>
#include <vector>
#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace basic_astrodynamics
{

//! Base class for determining the rate of change of a body's mass, to be used in numerical integration.
/*!
 *  Base class for determining the rate of change of a body's mass, to be used in numerical integration. Specific types
 *  of mass rate models are to be implemented in derived classes
 */
class MassRateModel
{
public:

    //! Constructor
    MassRateModel( ):
        currentTime_( TUDAT_NAN ), currentMassRate_( TUDAT_NAN ){ }

    //! Destructor
    virtual ~MassRateModel( ) { }

    //! Function to retrieve current mass rate
    /*!
     * Function to retrieve current mass rate, as set by last call to updateMembers (implemented in derived class)
     * \return Current mass rate
     */
    virtual double getMassRate( )
    {
        return currentMassRate_;
    }

    //! Update member variables used by the mass rate model, and internally compute mass rate.
    /*!
     * Updates member variables used by the mass rate model, and internally compute mass rate. In the case of mass rate
     * models containing varying parameters, function-pointers returning such a parameter (for instance
     * the Cartesian state of a body) will be set as a member variable.
     * This function evaluates such function-pointers and updates member variables to the 'current'
     * values of these parameters. Only these current values, not the function-pointers are then
     * used for the actual computation function.
     * \param currentTime Time at which acceleration model is to be updated.
     */
    virtual void updateMembers( const double currentTime = TUDAT_NAN ) = 0;

    //! Function to reset the current time
    /*!
     * Function to reset the current time of the acceleration model.
     * \param currentTime Current time (default NaN).
     */
    virtual void resetCurrentTime( )
    {
        currentTime_ = TUDAT_NAN;
    }

protected:

    //! Previous time to which mass rate model was updated.
    double currentTime_;

    //! Current mass rate, as set by last call to updateMembers (implemented in derived class)
    double currentMassRate_;

private:
};

//! Derived class for determining the rate of change of a body's mass, user-defined by a function pointer.
/*!
 *  Derived class for determining the rate of change of a body's mass, user-defined by a function pointer.
 *  This class can be used for any kind of mass rate model for which the user can define the dependency as a function
 *  of time.
 */
class CustomMassRateModel: public MassRateModel
{
public:

    //! Constructor.
    /*!
     * Constructor
     * \param massRateFunction Function returning mass rate as a function of time.
     */
    CustomMassRateModel(
            const std::function< double( const double ) > massRateFunction ):
    massRateFunction_( massRateFunction ){ }

    //! Destructor.
    ~CustomMassRateModel( ){ }

    //! Update member variables used by the mass rate model and compute the mass rate
    /*!
     * Update member variables used by the mass rate model and compute the mass rate
     * \param currentTime Time at which acceleration model is to be updated.
     */
    void updateMembers( const double currentTime = TUDAT_NAN )
    {
        // Check if update is needed.
        if( !( currentTime_ == currentTime ) )
        {
            currentMassRate_ = massRateFunction_( currentTime );
        }
    }

private:

    //! Function returning mass rate as a function of time.
    std::function< double( const double ) > massRateFunction_;

};

//! Typedef for the massrate model map.
typedef std::map< std::string, std::vector< std::shared_ptr< MassRateModel > > > MassRateModelMap;


} // namespace basic_astrodynamics

} // namespace tudat

#endif // TUDAT_MASSRATEMODEL_H
