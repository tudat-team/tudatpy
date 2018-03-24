/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
#ifndef TUDAT_ENTRYGUIDANCE_H
#define TUDAT_ENTRYGUIDANCE_H

#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"

namespace tudat
{

namespace aerodynamics
{

//! Base class in which the aerodynamic angles (angle of attack, sideslip and bank angle) are computed. A derived class
//! implementing a specific guidance law is to be implemented by the user.
class AerodynamicGuidance
{
public:

    //! Constructor.
    AerodynamicGuidance( ):
        currentAngleOfAttack_( TUDAT_NAN ), currentAngleOfSideslip_( TUDAT_NAN ), currentBankAngle_( TUDAT_NAN ){ }

    //! Destructor
    virtual ~AerodynamicGuidance( ){ }

    //! Pure virtual function to update the guidance to the current time and state
    /*!
     *  Pure virtual function to update the guidance to the current time and state
     *  \param currentTime Time to which the guidance object is to be updated.
     */
    virtual void updateGuidance( const double currentTime ) = 0;

    //! Function to retrieve the current angle of attack, as set by last call to updateGuidance function
    /*!
     * Function to retrieve the current angle of attack, as set by last call to updateGuidance function
     * \return Current angle of attack, as set by last call to updateGuidance function
     */
    double getCurrentAngleOfAttack( )
    {
        return currentAngleOfAttack_;
    }

    //! Function to retrieve the current angle of sideslip, as set by last call to updateGuidance function
    /*!
     * Function to retrieve the current angle of sideslip, as set by last call to updateGuidance function
     * \return Current angle of sideslip, as set by last call to updateGuidance function
     */
    double getCurrentAngleOfSideslip( )
    {
        return currentAngleOfSideslip_;
    }

    //! Function to retrieve the current bank angle, as set by last call to updateGuidance function
    /*!
     * Function to retrieve the current bank angle, as set by last call to updateGuidance function
     * \return Current bank angle, as set by last call to updateGuidance function
     */
    double getCurrentBankAngle( )
    {
        return currentBankAngle_;
    }

protected:

    //! Current angle of attack, as set by last call to updateGuidance function
    double currentAngleOfAttack_;

    //! Current angle of sideslip, as set by last call to updateGuidance function
    double currentAngleOfSideslip_;

    //! Current bank angle, as set by last call to updateGuidance function
    double currentBankAngle_;

};

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_ENTRYGUIDANCE_H
