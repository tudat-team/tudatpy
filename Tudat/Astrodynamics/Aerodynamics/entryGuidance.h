#ifndef TUDAT_ENTRYGUIDANCE_H
#define TUDAT_ENTRYGUIDANCE_H

#include <boost/bind.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/ReferenceFrames/aerodynamicAngleCalculator.h"

namespace tudat
{

namespace aerodynamics
{

//! Base class in which the aerodynamic angles (angle of attack, sideslip and bank angle) are computed. A derived class
//! implementing a specific guidance law is to be implemented by the user.
class EntryGuidance
{
public:

    //! Constructor.
    EntryGuidance( ):
        currentAngleOfAttack_( TUDAT_NAN ), currentAngleOfSideslip_( TUDAT_NAN ), currentBankAngle_( TUDAT_NAN ){ }

    //! Destructor
    virtual ~EntryGuidance( ){ }

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

//! Function that must be called to link the EntryGuidance object to the simulation
/*!
 * Function that must be called to link the EntryGuidance object to the simulation
 * \param entryGuidance Object computing the current aerodynamic angles.
 * \param angleCalculator Object that handles all aerodynamic angles in the numerical propagation
 */
void setGuidanceAnglesFunctions(
        const boost::shared_ptr< EntryGuidance > entryGuidance,
        const boost::shared_ptr< reference_frames::AerodynamicAngleCalculator > angleCalculator );

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_ENTRYGUIDANCE_H
