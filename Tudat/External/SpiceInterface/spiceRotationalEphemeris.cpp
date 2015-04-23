#include <iostream>
#include <stdexcept>

#include <boost/lexical_cast.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/External/SpiceInterface/spiceInterface.h"
#include "Tudat/External/SpiceInterface/spiceRotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{


//! Function to calculate the rotation quaternion from target frame to original frame.
Eigen::Quaterniond SpiceRotationalEphemeris::getRotationToBaseFrame(
        const double ephemerisTime, const double julianDayAtEpoch )
{
    if( julianDayAtEpoch != basic_astrodynamics::JULIAN_DAY_ON_J2000 )
    {
        boost::throw_exception(
                    boost::enable_error_info(
                        std::runtime_error(
                            "Spice rotational ephemeris must take J2000 as"
                            " reference input time.") ) );
    }

    // Get rotational quaternion from spice wrapper function
    Eigen::Quaterniond rotationalState = spice_interface::computeRotationQuaternionBetweenFrames(
                targetFrameOrientation_, baseFrameOrientation_, ephemerisTime );

    return rotationalState;
}


}

}
