/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_GCRSTOITRSROTATIONMODEL_H
#define TUDAT_GCRSTOITRSROTATIONMODEL_H

// #if TUDAT_BUILD_WITH_SOFA_INTERFACE

#include "tudat/math/basic/linearAlgebra.h"
#include "tudat/math/interpolators/interpolator.h"
#include "tudat/astro/ephemerides/rotationalEphemeris.h"
#include "tudat/astro/earth_orientation/earthOrientationCalculator.h"
#include "tudat/interface/spice/spiceInterface.h"

namespace tudat
{

namespace ephemerides
{

//! Intermediate frames in the IERS Earth-orientation rotation chain.
enum class EarthOrientationIntermediateFrame { itrs, tirs, cirs, icrs };

//! Class for rotation from ITRS to GCRS, according to IERS 2010 models.
/*!
 *  Class for rotation from ITRS to GCRS, according to IERS 2010 models and rotation angle corrections. Angles may be provided by an
 * interpolator to prevent this cllass becoming a computational bottleneck.
 */
class GcrsToItrsRotationModel : public RotationalEphemeris
{
public:
    //! Constructor taking class calculating earth orientation angles directly
    /*!
     *  Constructor taking class calculating earth orientation angles directly
     *  \param anglesCalculator Class performing calculation to obtain earth orientation angle.
     *  \param timeScale Time scale in which input to this class (in getRotationToBaseFrame, getDerivativeOfRotationToBaseFrame) is
     * provided, needed for correct input to EarthOrientationAnglesCalculator::getRotationAnglesFromItrsToGcrs.
     *  \param baseFrame Base axes: GCRS, J2000 (SPICE/ICRF convention), ICRF,
     *  ECLIPJ2000 (SPICE ecliptic axes), or EME2000 (dynamical mean equator/equinox at J2000).
     *  J2000 and ICRF use the GCRS axis orientation. This does not equate frame origins,
     *  time scales, or relativistic coordinate scaling. EME2000 applies the SOFA frame bias.
     */
    GcrsToItrsRotationModel( const std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > anglesCalculator,
                             const basic_astrodynamics::TimeScales inputTimeScale = basic_astrodynamics::tdb_scale,
                             const std::string& baseFrame = "GCRS" ):
        RotationalEphemeris( baseFrame, "ITRS" ), anglesCalculator_( anglesCalculator ), inputTimeScale_( inputTimeScale ),
        rotationFromGcrsToBaseFrame_( Eigen::Matrix3d::Identity( ) )

    {
        // Modern SPICE J2000 ephemerides already have ICRF axes. Applying
        // SOFA's dynamical mean-J2000 bias to these vectors adds an unwanted rotation.
        // The previous J2000 bias behavior is now explicitly named EME2000.
        if( baseFrame == "EME2000" )
        {
            rotationFromGcrsToBaseFrame_ = sofa_interface::getFrameBias(
                    0.0, anglesCalculator->getPrecessionNutationCalculator( )->getPrecessionNutationTheory( ) );
        }
        else if( baseFrame == "ECLIPJ2000" )
        {
            rotationFromGcrsToBaseFrame_ = spice_interface::getRotationFromJ2000ToEclipJ2000( );
        }
        else if( baseFrame != "GCRS" && baseFrame != "J2000" && baseFrame != "ICRF" )
        {
            throw std::runtime_error( "Error in GCRS<->ITRS model, base frame not recognized: " + baseFrame +
                                      ". Supported frames are GCRS, J2000, ICRF, ECLIPJ2000 and EME2000." );
        }
    }

    //! Function to calculate the rotation quaternion from ITRS to base frame
    /*!
     *  Function to calculate the rotation quaternion from ITRS to base frame at specified time.
     *  \param ephemerisTime Time at which rotation is to be calculated.
     *  \return Rotation quaternion from ITRS to base frame at specified time.
     */
    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime ) override
    {
        return Eigen::Quaterniond( rotationFromGcrsToBaseFrame_ ) *
                getRotationFromOrientation( anglesCalculator_->getRotationAnglesAndTimesFromItrsToGcrs( ephemerisTime, inputTimeScale_ ) );
    }

    //! Function to calculate the rotation quaternion from ITRS to base frame
    /*!
     *  Function to calculate the rotation quaternion from ITRS to base frame at specified time, in extended (e.g. Time class)
     *  format.
     *  \param ephemerisTime Time at which rotation is to be calculated.
     *  \return Rotation quaternion from ITRS to base frame at specified time.
     */
    Eigen::Quaterniond getRotationToBaseFrameFromExtendedTime( const Time ephemerisTime ) override
    {
        return Eigen::Quaterniond( rotationFromGcrsToBaseFrame_ ) *
                getRotationFromOrientation( anglesCalculator_->getRotationAnglesAndTimesFromItrsToGcrs( ephemerisTime, inputTimeScale_ ) );
    }

    //! Function to calculate the rotation quaternion from base frame to ITRS
    /*!
     *  Function to calculate the rotation quaternion from base frame to ITRS at specified time.
     *  \param ephemerisTime Time at which rotation is to be calculated.
     *  \return Rotation quaternion from base frame to ITRS at specified time.
     */
    Eigen::Quaterniond getRotationToTargetFrame( const double ephemerisTime ) override
    {
        return getRotationToBaseFrame( ephemerisTime ).inverse( );
    }

    //! Function to calculate the rotation quaternion from base frame to ITRS
    /*!
     *  Function to calculate the rotation quaternion from base frame to ITRS at specified time, in extended (e.g. Time class)
     *  formats
     *  \param ephemerisTime Time at which rotation is to be calculated.
     *  \return Rotation quaternion from base frame to ITRS at specified time.
     */
    Eigen::Quaterniond getRotationToTargetFrameFromExtendedTime( const Time ephemerisTime ) override
    {
        return getRotationToBaseFrameFromExtendedTime( ephemerisTime ).inverse( );
    }

    //! Function to calculate the derivative of the rotation matrix from ITRS to base frame
    /*!
     *  Function to calculate the derivative of the rotation matrix from ITRS to base frame at specified time,
     *  \param ephemerisTime Time at which derivative of rotation is to be calculated.
     *  \return Derivative of rotation from ITRS to base frame at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double ephemerisTime ) override
    {
        return rotationFromGcrsToBaseFrame_ *
                earth_orientation::calculateRotationRateFromItrsToGcrs(
                        anglesCalculator_->getRotationAnglesAndRatesFromItrsToGcrs( ephemerisTime, inputTimeScale_ ) );
    }

    Eigen::Matrix3d getDerivativeOfRotationToBaseFrameFromExtendedTime( const Time ephemerisTime ) override
    {
        return rotationFromGcrsToBaseFrame_ *
                earth_orientation::calculateRotationRateFromItrsToGcrs(
                        anglesCalculator_->getRotationAnglesAndRatesFromItrsToGcrs( ephemerisTime, inputTimeScale_ ) );
    }

    Eigen::Matrix3d getDerivativeOfRotationToTargetFrameFromExtendedTime( const Time ephemerisTime ) override
    {
        return getDerivativeOfRotationToBaseFrameFromExtendedTime( ephemerisTime ).transpose( );
    }

    //! Function to calculate the derivative of the rotation matrix from base frame to ITRS
    /*!
     *  Function to calculate the derivative of the rotation matrix from base frame to ITRS at specified time,
     *  \param ephemerisTime Time at which derivative of rotation is to be calculated.
     *  \return Derivative of rotation from base frame to ITRS at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double ephemerisTime ) override
    {
        return getDerivativeOfRotationToBaseFrame( ephemerisTime ).transpose( );
    }

    //! Function to calculate the rotation quaternion between intermediate Earth-orientation frames.
    /*!
     *  Function to calculate the rotation quaternion from an original intermediate frame to a target intermediate frame at a
     *  specified time. The ICRS enum value denotes the GCRS-side endpoint of the IERS 2010 ITRS->TIRS->CIRS->GCRS rotation
     *  sequence used by this model.
     *  \param originalFrame Frame from which the rotation is calculated.
     *  \param targetFrame Frame to which the rotation is calculated.
     *  \param ephemerisTime Time at which rotation is to be calculated.
     *  \return Rotation quaternion from originalFrame to targetFrame at specified time.
     */
    Eigen::Quaterniond getRotationBetweenIntermediateFrames( const EarthOrientationIntermediateFrame originalFrame,
                                                             const EarthOrientationIntermediateFrame targetFrame,
                                                             const double ephemerisTime )
    {
        const auto orientation = anglesCalculator_->getRotationAnglesAndTimesFromItrsToGcrs( ephemerisTime, inputTimeScale_ );

        auto getFrameIndex = []( const EarthOrientationIntermediateFrame frame ) {
            switch( frame )
            {
                case EarthOrientationIntermediateFrame::itrs:
                    return 0;
                case EarthOrientationIntermediateFrame::tirs:
                    return 1;
                case EarthOrientationIntermediateFrame::cirs:
                    return 2;
                case EarthOrientationIntermediateFrame::icrs:
                    return 3;
                default:
                    throw std::runtime_error( "Error when getting Earth rotation between intermediate frames, frame not recognized" );
            }
        };

        auto getForwardRotationForChainLink = [ & ]( const int lowerFrameIndex ) {
            switch( lowerFrameIndex )
            {
                case 0:
                    return earth_orientation::calculateRotationFromItrsToTirs(
                            orientation.angles[ 3 ],
                            orientation.angles[ 4 ],
                            earth_orientation::getApproximateTioLocator( static_cast< double >( orientation.tt ) ) );
                case 1:
                    return earth_orientation::calculateRotationFromTirsToCirs(
                            sofa_interface::calculateEarthRotationAngleTemplated< Time >( orientation.ut1 ) );
                case 2:
                    return earth_orientation::calculateRotationFromCirsToGcrs(
                            orientation.angles[ 0 ], orientation.angles[ 1 ], orientation.angles[ 2 ] );
                default:
                    throw std::runtime_error( "Error when getting Earth rotation between intermediate frames, chain link not recognized" );
            }
        };

        int originalFrameIndex = getFrameIndex( originalFrame );
        int targetFrameIndex = getFrameIndex( targetFrame );
        Eigen::Quaterniond rotationBetweenFrames = Eigen::Quaterniond::Identity( );

        if( originalFrameIndex < targetFrameIndex )
        {
            for( int currentFrameIndex = originalFrameIndex; currentFrameIndex < targetFrameIndex; currentFrameIndex++ )
            {
                rotationBetweenFrames = getForwardRotationForChainLink( currentFrameIndex ) * rotationBetweenFrames;
            }
        }
        else
        {
            for( int currentFrameIndex = originalFrameIndex - 1; currentFrameIndex >= targetFrameIndex; currentFrameIndex-- )
            {
                rotationBetweenFrames = getForwardRotationForChainLink( currentFrameIndex ).inverse( ) * rotationBetweenFrames;
            }
        }

        return rotationBetweenFrames;
    }

    //! Function to retrieve object responsible for computing the various rotation angles (precession, nutation, polar motion, etc.)
    /*!
     * Function to retrieve object responsible for computing the various rotation angles (precession, nutation, polar motion, etc.)
     * \return object responsible for computing the various rotation angles (precession, nutation, polar motion, etc.)
     */
    std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > getAnglesCalculator( )
    {
        return anglesCalculator_;
    }

    //! Function to retrieve time scale in which the input time for class functions are interpreted
    /*!
     * Function to retrieve time scale in which the input time for class functions are interpreted
     * \return Time scale in which the input time for class functions are interpreted
     */
    basic_astrodynamics::TimeScales getInputTimeScale( )
    {
        return inputTimeScale_;
    }

private:
    static Eigen::Quaterniond getRotationFromOrientation( const earth_orientation::EarthOrientationAngles& orientation )
    {
        return earth_orientation::calculateRotationFromItrsToGcrs< Time >( std::make_pair( orientation.angles, orientation.ut1 ),
                                                                           static_cast< double >( orientation.tt ) );
    }

    //! Object responsible for computing the various rotation angles (precession, nutation, polar motion, etc.)
    /*!
     *  Object responsible for computing the various rotation angles (precession, nutation, polar motion, etc.), as per IERS 2010
     *  conventions.
     */
    std::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > anglesCalculator_;

    //! Time scale in which the input time for class functions are interpreted
    basic_astrodynamics::TimeScales inputTimeScale_;

    //! Constant rotation from GCRS axes to the selected base-frame axes.
    /*!
     * Identity for GCRS, SPICE J2000 and ICRF; the SPICE ecliptic rotation for
     * ECLIPJ2000; or the SOFA dynamical mean-J2000 frame bias for EME2000.
     */
    Eigen::Matrix3d rotationFromGcrsToBaseFrame_;
};

}  // namespace ephemerides

}  // namespace tudat

// #endif

#endif  // TUDAT_GCRSTOITRSROTATIONMODEL_H
