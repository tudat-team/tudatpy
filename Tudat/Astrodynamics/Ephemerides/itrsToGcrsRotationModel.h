/*    Copyright (c) 2010-2018, Delft University of Technology
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

#if USE_SOFA

#include <boost/bind.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Mathematics/Interpolators/interpolator.h"
#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"
#include "Tudat/Astrodynamics/EarthOrientation/earthOrientationCalculator.h"

namespace tudat
{

namespace ephemerides
{


//! Class for rotation from ITRS to GCRS, according to IERS 2010 models.
/*!
 *  Class for rotation from ITRS to GCRS, according to IERS 2010 models and rotation angle corrections. Angles may be provided by an interpolator
 *  to prevent this cllass becoming a computational bottleneck.
 */
class GcrsToItrsRotationModel: public RotationalEphemeris
{
public:

    //    //! Constructor taking interpolator providing the earth orientation angles.
    //    /*!
    //     *  Constructor taking interpolator providing the earth orientation angles.
    //     *  \param anglesInterpolator Interpolator providing the earth orientation angles (dependent variable) as a function of time (independent
    //     *  variable) The return vector of the interpolator provides the values for X-nutation correction, Y-nutation correction,
    //     *  CIO-locator, earth orientation angle, x-component polar motion, y-component polar motion.
    //     */
    //    GcrsToItrsRotationModel( const boost::shared_ptr< interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d > >
    //                        anglesInterpolator,
    //                        const basic_astrodynamics::TimeScales inputTimeScale  = basic_astrodynamics::tdb_scale ):
    //        RotationalEphemeris( "GCRS", "ITRS" ), anglesCalculator_( NULL ), inputTimeScale_( inputTimeScale )
    //    {
    //        using namespace interpolators;

    //        // Set function binding to interpolator.
    //        functionToGetRotationAngles = boost::bind(
    //                    static_cast< Eigen::Vector6d(
    //                        interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d >::* )( const double )>
    //                    ( &interpolators::OneDimensionalInterpolator< double, Eigen::Vector6d >::interpolate ), anglesInterpolator, _1 );
    //    }

    //! Constructor taking class calculating earth orientation angles directly
    /*!
     *  Constructor taking class calculating earth orientation angles directly
     *  \param anglesCalculator Class performing calculation to obtain earth orientation angle.
     *  \param timeScale Time scale in which input to this class (in getRotationToBaseFrame, getDerivativeOfRotationFromFrame) is provided,
     *  needed for correct input to EarthOrientationAnglesCalculator::getRotationAnglesFromItrsToGcrs.
     */
    GcrsToItrsRotationModel( const boost::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > anglesCalculator,
                             const basic_astrodynamics::TimeScales inputTimeScale  = basic_astrodynamics::tdb_scale ):
        RotationalEphemeris( "GCRS", "ITRS" ), anglesCalculator_( anglesCalculator ), inputTimeScale_( inputTimeScale )

    {
        functionToGetRotationAngles = boost::bind(
                    &earth_orientation::EarthOrientationAnglesCalculator::getRotationAnglesFromItrsToGcrs< double >,
                    anglesCalculator, _1, inputTimeScale );
    }

    //! Function to calculate the rotation quaternion from ITRS to GCRS
    /*!
     *  Function to calculate the rotation quaternion from ITRS to GCRS at specified time, to be implemented by derived class.
     *  \param ephemerisTime Time at which rotation is to be calculated.
     *  \return Rotation quaternion from ITRS to GCRS at specified time.
     */
    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime )
    {
        return earth_orientation::calculateRotationFromItrsToGcrs< double >(
                    anglesCalculator_->getRotationAnglesFromItrsToGcrs< double >( ephemerisTime, inputTimeScale_ ),
                    ephemerisTime );
    }

    Eigen::Quaterniond getRotationToBaseFrameFromExtendedTime( const Time ephemerisTime )
    {
        return earth_orientation::calculateRotationFromItrsToGcrs< Time >(
                    anglesCalculator_->getRotationAnglesFromItrsToGcrs< Time >( ephemerisTime, inputTimeScale_ ),
                    ephemerisTime );
    }



    Eigen::Quaterniond getRotationToTargetFrame( const double ephemerisTime )
    {
        return getRotationToBaseFrame( ephemerisTime ).inverse( );
    }

    Eigen::Quaterniond getRotationToTargetFrameFromExtendedTime( const Time ephemerisTime )
    {
        return getRotationToBaseFrameFromExtendedTime( ephemerisTime ).inverse( );
    }



    //! Function to calculate the derivative of the rotation matrix from ITRS to GCRS
    /*!
     *  Function to calculate the derivative of the rotation matrix from ITRS to GCRS at specified time,
     *  \param ephemerisTime Time at which derivative of rotation is to be calculated.
     *  \return Derivative of rotation from ITRS to GCRS at specified time.
     */
    Eigen::Matrix3d getDerivativeOfRotationToBaseFrame( const double ephemerisTime )
    {
        return earth_orientation::calculateRotationRateFromItrsToGcrs< double >( functionToGetRotationAngles( ephemerisTime ),
                                                                                 ephemerisTime );
    }

    Eigen::Matrix3d getDerivativeOfRotationToTargetFrame( const double ephemerisTime )
    {
        return getDerivativeOfRotationToBaseFrame( ephemerisTime ).transpose( );
    }

    boost::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > getAnglesCalculator( )
    {
        return anglesCalculator_;
    }

    basic_astrodynamics::TimeScales getInputTimeScale( )
    {
        return inputTimeScale_;
    }


private:

    //! Function providing the earth orientation angles as a function of time
    /*!
     * Function providing the earth orientation angles as a function of time.
     *  The return vector of the interpolator provides the values for X-nutation correction, Y-nutation correction,
     *  CIO-locator, earth orientation angle, x-component polar motion, y-component polar motion.
     */
    boost::function< std::pair< Eigen::Vector5d, double >( const double& ) > functionToGetRotationAngles;

    boost::shared_ptr< earth_orientation::EarthOrientationAnglesCalculator > anglesCalculator_;

    basic_astrodynamics::TimeScales inputTimeScale_;
};

}

}

#endif

#endif // TUDAT_GCRSTOITRSROTATIONMODEL_H
