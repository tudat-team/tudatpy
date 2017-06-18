#ifndef COMPOSITEROTATIONALEPHEMERIS_H
#define COMPOSITEROTATIONALEPHEMERIS_H

#include <vector>

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/Ephemerides/rotationalEphemeris.h"

namespace tudat
{

namespace ephemerides
{

//! Class for combination of rotational ephemerides.
class CompositeRotationalEphemeris : public RotationalEphemeris
{
public:
    //! Constructor, taking list of rotations to apply.
    /*!
     *  Constructor, taking list of rotations to apply (index 0 applied last to vector)
     *  \param subRotations List of rotations.
     */
    CompositeRotationalEphemeris( std::vector< boost::shared_ptr< RotationalEphemeris > > subRotations );

    //! Destructor
    /*!
     *  Destructor.
     */
    ~CompositeRotationalEphemeris( ) { }

    Eigen::Quaterniond getRotationToBaseFrame( const double ephemerisTime );

    Eigen::Matrix3d getDerivativeOfRotationFromFrame( const double ephemerisTime );

    std::vector< boost::shared_ptr< RotationalEphemeris > >  getSubRotations( )
    {
        return subRotations_;
    }

private:
    std::vector< boost::shared_ptr< RotationalEphemeris > > subRotations_;
};

}

}

#endif // COMPOSITEROTATIONALEPHEMERIS_H
