#ifndef POLETIDEDEFORMATION_H
#define POLETIDEDEFORMATION_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"

#include "Astrodynamics/EarthOrientation/polarMotionCalculator.h"
#include "External/SofaInterface/sofaTimeConversions.h"
#include "Astrodynamics/GroundStations/groundStation.h"
#include "Astrodynamics/SiteDisplacements/bodyDeformationModel.h"

namespace tudat
{

namespace site_displacements
{

//! This class is used to calculate the displacement of a site on Earth due to pole tide deformation.
/*!
 *  This class is used to calculate the displacement of a site on Earth due to pole tide deformation,
 *  as described in Section 7.1.4 of the IERS 2010 Conventions.
 */
class PoleTideDeformation: public BodyDeformationModel
{
public:

    //! Constructor creating pole tide deformation object from polar motion model.
    /*!
     *  Constructor creating pole tide deformation object from polar motion model.
     *  \param polarMotionCalculator Object capable of calculating pole position in ITRS at a given time instant.
     */
    PoleTideDeformation( const  boost::function< Eigen::Vector2d( const double ) > polarMotionFunction ):
        polarMotionFunction_( polarMotionFunction ){ }

    //! Function to calculate site displacement due to pole tide deformation.
    /*!
     *  Function to calculate site displacement due to pole tide deformation.
     *  \param terrestrialTime Time (TT scale) in seconds since J2000.
     *  \param nominalSiteState Data structure containing the nominal state of the site.
     *  \return Displacement of site in ITRS at given site due to pole tide.
     */
    Eigen::Vector3d calculateSiteDisplacement( const double terrestrialTime,
                                               const boost::shared_ptr< NominalGroundStationState > nominalSiteState );

    Eigen::Vector3d calculateSiteDisplacement(
            const double time,
            const boost::shared_ptr< NominalGroundStationState > siteState,
            const basic_mathematics::Vector6d doodsonArguments )
    {
        return calculateSiteDisplacement( time, siteState );
    }

private:

    //! Pointer to object to calculate the polar motion at a given time.
    /*!
     *  Pointer to object to calculate the polar motion at a given time.
     */
    boost::function< Eigen::Vector2d( const double ) > polarMotionFunction_;
};


}

}
#endif // POLETIDEDEFORMATION_H
