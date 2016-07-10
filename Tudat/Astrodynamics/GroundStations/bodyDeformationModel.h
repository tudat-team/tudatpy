#ifndef BODYDEFORMATIONMODEL_H
#define BODYDEFORMATIONMODEL_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/GroundStations/nominalGroundStationState.h"

namespace tudat
{

namespace site_displacements
{

//! Base class for body deformation models, such as tides, ocean tides, etc.
/*!
 *  Base class for body deformation models, such as tides, ocean tides, etc. Each derived class implements a single type of deformation.
 *  A list of objects of this type are provided to the CelestialBody
 */
class BodyDeformationModel
{
public:

    BodyDeformationModel( ){ }

    virtual ~BodyDeformationModel( ){ }

    virtual Eigen::Vector3d calculateSiteDisplacement(
            const double time,
            const boost::shared_ptr< ground_stations::NominalGroundStationState > siteState ) = 0;
protected:
};

std::map< double, Eigen::Vector3d > calculateSiteDisplacements(
        const std::vector< double >& times,
        const boost::shared_ptr< ground_stations::NominalGroundStationState > nominalSiteState,
        const boost::shared_ptr< BodyDeformationModel > deformationModel );


}

}

#endif // BODYDEFORMATIONMODEL_H
