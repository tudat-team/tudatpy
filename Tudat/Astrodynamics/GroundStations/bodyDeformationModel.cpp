#include "Tudat/Astrodynamics/GroundStations/bodyDeformationModel.h"

namespace tudat
{

namespace site_displacements
{

std::map< double, Eigen::Vector3d > calculateSiteDisplacements(
        const std::vector< double >& times,
        const boost::shared_ptr< ground_stations::NominalGroundStationState > nominalSiteState,
        const boost::shared_ptr< BodyDeformationModel > deformationModel )
{
    std::map< double, Eigen::Vector3d > displacementMap;
    for( unsigned int i = 0; i < times.size( ); i++ )
    {
        displacementMap[ times[ i ] ] = deformationModel->calculateSiteDisplacement(
                    times[ i ], nominalSiteState );
    }
    return displacementMap;
}

}

}
