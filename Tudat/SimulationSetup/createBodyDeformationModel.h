#ifndef CREATEBODYDEFORMATIONMODEL_H
#define CREATEBODYDEFORMATIONMODEL_H

#include <string>
#include <vector>

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Gravitation/gravityFieldVariations.h"
#include "Tudat/Astrodynamics/GroundStations/bodyDeformationModel.h"
#include "Tudat/Astrodynamics/GroundStations/nominalGroundStationState.h"

namespace tudat
{

namespace simulation_setup
{

class BodyDeformationSettings
{
public:
    BodyDeformationSettings( gravitation::BodyDeformationTypes bodyDeformationType ):bodyDeformationType_( bodyDeformationType ){ }

    virtual ~BodyDeformationSettings( ){ }

    gravitation::BodyDeformationTypes getBodyDeformationType( ){ return bodyDeformationType_; }

protected:
    gravitation::BodyDeformationTypes bodyDeformationType_;
};

class BasicSolidBodyDeformationSettings: public BodyDeformationSettings
{
public:
    BasicSolidBodyDeformationSettings( const std::vector< std::string > deformingBodies,
                                       const std::vector< int > maximumDeformationOrders,
                                       const std::vector< double > displacementLoveNumbers,
                                       const std::vector< double > displacementShidaNumbers,
                                       const double bodyReferenceRadius,
                                       const std::vector< Eigen::Vector3d > meanVectorsToDeformingBodies =
            std::vector< Eigen::Vector3d >( ) ):
        BodyDeformationSettings( gravitation::basic_solid_body ), deformingBodies_( deformingBodies ),
        maximumDeformationOrders_( maximumDeformationOrders ), displacementLoveNumbers_( displacementLoveNumbers ),
        displacementShidaNumbers_( displacementShidaNumbers ), bodyReferenceRadius_( bodyReferenceRadius ),
        meanVectorsToDeformingBodies_( meanVectorsToDeformingBodies )
    { }

    std::vector< std::string > getDeformingBodies( ){ return deformingBodies_;}
    std::vector< int > getMaximumDeformationOrders( ){ return maximumDeformationOrders_; }
    std::vector< double > getDisplacementLoveNumbers( ){ return displacementLoveNumbers_; }
    std::vector< double > getDisplacementShidaNumbers( ){ return displacementShidaNumbers_; }
    double getBodyReferenceRadius( ){ return bodyReferenceRadius_; }
    std::vector< Eigen::Vector3d > getMeanVectorsToDeformingBodies( )
    {
        return meanVectorsToDeformingBodies_;
    }


protected:
    std::vector< std::string > deformingBodies_;
    std::vector< int > maximumDeformationOrders_;
    std::vector< double > displacementLoveNumbers_;
    std::vector< double > displacementShidaNumbers_;
    double bodyReferenceRadius_;
    std::vector< Eigen::Vector3d > meanVectorsToDeformingBodies_;
};


boost::shared_ptr< site_displacements::BodyDeformationModel > createBodyDeformationModel(
        boost::shared_ptr< BodyDeformationSettings > bodyDeformationSettings,
        const std::string body,
        const NamedBodyMap bodyMap );
}

}

#endif // CREATEBODYDEFORMATIONMODEL_H
