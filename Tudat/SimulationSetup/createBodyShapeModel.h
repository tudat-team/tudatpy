#ifndef CREATEBODYSHAPEMODEL_H
#define CREATEBODYSHAPEMODEL_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/BasicAstrodynamics/bodyShapeModel.h"
#include "Tudat/SimulationSetup/body.h"


namespace tudat
{

namespace simulation_setup
{

enum BodyShapeTypes
{
    spherical,
    spherical_spice,
    oblate_spheroid
};

class BodyShapeSettings
{
public:
    BodyShapeSettings( BodyShapeTypes bodyShapeType ):bodyShapeType_( bodyShapeType ){ }

    virtual ~BodyShapeSettings( ){ }

    BodyShapeTypes getBodyShapeType( ){ return bodyShapeType_; }

protected:
    BodyShapeTypes bodyShapeType_;
};

class SphericalBodyShapeSettings: public BodyShapeSettings
{
public:
    SphericalBodyShapeSettings( const double radius ):BodyShapeSettings( spherical ), radius_( radius ){ }

    double getRadius( ){ return radius_; }

private:
    double radius_;
};

class OblateSphericalBodyShapeSettings: public BodyShapeSettings
{
public:
    OblateSphericalBodyShapeSettings( const double equatorialRadius,
                                      const double flattening ):
        BodyShapeSettings( oblate_spheroid ), equatorialRadius_( equatorialRadius ), flattening_( flattening ){ }

    double getEquatorialRadius( ){ return equatorialRadius_; }

    double getFlattening( ){ return flattening_; }

private:
    double equatorialRadius_;

    double flattening_;
};

boost::shared_ptr< basic_astrodynamics::BodyShapeModel > createBodyShapeModel(
        const boost::shared_ptr< BodyShapeSettings > shapeSettings,
        const std::string& body );


}

}

#endif // CREATEBODYSHAPEMODEL_H
