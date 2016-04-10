#ifndef RADIATIONPRESSURECOEFFICIENT_H
#define RADIATIONPRESSURECOEFFICIENT_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

namespace tudat
{

namespace estimatable_parameters
{


class RadiationPressureCoefficient: public EstimatableParameter< double >
{

public:
    RadiationPressureCoefficient(
            boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface,
            std::string& associatedBody ):
        EstimatableParameter< double >( radiation_pressure_coefficient, associatedBody ),
        radiationPressureInterface_( radiationPressureInterface )
    { }

    ~RadiationPressureCoefficient( ) { }

    double getParameterValue( )
    { return radiationPressureInterface_->getRadiationPressureCoefficient( ); }

    void setParameterValue( double parameterValue )
    {
        radiationPressureInterface_->resetRadiationPressureCoefficient( parameterValue );
    }

    int getParameterSize( ){ return 1; }

protected:

private:
    boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface_;
};

}

}

#endif // RADIATIONPRESSURECOEFFICIENT_H
