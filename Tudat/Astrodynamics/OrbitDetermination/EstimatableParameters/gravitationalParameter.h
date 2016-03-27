#ifndef GRAVITATIONALPARAMETER_H
#define GRAVITATIONALPARAMETER_H

#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

using std::string;

namespace tudat
{

namespace estimatable_parameters
{

class GravitationalParameter: public EstimatableParameter< double >
{

public:
    GravitationalParameter( boost::shared_ptr< gravitation::GravityFieldModel >
                            gravityFieldModel, std::string associatedBody ):
        EstimatableParameter< double >( gravitational_parameter, associatedBody )
    { gravityFieldModel_ = gravityFieldModel ; }

    ~GravitationalParameter( ) { }

    double getParameterValue( )
    { return gravityFieldModel_->getGravitationalParameter( ); }

    void setParameterValue( double parameterValue )
    { gravityFieldModel_->resetGravitationalParameter( parameterValue ); }

    int getParameterSize( ){ return 1; }

protected:

private:
    boost::shared_ptr< gravitation::GravityFieldModel > gravityFieldModel_;

};

}

}


#endif // GRAVITATIONALPARAMETER_H
