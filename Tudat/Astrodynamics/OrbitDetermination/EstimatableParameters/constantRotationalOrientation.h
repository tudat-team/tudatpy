#ifndef CONSTANTROTATIONALORIENTATION_H
#define CONSTANTROTATIONALORIENTATION_H


#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"

using std::string;

namespace tudat
{

namespace estimatable_parameters
{

// pole right ascension and declination
class ConstantRotationalOrientation: public EstimatableParameter< Eigen::VectorXd >
{

public:
    ConstantRotationalOrientation(
            boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( rotation_pole_position, associatedBody ),
        rotationModel_( rotationModel ) { }

    ~ConstantRotationalOrientation( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        return rotationModel_->getInitialEulerAngles( ).segment( 0, 2 ); }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        rotationModel_->resetInitialPoleRightAscensionAndDeclination( parameterValue.x( ),
                                                                      parameterValue.y( ) );
    }

    int getParameterSize( ){ return 2; }

protected:

private:
    boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel_;
};

}

}

#endif // CONSTANTROTATIONALORIENTATION_H
