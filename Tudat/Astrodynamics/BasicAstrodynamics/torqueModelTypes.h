#ifndef TORQUEMODELTYPES_H
#define TORQUEMODELTYPES_H

#include "Tudat/Astrodynamics/BasicAstrodynamics/torqueModel.h"

namespace tudat
{

namespace basic_astrodynamics
{


enum AvailableTorque
{
    underfined_torque = -1
};


class TorqueSettings
{
public:

    TorqueSettings( const AvailableTorque torqueType ):torqueType_( torqueType ){ }

    virtual ~TorqueSettings( ){ }

    AvailableTorque torqueType_;

};
AvailableTorque getTorqueModelType(
        boost::shared_ptr< basic_astrodynamics::TorqueModel > torqueModel );

}

}


#endif // TORQUEMODELTYPES_H
