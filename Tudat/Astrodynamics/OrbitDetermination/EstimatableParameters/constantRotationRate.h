#ifndef CONSTANTROTATIONRATE_H
#define CONSTANTROTATIONRATE_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Class for the body rotation rate parameter.
/*!
 *  Class for the body rotation rate parameter. Specific parameter is the z-component of the rotation rate vector in SimpleRotationalEphemeris,
 *  assuming no x- and y- component.
 */
class RotationRate: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param rotationModel SimpleRotationalEphemeris object of which parameter is a property
     *  \param associatedBody Body of which parameter is a property.
     */
    RotationRate( boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel, const std::string& associatedBody ):
        EstimatableParameter< double >( constant_rotation_rate, associatedBody ),
        rotationModel_( rotationModel ){ }

    //! Destructor
    /*!
     *  Destructor
     */
    ~RotationRate( ) { }

    //! Get value of rotation rate.
    /*!
     *  Get value of rotation rate.
     */
    double getParameterValue( )
    {
        return rotationModel_->getRotationRate( );
    }

    //! Reset value of rotation rate.
    /*!
     *  Reset value of rotation rate.
     */
    void setParameterValue( double parameterValue )
    {
        rotationModel_->resetRotationRate( parameterValue );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value, 1 for this parameter
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! SimpleRotationalEphemeris object of which parameter is a property
    /*!
     *  SimpleRotationalEphemeris object of which parameter is a property
     */
    boost::shared_ptr< ephemerides::SimpleRotationalEphemeris > rotationModel_;
};

}

}

#endif // CONSTANTROTATIONRATE_H
