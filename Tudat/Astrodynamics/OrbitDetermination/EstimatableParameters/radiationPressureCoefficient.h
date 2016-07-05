#ifndef RADIATIONPRESSURECOEFFICIENT_H
#define RADIATIONPRESSURECOEFFICIENT_H

#include <Eigen/Core>

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a radiation pressure coefficient
class RadiationPressureCoefficient: public EstimatableParameter< double >
{

public:
    //! Constructor.
    /*!
     * Constructor
     * \param radiationPressureInterface Object containing the radiation pressure coefficient to be estimated.
     * \param associatedBody Name of body containing the radiationPressureInterface object
     */
    RadiationPressureCoefficient(
            boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface,
            std::string& associatedBody ):
        EstimatableParameter< double >( radiation_pressure_coefficient, associatedBody ),
        radiationPressureInterface_( radiationPressureInterface )
    { }

    //! Destructor.
    ~RadiationPressureCoefficient( ) { }

    //! Function to get the current value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to get the current value of the radiation pressure coefficient that is to be estimated.
     * \return Current value of the radiation pressure coefficient that is to be estimated.
     */
    double getParameterValue( )
    {
        return radiationPressureInterface_->getRadiationPressureCoefficient( );
    }

    //! Function to reset the value of the radiation pressure coefficient that is to be estimated.
    /*!
     * Function to reset the value of the radiation pressure coefficient that is to be estimated.
     * \param parameterValue New value of the radiation pressure coefficient that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        radiationPressureInterface_->resetRadiationPressureCoefficient( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( ){ return 1; }

protected:

private:

    //! Object containing the radiation pressure coefficient to be estimated.
    boost::shared_ptr< electro_magnetism::RadiationPressureInterface > radiationPressureInterface_;
};

}

}

#endif // RADIATIONPRESSURECOEFFICIENT_H
