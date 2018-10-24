/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MEANMOMENTOFINERTIAPARAMETER_H
#define TUDAT_MEANMOMENTOFINERTIAPARAMETER_H

#include "Tudat/Astrodynamics/Gravitation/gravityFieldModel.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a mean moment of inertia parameter
class MeanMomentOfInertiaParameter: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param getMeanMomentOfInertia Function that returns the current mean moment of inertia
     * \param setMeanMomentOfInertia Function that resets the mean moment of inertia
     * \param associatedBody Name of body containing the gravityFieldModel object
     */
    MeanMomentOfInertiaParameter(
            const std::function< double( ) > getMeanMomentOfInertia,
            const std::function< void( const double ) > setMeanMomentOfInertia,
            const std::string& associatedBody ):
        EstimatableParameter< double >( mean_moment_of_inertia, associatedBody ),
        getMeanMomentOfInertia_( getMeanMomentOfInertia ), setMeanMomentOfInertia_( setMeanMomentOfInertia ){ }

    //! Destructor
    ~MeanMomentOfInertiaParameter( ) { }

    //! Function to get the current value of the mean moment of inertia that is to be estimated.
    /*!
     * Function to get the current value of the mean moment of inertia that is to be estimated.
     * \return Current value of the mean moment of inertia that is to be estimated.
     */
    double getParameterValue( )
    {
        return getMeanMomentOfInertia_( );
    }

    //! Function to reset the value of the mean moment of inertia that is to be estimated.
    /*!
     * Function to reset the value of the mean moment of inertia that is to be estimated.
     * \param parameterValue New value of the mean moment of inertia that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        setMeanMomentOfInertia_( parameterValue );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! Function that returns the current mean moment of inertia
    std::function< double( ) > getMeanMomentOfInertia_;

    //! Function that resets the mean moment of inertia
    std::function< void( const double ) > setMeanMomentOfInertia_;

};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_MEANMOMENTOFINERTIAPARAMETER_H
