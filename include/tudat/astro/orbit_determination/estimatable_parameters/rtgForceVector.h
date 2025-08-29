/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_RTGFORCEVECTOR_H
#define TUDAT_RTGFORCEVECTOR_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/basic_astro/empiricalAcceleration.h"
#include "tudat/math/interpolators/piecewiseConstantInterpolator.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's time-independent empirical accelerations
/*!
 * Interface class for estimation of a body's time-independent empirical accelerations. Interfaces the estimation with the
 * acceleration components in the EmpiricalAcceleration class
 */
class RTGForceVector : public EstimatableParameter< Eigen::VectorXd >
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param rtgAccelerationModel Class defining properties of rtg acceleration used in propagation.
     *  \param associatedBody Body for which empirical accelerations are estimated
     */
    RTGForceVector(
            const std::shared_ptr< system_models::RTGAccelerationModel > rtgAccelerationModel,
            const std::string& associatedBody):
        EstimatableParameter< Eigen::VectorXd >( rtg_force_vector, associatedBody),
        rtgAccelerationModel_( rtgAccelerationModel ),
        parameterSize_(3)
    { }

    //! Destructor
    ~RTGForceVector( ) { }

    //! Get value of rtg acceleration components
    /*!
     *  Get value of rtg acceleration components
     *  \return Value of rtg acceleration components
     */
    Eigen::VectorXd getParameterValue( )
    {
        Eigen::Vector3d parameter = rtgAccelerationModel_->getbodyFixedForceVectorAtReferenceEpoch();
        return parameter;
    }

    //! Reset value of rtg acceleration components
    /*!
     *  Reset value of rtg acceleration components
     *  \param parameterValue New value of rtg acceleration components
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        // test size of Xd
        if( parameterValue.size() != parameterSize_ )
        {
            throw std::runtime_error( "Error when getting rtg force parameter size; inconsistent sizes found." );
        }

        // Reset components in acceleration model
        rtgAccelerationModel_->resetForceVectorAtReferenceEpoch( parameterValue );

    }

    //! Function to retrieve the size of the rtg force vector parameter (always 3)

    int getParameterSize( )
    {
        return parameterSize_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription = "RTG force vector (body-fixed) at user-defined reference epoch";
        return parameterDescription;
    }

    //! Function to retrieve list of components in rtg accelerations that are to be estimated (always 0, 1, 2).a
protected:
private:
    //! Class defining properties of rtg acceleration used in propagation.
    std::shared_ptr< system_models::RTGAccelerationModel > rtgAccelerationModel_;

    //! Number of rtg acceleration components that are to be estimated.
    int parameterSize_;

    //! List of component indices in rtg accelerations that are to be estimated.
    Eigen::Vector3i accelerationIndices_;
};




//! Interface class for estimation of a body's time-independent empirical accelerations
/*!
 * Interface class for estimation of a body's time-independent empirical accelerations. Interfaces the estimation with the
 * acceleration components in the EmpiricalAcceleration class
 */

class RTGForceVectorMagnitude : public EstimatableParameter< double >
{
public:
    //! Constructor
    /*!
     *  Constructor
     *  \param rtgAccelerationModel Class defining properties of rtg acceleration used in propagation.
     *  \param associatedBody Body for which empirical accelerations are estimated
     */
    RTGForceVectorMagnitude(
            const std::shared_ptr< system_models::RTGAccelerationModel > rtgAccelerationModel,
            const std::string& associatedBody):
        EstimatableParameter< double >( rtg_force_vector_magnitude, associatedBody),
        rtgAccelerationModel_( rtgAccelerationModel ),
        parameterSize_(1)
    { }

    //! Destructor
    ~RTGForceVectorMagnitude( ) { }

    //! Get value of rtg acceleration magnitude
    /*!
     *  Get value of rtg acceleration components
     *  \return Value of rtg acceleration components
     */
    double getParameterValue( )
    {
        return rtgAccelerationModel_->getForceVectorMagnitudeAtReferenceEpoch();
    }

    //! Reset value of rtg force magnitude
    /*!
     *  Reset value of rtg force magnitude
     *  \param parameterValue New value of rtg force magnitude
     */
    void setParameterValue( double parameterValue )
    {
        // Reset value of rtg force magnitude
        rtgAccelerationModel_->resetForceMagnitudeAtReferenceEpoch( parameterValue );

    }

    //! Function to retrieve the size of the rtg force vector parameter (always 3)

    int getParameterSize( )
    {
        return parameterSize_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription = "RTG force magnitude at user-defined reference epoch";
        return parameterDescription;
    }

    //! Function to retrieve list of components in rtg accelerations that are to be estimated (always 0, 1, 2).
    Eigen::Vector3i getIndices( )
    {
        return accelerationIndices_;
    }

protected:
private:
    //! Class defining properties of rtg acceleration used in propagation.
    std::shared_ptr< system_models::RTGAccelerationModel > rtgAccelerationModel_;

    //! Number of rtg acceleration components that are to be estimated.
    int parameterSize_;

    //! List of component indices in rtg accelerations that are to be estimated.
    Eigen::Vector3i accelerationIndices_;
};


}  // namespace estimatable_parameters

}  // namespace tudat

#endif  // TUDAT_RTGFORCEVECTOR_H
