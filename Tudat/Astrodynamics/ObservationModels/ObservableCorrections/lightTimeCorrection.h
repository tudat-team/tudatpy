/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_LIGHTTIMECORRECTION_H
#define TUDAT_LIGHTTIMECORRECTION_H

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

namespace tudat
{

namespace observation_models
{

//! Enum defining different types of light time corrections.
enum LightTimeCorrectionType
{
    first_order_relativistic,
    function_wrapper_light_time_correction
};


//! Base class for computing deviations from the light-time w.r.t. straight-line propagation at constant velocity (c).
/*!
 *  Base class for computing deviations from the light-time w.r.t. straight-line propagation at constant velocity (c).
 *  This base class is non-functional, and each time of light-time correction must be defined in a dedicated derived class
 */
class LightTimeCorrection
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param lightTimeCorrectionType Type of light-time correction represented by instance of class.
     */
    LightTimeCorrection( const LightTimeCorrectionType lightTimeCorrectionType ):
        lightTimeCorrectionType_( lightTimeCorrectionType ){ }

    //! Destructor
    virtual ~LightTimeCorrection( ){ }

    //! Pure virtual function to compute the light-time correction
    /*!
     * Pure virtual function to compute the light-time correction, function is to be implemented in derived class
     * for specific correction model. The input is the states and times of the two link ends (at the current iteration of
     * the solution of the implicit light-time equation), which is the information of the link already computed by the
     * light time calculator, ensuring that no double computations are performed.
     * \param transmitterState State of transmitted at transmission time
     * \param receiverState State of receiver at reception time
     * \param transmissionTime Time of signal transmission
     * \param receptionTime Time of singal reception
     * \return
     */
    virtual double calculateLightTimeCorrection(
            const basic_mathematics::Vector6d& transmitterState,
            const basic_mathematics::Vector6d& receiverState,
            const double transmissionTime,
            const double receptionTime ) = 0;

    //! Function to retrieve the type of light-time correction represented by instance of class.
    /*!
     *  Function to retrieve the type of light-time correction represented by instance of class.
     *  \return Type of light-time correction represented by instance of class.
     */
    LightTimeCorrectionType getLightTimeCorrectionType( )
    {
        return lightTimeCorrectionType_;
    }

protected:

    //! Type of light-time correction represented by instance of class.
    LightTimeCorrectionType lightTimeCorrectionType_;

};

} // namespace observation_models

} // namespace tudat

#endif // TUDAT_LIGHTTIMECORRECTION_H
