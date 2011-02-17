/*! \file convertMeanAnomalyBase.h
 *    This header file contains a base class to convert mean anomaly to
 *    eccentric and hyperbolic eccentric anomalies.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 1
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : E. Iorfida
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : elisabetta_iorfida@yahoo.it
 *
 *    Date created      : 10 February, 2011
 *    Last modified     : 10 February, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      110210    K. Kumar          First creation of code.
 */

#ifndef CONVERTMEANANOMALYBASE_H
#define CONVERTMEANANOMALYBASE_H

// Include statements.
#include <ctime>
#include "basicMathematicsFunctions.h"
#include "body.h"
#include "celestialBody.h"
#include "keplerianElements.h"
#include "newtonRaphson.h"
#include "newtonRaphsonAdaptor.h"
#include "orbitalElementConversions.h"

//! Orbital element conversions namespace.
/*!
 *  Orbital element conversions namespace.
 */
namespace orbital_element_conversions
{

class ConvertMeanAnomalyBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ConvertMeanAnomalyBase( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~ConvertMeanAnomalyBase( );

    //! Set Newton-Raphson method.
    /*!
     * Sets the Newton-Raphson method used.
     * \param pointerToNewtonRaphson Pointer to NewtonRaphson object.
     */
    void setNewtonRaphson( NewtonRaphson* pointerToNewtonRaphson );

    //! Set eccentricity.
    /*!
     * Sets eccentricity of orbit.
     * \param eccentricity Eccentricity.
     */
    void setEccentricity( const double& eccentricity );

    //! Set mean anomaly.
    /*!
     * Sets the mean anomaly.
     * \param meanAnomaly Mean anomaly.
     */
    void setMeanAnomaly( const double& meanAnomaly );

protected:

    //! Eccentricity.
    /*!
     * Eccentricity.
     */
    double eccentricity_;

    //! Mean anomaly.
    /*!
     * Mean anomaly.
     */
    double meanAnomaly_;

    //! Pointer to Newton-Raphson.
    /*!
     * Pointer to Newton-Raphson method.
     */
    NewtonRaphson* pointerToNewtonRaphson_;

private:
};

}

#endif // CONVERTMEANANOMALYBASE_H

// End of file.
