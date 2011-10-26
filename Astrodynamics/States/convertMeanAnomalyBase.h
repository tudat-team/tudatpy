/*! \file convertMeanAnomalyBase.h
 *    This header file contains a base class to convert mean anomaly to eccentric and hyperbolic
 *    eccentric anomalies.
 *
 *    Path              : /Astrodynamics/States/
 *    Version           : 2
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
 *    Last modified     : 10 August, 2011
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
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
 *      110810    J. Leloux         Corrected doxygen documentation.
 */

#ifndef CONVERTMEANANOMALYBASE_H
#define CONVERTMEANANOMALYBASE_H

// Include statements.
#include <ctime>
#include "Astrodynamics/Bodies/body.h"
#include "Astrodynamics/Bodies/celestialBody.h"
#include "Astrodynamics/States/keplerianElements.h"
#include "Astrodynamics/States/orbitalElementConversions.h"
#include "Mathematics/basicMathematicsFunctions.h"
#include "Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

//! Tudat library namespace.
/*!
 * The Tudat library namespace.
 */
namespace tudat
{

//! Orbital element conversions namespace.
/*!
 * Orbital element conversions namespace.
 */
namespace orbital_element_conversions
{

//! Definition of mean anomaly conversion base class.
/*!
 * Definition of mean anomaly conversion base class.
 */
class ConvertMeanAnomalyBase
{
public:

    //! Default constructor.
    /*!
     * Default constructor.
     */
    ConvertMeanAnomalyBase( ) : eccentricity_( -1.0 ), meanAnomaly_( -1.0 ),
        pointerToNewtonRaphson_( NULL ) { }

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~ConvertMeanAnomalyBase( ) { }

    //! Set Newton-Raphson method.
    /*!
     * Sets the Newton-Raphson method used.
     * \param pointerToNewtonRaphson Pointer to NewtonRaphson object.
     */
    void setNewtonRaphson( NewtonRaphson* pointerToNewtonRaphson )
    { pointerToNewtonRaphson_ = pointerToNewtonRaphson; }

    //! Set eccentricity.
    /*!
     * Sets eccentricity of orbit.
     * \param eccentricity Eccentricity.
     */
    void setEccentricity( const double& eccentricity ) { eccentricity_ = eccentricity; }

    //! Set mean anomaly.
    /*!
     * Sets the mean anomaly.
     * \param meanAnomaly Mean anomaly.
     */
    void setMeanAnomaly( const double& meanAnomaly ) { meanAnomaly_ = meanAnomaly; }

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

}

#endif // CONVERTMEANANOMALYBASE_H

// End of file.
