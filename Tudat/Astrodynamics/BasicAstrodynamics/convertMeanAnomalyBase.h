/*    Copyright (c) 2010-2012 Delft University of Technology.
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
 *
 *    References
 *
 */

#ifndef TUDAT_CONVERT_MEANANOMALY_BASE_H
#define TUDAT_CONVERT_MEANANOMALY_BASE_H

#include <ctime>
#include "Tudat/Astrodynamics/Bodies/body.h"
#include "Tudat/Astrodynamics/Bodies/celestialBody.h"
#include "Tudat/Astrodynamics/States/keplerianElements.h"
#include <TudatCore/Mathematics/BasicMathematics/basicMathematicsFunctions.h>
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphson.h"
#include "Tudat/Mathematics/RootFindingMethods/newtonRaphsonAdaptor.h"

namespace tudat
{
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
    void setEccentricity( double eccentricity ) { eccentricity_ = eccentricity; }

    //! Set mean anomaly.
    /*!
     * Sets the mean anomaly.
     * \param meanAnomaly Mean anomaly.
     */
    void setMeanAnomaly( double meanAnomaly ) { meanAnomaly_ = meanAnomaly; }

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

} // namespace orbital_element_conversions
} // namespace tudat

#endif // TUDAT_CONVERT_MEANANOMALY_BASE_H
