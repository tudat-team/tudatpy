/*! \file centralGravityField.h
 *    Header file that defines the central gravity field class in Tudat.
 *
 *    Path              : /Astrodynamics/EnvironmentModels/
 *    Version           : 2
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Checker           : J.C.P Melman
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : J.C.P.Melman@tudelft.nl
 *
 *    Date created      : 23 June, 2011
 *    Last modified     : 1 July, 2011
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
 *      110623    K. Kumar          File created.
 *      110701    K. Kumar          Added Mercury, Saturn, Neptune.
 */

#ifndef CENTRALGRAVITYFIELD_H
#define CENTRALGRAVITYFIELD_H

// Include statements.
#include "sphericalHarmonicsGravityField.h"

//! Central gravity field class.
/*!
 * Class that defines a central gravity field.
 */
class CentralGravityField : public SphericalHarmonicsGravityField
{
public:

    //! Bodies with predefined central gravity fields.
    /*!
     * Bodies with predefined central gravity fields.
     */
    enum BodiesWithPredefinedCentralGravityFields
    {
        sun,
        mercury,
        venus,
        earth,
        moon,
        mars,
        jupiter,
        saturn,
        uranus,
        neptune
    };

    //! Default constructor.
    /*!
     * Default constructor.
     */
    CentralGravityField( );

    //! Default destructor.
    /*!
     * Default destructor.
     */
    ~CentralGravityField( );

    //! Set predefined central gravity field settings.
    /*!
     * Sets predefined central gravity field settings.
     * \param bodyWithPredefinedCentralGravityField Body with predefined
     *          central gravity field.
     */
    void setPredefinedCentralGravityFieldSettings(
        BodiesWithPredefinedCentralGravityFields
        bodyWithPredefinedCentralGravityField );

protected:

private:
};

#endif // CENTRALGRAVITYFIELD_H

// End of file.
