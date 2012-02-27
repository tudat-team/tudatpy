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
 *      101111    E. Iorfida        First creation of code.
 *
 *    References
 *
 */

#ifndef TUDAT_TRAJECTORY_DESIGN_METHOD_H
#define TUDAT_TRAJECTORY_DESIGN_METHOD_H

namespace tudat
{

//! Trajectory design method base class.
/*!
 * Trajectory design method class.
 */
class TrajectoryDesignMethod
{
public:

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~TrajectoryDesignMethod( ) { }

    //! Execute trajectory design method.
    /*!
     * Execute trajectory design method.
     */
    virtual void execute( ) = 0;

protected:

private:
};

} // namespace tudat

#endif // TUDAT_TRAJECTORY_DESIGN_METHOD_H
