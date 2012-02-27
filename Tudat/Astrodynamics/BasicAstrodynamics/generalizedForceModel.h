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
 *      120204    D. Dirkx        File created.
 *
 *    References
 *
 */

#ifndef GENERALIZEDFORCEMODEL_H
#define GENERALIZEDFORCEMODEL_H

namespace tudat
{

//! Base class for all generalized forces.
/*!
 *  Base class for all generalized forces (i.e. forces, moments, etc.). Base class only has
 *  the interface for retrieving the generlaized force, which can be wrapped by a more
 *  specific function (for instance getMoment) in derived classed. Functions to calculate
 *  generalized forces have different interfaces for different derived classes and therefore
 *  has no base class interface.
 */
template< typename GeneralizedForceType, int sizeOfGeneralizedForce >
class GeneralizedForceModel
{
public:

    //! Virtual destructor.
    /*!
     *  Virtual destructor.
     */
    virtual ~GeneralizedForceModel( ){ }

    //! Function to get generalized force.
    /*!
     *  Function to get generalized force. Needs to be implemented in all derived classes
     *  to return the actual type of generalized force that is modelled.
     */
    virtual GeneralizedForceType getGeneralizedForce( ) = 0;

protected:

private:

};

} // namespace tudat

#endif // GENERALIZEDFORCEMODEL_H
