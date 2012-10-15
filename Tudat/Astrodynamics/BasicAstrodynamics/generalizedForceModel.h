/*    Copyright (c) 2010-2012, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      120204    D. Dirkx        File created.
 *
 *    References
 *
 */

#ifndef TUDAT_GENERALIZED_FORCE_MODEL_H
#define TUDAT_GENERALIZED_FORCE_MODEL_H

namespace tudat
{

//! Base class for all generalized forces.
/*!
 * Base class for all generalized forces (i.e. forces, moments, etc.). Base class only has
 * the interface for retrieving the generlaized force, which can be wrapped by a more
 * specific function (for instance getMoment) in derived classed. Functions to calculate
 * generalized forces have different interfaces for different derived classes and therefore
 * has no base class interface.
 */
template< typename GeneralizedForceType, int sizeOfGeneralizedForce >
class GeneralizedForceModel
{
public:

    //! Virtual destructor.
    /*!
     *  Virtual destructor.
     */
    virtual ~GeneralizedForceModel( ) { }

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

#endif // TUDAT_GENERALIZED_FORCE_MODEL_H
