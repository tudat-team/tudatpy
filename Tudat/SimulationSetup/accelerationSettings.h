#ifndef TUDAT_ACCELERATIONSETTINGS_H
#define TUDAT_ACCELERATIONSETTINGS_H

/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      150501    D. Dirkx          Ported from personal code
 *
 *    References
 *
 *    Notes
 *
 */

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/cannonBallRadiationPressureAcceleration.h"
#include "Tudat/Astrodynamics/Gravitation/centralGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/sphericalHarmonicsGravityModel.h"
#include "Tudat/Astrodynamics/Gravitation/thirdBodyPerturbation.h"
#include "Tudat/Astrodynamics/Aerodynamics/aerodynamicAcceleration.h"


namespace tudat
{

namespace simulation_setup
{

//! Class for providing settings for acceleration model.
/*!
 *  Class for providing settings for acceleration model. This class is a functional (base) class for
 *  settings of acceleration models that  require no information in addition to their type.
 *  Classes defining settings for acceleration models requiring additional information must be
 *  derived from this class.
 *  Bodies exerting and undergong acceleration are set externally from this class.
 *  This class can be used for the easy setup of acceleration models
 *  (see createAccelerationModels.h), but users may also chose to do so manually.
 *  (Derived) Class members are all public, for ease of access and modification.
 */
class AccelerationSettings
{
public:
    //! Constructor, sets type of acceleration.
    /*!
     *  Constructor, sets type of acceleration.
     *  \param accelerationType Type of acceleration from AvailableAcceleration enum.
     */
    AccelerationSettings( const basic_astrodynamics::AvailableAcceleration accelerationType ):
        accelerationType_( accelerationType ){ }

    //! Destructor.
    virtual ~AccelerationSettings( ){ }

    //! Type of acceleration from AvailableAcceleration enum.
    basic_astrodynamics::AvailableAcceleration accelerationType_;

};

//! Class for providing settings for spherical harmonics acceleration model.
/*!
 *  Class for providing settings for spherical harmonics acceleration model,
 *  specifically the maximum degree and order up to which the field is to be expanded. Note that
 *  the minimum degree and order are currently always set to zero.
 */
class SphericalHarmonicAccelerationSettings: public AccelerationSettings
{
public:
    //! Constructor to set maximum degree and order that is to be taken into account.
    /*!
     *  Constructor to set maximum degree and order that is to be taken into account.
     *  \param maximumDegree Maximum degree
     *  \param maximumOrder Maximum order
     */
    SphericalHarmonicAccelerationSettings( const int maximumDegree,
                                           const int maximumOrder ):
        AccelerationSettings( basic_astrodynamics::spherical_harmonic_gravity ), maximumDegree_( maximumDegree ),
        maximumOrder_( maximumOrder ){ }

    //! Maximum degree that is to be used for spherical harmonic acceleration
    int maximumDegree_;

    //! Maximum order that is to be used for spherical harmonic acceleration
    int maximumOrder_;
};

//! Typedef defining a list of acceleration settings, set up in the same manner as the
//! AccelerationMap typedef.
typedef std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr<
AccelerationSettings > > > > SelectedAccelerationMap;

}

}

#endif // TUDAT_ACCELERATIONSETTINGS_H
