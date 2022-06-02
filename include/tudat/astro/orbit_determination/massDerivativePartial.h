/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_MASSDERIVATIVEPARTIAL_H
#define TUDAT_MASSDERIVATIVEPARTIAL_H

#include "tudat/astro/orbit_determination/stateDerivativePartial.h"


namespace tudat
{

namespace orbit_determination
{

class MassDerivativePartial: public StateDerivativePartial
{
public:
    MassDerivativePartial( const std::string& body,
                         const basic_astrodynamics::AvailableMassRateModels massRateType ):
        StateDerivativePartial( propagators::body_mass_state, std::make_pair( body, "" ) ),
        body_( body ), massRateType_( massRateType ) { }

    //! Virtual destructor.
    virtual ~MassDerivativePartial( ) { }

    std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >
    getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        // Initialize to empty function; 0 parameter size.
        std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int >
                partialFunction = std::make_pair( std::function< void( Eigen::Block< Eigen::MatrixXd > ) >( ), 0 );

        // Check if state dependency exists
        switch( integratedStateType )
        {
        case propagators::translational_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting mass rate derivative model, cannot have reference point on body for translational dynamics" );
            }
            // Check if propagated body corresponds to accelerated, accelerating, ro relevant third body.
            else if( stateReferencePoint.first == acceleratedBody_ )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtStateOfAcceleratedBody, this, std::placeholders::_1 ), 3 );
            }
            else if( stateReferencePoint.first == acceleratingBody_ )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtStateOfAcceleratingBody, this, std::placeholders::_1 ), 3 );
            }
            else if( isAccelerationPartialWrtAdditionalBodyNonnullptr( stateReferencePoint.first ) )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtStateOfAdditionalBody,
                                                             this, std::placeholders::_1, stateReferencePoint.first ), 3 );
            }
            break;
        }
        case propagators::rotational_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting state derivative partial acceleration model, cannot have reference point on body for body mass" );
            }
            else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtNonTranslationalStateOfAdditionalBody,
                                                             this, std::placeholders::_1, stateReferencePoint, integratedStateType, true ), 1 );
            }
            break;
        }
        case propagators::body_mass_state:
        {
            // Check if reference id is consistent.
            if( stateReferencePoint.second != "" )
            {
                throw std::runtime_error( "Error when getting state derivative partial acceleration model, cannot have reference point on body for body mass" );
            }
            else if( isStateDerivativeDependentOnIntegratedAdditionalStateTypes( stateReferencePoint, integratedStateType ) )
            {
                partialFunction = std::make_pair( std::bind( &AccelerationPartial::wrtNonTranslationalStateOfAdditionalBody,
                                                             this, std::placeholders::_1, stateReferencePoint, integratedStateType, true ), 1 );
            }
            break;
        }
        case propagators::custom_state:
        {
            break;
        }
        default:
            std::string errorMessage =
                    "Error when getting state derivative partial acceleration model, dynamics type " +
                    std::to_string( integratedStateType ) + "not recognized" ;
            throw std::runtime_error( errorMessage );
            break;
        }


        return partialFunction;
    }

//    virtual bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
//            const std::pair< std::string, std::string >& stateReferencePoint,
//            const propagators::IntegratedStateType integratedStateType )
//    {
//        return false;
//    }

    virtual bool isAccelerationPartialWrtAdditionalBodyNonnullptr( const std::string& bodyName )
    {
        return 0;
    }

    std::string getBody( ) { return body_; }


    basic_astrodynamics::AvailableAcceleration getMassRateType( )
    {
        return massRateType_;
    }

protected:

    std::string body_;

    basic_astrodynamics::AvailableMassRateModels massRateType_;
};


} // namespace orbit_determination

} // namespace tudat

#endif // TUDAT_MASSDERIVATIVEPARTIAL_H
