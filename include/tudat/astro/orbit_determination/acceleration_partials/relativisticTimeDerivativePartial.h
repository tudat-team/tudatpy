/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license.
 */

#ifndef TUDAT_RELATIVISTICTIMEDERIVATIVEPARTIAL_H
#define TUDAT_RELATIVISTICTIMEDERIVATIVEPARTIAL_H

#include <functional>
#include <map>
#include <string>

#include <Eigen/Core>

#include "tudat/astro/orbit_determination/stateDerivativePartial.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/simulation/propagation_setup/propagationSettings.h"

namespace tudat
{

namespace orbit_determination
{

namespace partial_derivatives
{

//! Base class for partial derivatives of relativistic time state-derivative models.
class RelativisticTimeDerivativePartial : public StateDerivativePartial
{
public:
    //! Constructor.
    /*!
     *  \param referencePoint Reference point for proper-time state (body, optional point identifier).
     *  \param relativisticStateDerivativeType Relativistic time-derivative model type.
     */
    RelativisticTimeDerivativePartial(
            const std::pair< std::string, std::string >& referencePoint,
            const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType ):
        StateDerivativePartial( propagators::proper_time, referencePoint ),
        relativisticStateDerivativeType_( relativisticStateDerivativeType )
    { }

    ~RelativisticTimeDerivativePartial( ) override = default;

    //! Retrieve function for partial w.r.t. state of an integrated body.
    /*!
     *  \param stateReferencePoint Body/reference-point identifier of integrated state.
     *  \param integratedStateType Integrated state type.
     *  \return Pair of (partial function, number of output columns).
     */
    std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) override
    {
        std::function< void( Eigen::Block< Eigen::MatrixXd > ) > partialFunction;
        int numberOfColumns = 0;
        if( integratedStateType == propagators::translational_state )
        {
            if( !stateReferencePoint.second.empty( ) )
            {
                throw std::runtime_error(
                            "Error when requesting relativistic time derivative partials: "
                            "translational reference points may not specify a body-fixed point." );
            }
            else if( isStateDerivativePartialWrtTranslationalStateNonNull( stateReferencePoint.first ) )
            {
                const auto matrixFunction =
                        getDerivativeFunctionWrtTranslationalStateOfBody( stateReferencePoint.first );
                numberOfColumns = 6;
                partialFunction = [ matrixFunction ]( Eigen::Block< Eigen::MatrixXd > partialBlock )
                {
                    partialBlock = matrixFunction( );
                };
            }
        }

        return std::make_pair( partialFunction, numberOfColumns );
    }

    //! Check whether this partial depends on an integrated additional-state type.
    /*!
     *  \param stateReferencePoint Body/reference-point identifier of integrated state.
     *  \param integratedStateType Integrated state type.
     *  \return True if a non-zero dependency exists.
     */
    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType ) override
    {
        if( integratedStateType == propagators::translational_state )
        {
            if( !stateReferencePoint.second.empty( ) )
            {
                throw std::runtime_error(
                            "Error when checking relativistic time derivative partial dependency: "
                            "translational reference points may not specify a body-fixed point." );
            }
            return isStateDerivativePartialWrtTranslationalStateNonNull( stateReferencePoint.first );
        }

        return false;
    }

    //! Retrieve function for partial w.r.t. translational state of a body.
    /*!
     *  \param bodyName Name of body for which partial is requested.
     *  \return Function returning a 1x6 partial matrix.
     */
    virtual std::function< Eigen::Matrix< double, 1, 6 >( ) > getDerivativeFunctionWrtTranslationalStateOfBody(
            const std::string& bodyName ) = 0;

    //! Check whether translational-state partial for a body is non-zero.
    /*!
     *  \param bodyName Name of body for which partial is requested.
     *  \return True if partial is non-zero.
     */
    virtual bool isStateDerivativePartialWrtTranslationalStateNonNull( const std::string& bodyName ) = 0;

protected:
    propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType_;
};

}  // namespace partial_derivatives

}  // namespace orbit_determination

}  // namespace tudat

#endif  // TUDAT_RELATIVISTICTIMEDERIVATIVEPARTIAL_H
