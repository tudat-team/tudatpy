/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_EIHPARTIALS_H
#define TUDAT_EIHPARTIALS_H

#include "tudat/astro/relativity/einsteinInfeldHoffmannEquations.h"

#include "tudat/astro/orbit_determination/acceleration_partials/centralGravityAccelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

Eigen::Matrix< double, 1, 3 > calculatePartialOfPointMassPotentialWrtBodyPosition(
    const Eigen::Vector3d& relativePosition,
    const double gravitationalParameter );

class EihEquationsPartials
{
public:


    EihEquationsPartials(
            const std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations );

    void update( const double currentTime );

    void addSingleScalarTermWrtPositionPartial(
        Eigen::Matrix< double, 1, 3 >& scalarTermWrtPosition, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex );

    Eigen::Matrix< double, 1, 3 > getSingleScalarCrossTermWrtPositionPartial(
        const int bodyUndergoing, const int bodyExerting, const int bodyPartial, const int termIndex );

    void addSingleScalarTermWrtVelocityPartial(
        Eigen::Matrix< double, 1, 3 >& scalarTermWrtVelocity, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex );

    void addSingleVectorTermWrtPositionPartial(
        Eigen::Matrix3d& vectorTermWrtPosition, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex );

    Eigen::Matrix< double, 3, 3 > getSingleVectorCrossTermWrtPositionPartial(
        const int bodyUndergoing, const int bodyExerting, const int bodyPartial, const int termIndex );

    void addSingleVectorTermWrtVelocityPartial(
        Eigen::Matrix3d& vectorTermWrtVelocity, const int bodyUndergoing, const int bodyExerting, const bool wrtExerting, const int termIndex );

    std::vector< std::vector< Eigen::Matrix3d > > getCurrentTotalAccelerationsWrtPosition( )
    {
        return currentTotalAccelerationsWrtPosition_;
    }

    Eigen::Matrix3d getCurrentTotalAccelerationWrtPosition( const int bodyUndergoing, const int bodyExerting )
    {
        return currentTotalAccelerationsWrtPosition_.at( bodyUndergoing ).at( bodyExerting );
    }

    std::vector< std::vector< Eigen::Matrix3d > > getCurrentTotalAccelerationsWrtVelocity( )
    {
        return currentTotalAccelerationsWrtVelocity_;
    }

    Eigen::Matrix3d getCurrentTotalAccelerationWrtVelocity( const int bodyUndergoing, const int bodyExerting )
    {
        return currentTotalAccelerationsWrtVelocity_.at( bodyUndergoing ).at( bodyExerting );
    }

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > getCurrentTotalPotentialWrtPosition( )
    {
        return currentTotalPotentialWrtPosition_;
    }

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > getCurrentLocalPotentialWrtPosition( )
    {
        return currentLocalPotentialsWrtPosition_;
    }

    std::vector< std::vector< Eigen::Matrix3d > > getCurrentSinglePointMassAccelerationWrtPosition( )
    {
        return currentSinglePointMassAccelerationsWrtExertingPosition_;
    }

    std::vector< std::vector< Eigen::Matrix3d > > getCurrentTotalPointMassAccelerationWrtPosition( )
    {
        return currentTotalPointMassAccelerationsWrtPosition_;
    }

    std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > getEihEquations( )
    {
        return eihEquations_;
    }


    void getAccelerationWrtGamma( Eigen::MatrixXd& gammaPartial, const int bodyIndex );

    void getAccelerationWrtBeta( Eigen::MatrixXd& betaPartial, const int bodyIndex );

    void getAccelerationWrtGravitationalParameter( Eigen::MatrixXd& muPartial, const int bodyIndex, const int muIndex );

protected:

    std::shared_ptr< relativity::EinsteinInfeldHoffmannEquations > eihEquations_;

    int numberOfExertingBodies_;

    int numberOfUndergoingBodies_;

    double currentTime_;

    std::vector< std::vector< Eigen::Matrix3d > > currentTotalAccelerationsWrtPosition_;

    std::vector< std::vector< Eigen::Matrix3d > > currentTotalAccelerationsWrtPositionCrossTerms_;

    std::vector< std::vector< Eigen::Matrix3d > > currentTotalAccelerationsWrtVelocity_;


    std::vector< std::vector< double > > currentTotalAccelerationsWrtGravitationalParameter_;




    std::vector< std::vector< Eigen::Matrix3d > > currentSinglePointMassAccelerationsWrtExertingPosition_;

    std::vector< std::vector< Eigen::Matrix3d > > currentTotalPointMassAccelerationsWrtPosition_;

//    std::vector< std::vector< Eigen::Vector3d > > currentSingleAccelerationsWrtGravitationalParameter_;



    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentLocalPotentialsWrtPosition_;

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentTotalPotentialWrtPosition_;


//    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > > currentLocalPotentialsWrtGravitationalParameter_;


    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > >  currentTotalScalarTermWrtExertingPosition_;

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > >  currentTotalScalarTermWrtUndergoingPosition_;

    std::vector< std::vector< Eigen::Matrix3d > >  currentTotalVectorTermWrtExertingPosition_;

    std::vector< std::vector< Eigen::Matrix3d > >  currentTotalVectorTermWrtUndergoingPosition_;


    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > >  currentTotalScalarTermWrtExertingVelocity_;

    std::vector< std::vector< Eigen::Matrix< double, 1, 3  > > >  currentTotalScalarTermWrtUndergoingVelocity_;

    std::vector< std::vector< Eigen::Matrix3d > >  currentTotalVectorTermWrtExertingVelocity_;

    std::vector< std::vector< Eigen::Matrix3d > >  currentTotalVectorTermWrtUndergoingVelocity_;



};


//! Class to calculate the partials of the central gravitational acceleration w.r.t. parameters and states.
class EihAccelerationPartial: public AccelerationPartial
{
public:

    //! Default constructor.
    /*!
     *  Default constructor.
     *  \param gravitationalAcceleration Central gravitational acceleration w.r.t. which partials are to be taken.
     *  \param acceleratedBody Body undergoing acceleration.
     *  \param acceleratingBody Body exerting acceleration.
     */
    EihAccelerationPartial(
        const std::shared_ptr< acceleration_partials::EihEquationsPartials > fulEihPartials,
        const std::string acceleratedBody ):
        AccelerationPartial( acceleratedBody, "", basic_astrodynamics::einstein_infeld_hoffmann_acceleration ),
        fullEihPartials_( fulEihPartials )
    {
        std::vector< std::string > bodyList = fullEihPartials_->getEihEquations( )->getBodiesExertingAcceleration( );
        for( unsigned int i = 0; i < bodyList.size( ); i++ )
        {
            bodyIndices_[ bodyList.at( i ) ] = i;
        }
        acceleratedBodyIndex_ = bodyIndices_.at( acceleratedBody );
    }


    void wrtPositionOfAcceleratedBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += fullEihPartials_->getCurrentTotalAccelerationWrtPosition(
                acceleratedBodyIndex_, acceleratedBodyIndex_ );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= fullEihPartials_->getCurrentTotalAccelerationWrtPosition(
                acceleratedBodyIndex_, acceleratedBodyIndex_ );
        }
    }


    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        throw std::runtime_error( "Error when calculating EIH partial w.r.t. body exerting acceleration, no such single body exists" );
    }

    virtual void wrtPositionOfAdditionalBody(
        const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        int additionalBodyIndex = bodyIndices_.at( bodyName );
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += fullEihPartials_->getCurrentTotalAccelerationWrtPosition(
                acceleratedBodyIndex_, additionalBodyIndex );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= fullEihPartials_->getCurrentTotalAccelerationWrtPosition(
                acceleratedBodyIndex_, additionalBodyIndex );
        }
    }


    void wrtVelocityOfAcceleratedBody(
        Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += fullEihPartials_->getCurrentTotalAccelerationWrtVelocity(
                acceleratedBodyIndex_, acceleratedBodyIndex_ );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= fullEihPartials_->getCurrentTotalAccelerationWrtVelocity(
                acceleratedBodyIndex_, acceleratedBodyIndex_ );
        }
    }


    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        throw std::runtime_error( "Error when calculating EIH partial w.r.t. body exerting acceleration, no such single body exists" );
    }

    virtual void wrtVelocityOfAdditionalBody(
        const std::string& bodyName, Eigen::Block< Eigen::MatrixXd > partialMatrix,
        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        int additionalBodyIndex = bodyIndices_.at( bodyName );
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += fullEihPartials_->getCurrentTotalAccelerationWrtVelocity(
                acceleratedBodyIndex_, additionalBodyIndex );
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= fullEihPartials_->getCurrentTotalAccelerationWrtVelocity(
                acceleratedBodyIndex_, additionalBodyIndex );
        }
    }



    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes(
        const std::pair< std::string, std::string >& stateReferencePoint,
        const propagators::IntegratedStateType integratedStateType )
    {
        return 0;
    }


    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > partialFunctionPair;

        // Check dependencies.
        if( parameter->getParameterName( ).first ==  estimatable_parameters::gravitational_parameter )
        {
            if( fullEihPartials_->getEihEquations( )->getAcceleratingBodyMap( ).count( parameter->getParameterName( ).second.first ) > 0 )
            {
                // If parameter is gravitational parameter, check and create dependency function .
                partialFunctionPair =
                    std::make_pair( std::bind( &EihEquationsPartials::getAccelerationWrtGravitationalParameter, fullEihPartials_,
                                                 std::placeholders::_1, acceleratedBodyIndex_,
                                                 fullEihPartials_->getEihEquations( )->getAcceleratingBodyMap( ).at( parameter->getParameterName( ).second.first ) ), 1 );
            }
        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::ppn_parameter_gamma )
        {
            partialFunctionPair =
                std::make_pair( std::bind( &EihEquationsPartials::getAccelerationWrtGamma, fullEihPartials_,
                                           std::placeholders::_1, acceleratedBodyIndex_ ), 1 );
        }
        else if( parameter->getParameterName( ).first == estimatable_parameters::ppn_parameter_beta )
        {
            partialFunctionPair =
                std::make_pair( std::bind( &EihEquationsPartials::getAccelerationWrtBeta, fullEihPartials_,
                                           std::placeholders::_1, acceleratedBodyIndex_ ), 1 );
        }
        else
        {
            partialFunctionPair = std::make_pair( std::function< void( Eigen::MatrixXd& ) >( ), 0 );
        }

        return partialFunctionPair;
    }


    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
        std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }


    void update( const double currentTime = TUDAT_NAN )
    {
        fullEihPartials_->update( currentTime );
        currentTime_ = currentTime;
    }

    const std::shared_ptr< acceleration_partials::EihEquationsPartials > getFulEihPartials( )
    {
        return fullEihPartials_;
    }


protected:

    const std::shared_ptr< acceleration_partials::EihEquationsPartials > fullEihPartials_;

    std::map< std::string, int > bodyIndices_;

    int acceleratedBodyIndex_;
};


} // namespace acceleration_partials

} // namespace tudat

#endif // TUDAT_EIHPARTIALS_H
