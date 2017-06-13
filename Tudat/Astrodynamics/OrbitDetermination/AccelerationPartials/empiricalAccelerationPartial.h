#ifndef TUDAT_EMPIRICALACCELERATIONPARTIAL_H
#define TUDAT_EMPIRICALACCELERATIONPARTIAL_H

#include <iostream>

#include <boost/function.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/empiricalAccelerationCoefficients.h"
#include "Tudat/Mathematics/BasicMathematics/linearAlgebra.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/empiricalAcceleration.h"

namespace tudat
{

namespace acceleration_partials
{
//! Function determine the numerical partial derivative of the true anomaly wrt the elements of the Cartesian state
/*!
 *  unction determine the numerical partial derivative of the true anomaly wrt the elements of the Cartesian state,
 *  from the Cartesian state as input.
 *  A first-order central difference with a user-defined Cartesian state perturbation vector is used.
 *  \param cartesianElements Nominal Cartesian elements at which the partials are to be computed
 *  \param gravitationalParameter Gravitational parameter of central body around which Keplerian orbit is given
 *  \param perturbations Numerical perturbations of Cartesian state that are to be used
 *  \return Partial of Cartesian state wrt true anomaly of orbit.
 */
Eigen::Matrix< double, 1, 6 > calculateNumericalPartialOfTrueAnomalyWrtState(
        const Eigen::Vector6d& cartesianElements, const double gravitationalParameter,
        const Eigen::Vector6d& perturbations );

class EmpiricalAccelerationPartial: public AccelerationPartial
{
public:
    using AccelerationPartial::getParameterPartialFunction;

    EmpiricalAccelerationPartial(
            boost::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration,
            std::string acceleratedBody,
            std::string acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::empirical_acceleration ),
        empiricalAcceleration_( empiricalAcceleration ){ perturbations<<0.1, 0.1, 0.1, 0.001, 0.001, 0.001; }


    void wrtPositionOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPositionPartial_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPositionPartial_;
        }
    }

    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPositionPartial_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPositionPartial_;
        }
    }

    void wrtVelocityOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentVelocityPartial_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentVelocityPartial_;
        }
    }

    void wrtVelocityOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentVelocityPartial_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentVelocityPartial_;
        }
    }

    std::pair< boost::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            boost::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter );

    void update( const double currentTime = TUDAT_NAN );


    void wrtArcWiseEmpiricalAccelerationCoefficient(
            boost::shared_ptr< estimatable_parameters::ArcWiseEmpiricalAccelerationCoefficientsParameter > parameter,
            Eigen::MatrixXd& partialDerivativeMatrix )
    {
        boost::shared_ptr< interpolators::LookUpScheme< double > > currentArcIndexLookUp =
                parameter->getArcTimeLookupScheme( );

        partialDerivativeMatrix = Eigen::MatrixXd::Zero( 3, parameter->getParameterSize( ) );
        int currentArc = currentArcIndexLookUp->findNearestLowerNeighbour( currentTime_ );

        int singleArcParameterSize = parameter->getSingleArcParameterSize( );
        Eigen::MatrixXd partialWrtCurrentArcAccelerations;
        wrtEmpiricalAccelerationCoefficientFromIndices(
                    singleArcParameterSize, parameter->getIndices( ), partialWrtCurrentArcAccelerations );

        partialDerivativeMatrix.block(
                    0, currentArc * singleArcParameterSize, 3, singleArcParameterSize ) =
                partialWrtCurrentArcAccelerations;

    }

    void wrtEmpiricalAccelerationCoefficient(
            boost::shared_ptr< estimatable_parameters::EmpiricalAccelerationCoefficientsParameter > parameter,
            Eigen::MatrixXd& partial )
    {
        return wrtEmpiricalAccelerationCoefficientFromIndices( parameter->getParameterSize( ), parameter->getIndices( ), partial );
    }

    void wrtEmpiricalAccelerationCoefficientFromIndices(
            const int numberOfAccelerationComponents,
            const std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >& accelerationIndices,
            Eigen::MatrixXd& partial )
    {
        Eigen::Matrix3d fullMatrix = Eigen::Matrix3d( empiricalAcceleration_->getCurrentToInertialFrame( ) );


        partial = Eigen::Matrix< double, 3, Eigen::Dynamic >::Zero( 3, numberOfAccelerationComponents );

        int currentIndex = 0;
        double multiplier = 0.0;

        for( std::map<  basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator indexIterator = accelerationIndices.begin( );
             indexIterator != accelerationIndices.end( ); indexIterator++ )
        {
            switch( indexIterator->first )
            {
            case basic_astrodynamics::constant_empirical:
                multiplier = 1.0;
                break;
            case basic_astrodynamics::sine_empirical:
                multiplier = std::sin( empiricalAcceleration_->getCurrentTrueAnomaly( ) );
                break;
            case basic_astrodynamics::cosine_empirical:
                multiplier = std::cos( empiricalAcceleration_->getCurrentTrueAnomaly( ) );
                break;
            default:
                std::cerr<<"Error when calculating partial w.r.t. empirical accelerations, could not find functional shape "<<
                           indexIterator->first<<std::endl;
            }

            for( unsigned int i = 0; i < indexIterator->second.size( ); i++ )
            {
                partial.block( 0, currentIndex, 3, 1 ) = multiplier * fullMatrix.block ( 0, indexIterator->second.at( i ), 3, 1 );
                currentIndex++;
            }

        }
    }


private:

    boost::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration_;

    Eigen::Matrix3d currentPositionPartial_;

    Eigen::Matrix3d currentVelocityPartial_;

    Eigen::Matrix< double, 1, 6 > perturbations;

};

}
}

#endif // TUDAT_EMPIRICALACCELERATIONPARTIAL_H
