#ifndef EMPIRICALACCELERATIONCOEFFICIENTS_H
#define EMPIRICALACCELERATIONCOEFFICIENTS_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Astrodynamics/BasicAstrodynamics/empiricalAcceleration.h"
#include "Tudat/Mathematics/Interpolators/piecewiseConstantInterpolator.h"

namespace tudat
{

namespace estimatable_parameters
{

class EmpiricalAccelerationCoefficientsParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:
    EmpiricalAccelerationCoefficientsParameter(
            const boost::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration, std::string associatedBody,
            const std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate ):
        EstimatableParameter< Eigen::VectorXd >( empirical_acceleration_coefficients, associatedBody )
    {
        empiricalAcceleration_ = empiricalAcceleration;

        parameterSize_ = 0;

        for( std::map< basic_astrodynamics::EmpiricalAccelerationComponents, std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > >::
             const_iterator inputComponentIterator =  componentsToEstimate.begin( ); inputComponentIterator != componentsToEstimate.end( );
             inputComponentIterator++ )
        {
            int currentComponentIndex = static_cast< int >( inputComponentIterator->first );

            for( unsigned int i = 0; i < inputComponentIterator->second.size( ); i++ )
            {
                if( accelerationIndices_.count( inputComponentIterator->second.at( i ) ) > 0 )
                {
                    std::vector< int > currentIndices = accelerationIndices_.at( inputComponentIterator->second.at( i ) );
                    if( std::find( currentIndices.begin( ), currentIndices.end( ), currentComponentIndex )  != currentIndices.end( ) )
                    {
                        std::cerr<<"Error when creating empirical acceleration parameter object, found duplicate of component "<<
                                   inputComponentIterator->second.at( i )<<" "<<inputComponentIterator->first<<std::endl;
                    }
                }
                accelerationIndices_[ inputComponentIterator->second.at( i ) ].push_back( currentComponentIndex );
                parameterSize_++;
            }
        }

    }

    ~EmpiricalAccelerationCoefficientsParameter( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd parameter = Eigen::VectorXd::Zero( parameterSize_ );
        int currentIndex = 0;

        Eigen::Matrix3d accelerationComponents = empiricalAcceleration_->getAccelerationComponents( );
        for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator componentIterator =
             accelerationIndices_.begin( ); componentIterator != accelerationIndices_.end( ); componentIterator++ )
        {
            Eigen::Vector3d accelerations = accelerationComponents.block( 0, static_cast< int >( componentIterator->first ), 3, 1 );

            std::cout<<"Emp. acc. A "<<accelerations.transpose( )<<" "<<static_cast< int >( componentIterator->first )<<std::endl;
            for( unsigned int i = 0; i < componentIterator->second.size( ); i++ )
            {
                parameter( currentIndex ) = accelerations( componentIterator->second.at( i ) );
                currentIndex++;
                std::cout<<"Emp. acc. B "<<parameter.transpose( )<<" "<<currentIndex<<" "<<i<<std::endl;
            }
        }

        if( currentIndex != parameterSize_ )
        {
            std::cerr<<"Error when getting empirical parameter size; inconsistent sizes found "<<currentIndex<<" "<<parameterSize_<<std::endl;
        }

        std::cout<<"Current emp. acc. "<<parameter.transpose( )<<std::endl;

        return parameter;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        int currentIndex = 0;

        Eigen::Matrix3d accelerationComponents = empiricalAcceleration_->getAccelerationComponents( );
        for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator componentIterator =
             accelerationIndices_.begin( ); componentIterator != accelerationIndices_.end( ); componentIterator++ )
        {
            Eigen::Vector3d accelerations = accelerationComponents.block( 0, static_cast< int >( componentIterator->first ), 3, 1 );
            for( unsigned int i = 0; i < componentIterator->second.size( ); i++ )
            {
                accelerations( componentIterator->second.at( i ) ) = parameterValue( currentIndex );
                currentIndex++;
            }
            accelerationComponents.block( 0, static_cast< int >( componentIterator->first ), 3, 1 ) = accelerations;
        }

        empiricalAcceleration_->resetAccelerationComponents( accelerationComponents );

        if( currentIndex != parameterSize_ )
        {
            std::cerr<<"Error when getting empirical parameter size; inconsistent sizes found "<<currentIndex<<" "<<parameterSize_<<std::endl;
        }
    }

    int getParameterSize( ){ return parameterSize_; }

    std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > > getIndices( ){ return accelerationIndices_; }

protected:

private:
    boost::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration_;

    int parameterSize_;

    std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > > accelerationIndices_;
};

class ArcWiseEmpiricalAccelerationCoefficientsParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:
    ArcWiseEmpiricalAccelerationCoefficientsParameter(
            const boost::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration, std::string associatedBody,
            const std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate,
            const std::vector< double > arcStartTimeList ):
        EstimatableParameter< Eigen::VectorXd >( arc_wise_empirical_acceleration_coefficients, associatedBody ), arcStartTimeList_( arcStartTimeList )
    {
        empiricalAcceleration_ = empiricalAcceleration;
        singleArcParameterSize_ = 0;

        for( std::map< basic_astrodynamics::EmpiricalAccelerationComponents, std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > >::
             const_iterator inputComponentIterator =  componentsToEstimate.begin( ); inputComponentIterator != componentsToEstimate.end( );
             inputComponentIterator++ )
        {
            int currentComponentIndex = static_cast< int >( inputComponentIterator->first );

            for( unsigned int i = 0; i < inputComponentIterator->second.size( ); i++ )
            {
                if( accelerationIndices_.count( inputComponentIterator->second.at( i ) ) > 0 )
                {
                    std::vector< int > currentIndices = accelerationIndices_.at( inputComponentIterator->second.at( i ) );
                    if( std::find( currentIndices.begin( ), currentIndices.end( ), currentComponentIndex )  != currentIndices.end( ) )
                    {
                        std::cerr<<"Error when creating empirical acceleration parameter object, found duplicate of component "<<
                                   inputComponentIterator->second.at( i )<<" "<<inputComponentIterator->first<<std::endl;
                    }
                }
                accelerationIndices_[ inputComponentIterator->second.at( i ) ].push_back( currentComponentIndex );
                singleArcParameterSize_++;
            }
        }

        arcStartTimeList_.push_back( 1.0E300 );
        Eigen::Matrix3d currentTimeInvariantEmpiricalAccelerations = empiricalAcceleration->getAccelerationComponents( );

        for( unsigned int i = 0; i < arcStartTimeList_.size( ); i++ )
        {
            empiricalAccelerationList_.push_back( currentTimeInvariantEmpiricalAccelerations );
        }

        empiricalAccelerationInterpolator_ = boost::make_shared< interpolators::PiecewiseConstantInterpolator< double, Eigen::Matrix3d > >(
                    arcStartTimeList_, empiricalAccelerationList_ );
        typedef interpolators::OneDimensionalInterpolator< double, Eigen::Matrix3d > LocalInterpolator;

        empiricalAccelerationFunction_ = boost::bind( static_cast< Eigen::Matrix3d( LocalInterpolator::* )( const double ) >
                                                      ( &LocalInterpolator::interpolate ), empiricalAccelerationInterpolator_, _1 );

    }

    ~ArcWiseEmpiricalAccelerationCoefficientsParameter( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd parameter = Eigen::VectorXd::Zero( getParameterSize( ) );
        int currentIndex = 0;

        for( unsigned int i = 0; i < arcStartTimeList_.size( ) - 1; i++ )
        {
            Eigen::Matrix3d currentArcAccelerations = empiricalAccelerationList_.at( i );
            for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator componentIterator =
                 accelerationIndices_.begin( ); componentIterator != accelerationIndices_.end( ); componentIterator++ )
            {
                int currentBlockColumn = static_cast< int >( componentIterator->first );
                for( unsigned int j = 0; j < componentIterator->second.size( ); j++ )
                {
                    parameter( currentIndex ) = currentArcAccelerations( componentIterator->second.at( j ), currentBlockColumn );
                    currentIndex++;
                }
            }
        }

        if( currentIndex != getParameterSize( ) )
        {
            std::cerr<<"Error when getting arc-wise empirical parameter size; inconsistent sizes found "<<
                       currentIndex<<" "<<getParameterSize( )<<std::endl;
        }

        return parameter;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        int currentIndex = 0;

        for( unsigned int i = 0; i < arcStartTimeList_.size( ) - 1; i++ )
        {
            Eigen::Matrix3d currentArcAccelerations = empiricalAccelerationList_.at( i );
            for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator componentIterator =
                 accelerationIndices_.begin( ); componentIterator != accelerationIndices_.end( ); componentIterator++ )
            {
                int currentBlockColumn = static_cast< int >( componentIterator->first );
                for( unsigned int j = 0; j < componentIterator->second.size( ); j++ )
                {
                    currentArcAccelerations( componentIterator->second.at( j ), currentBlockColumn ) = parameterValue( currentIndex );
                    currentIndex++;
                }
            }
            empiricalAccelerationList_[ i ] = currentArcAccelerations;
        }
        empiricalAccelerationInterpolator_->resetDependentValues( empiricalAccelerationList_ );
        empiricalAcceleration_->resetAccelerationComponentsFunction( empiricalAccelerationFunction_ );

        if( currentIndex != getParameterSize( ) )
        {
            std::cerr<<"Error when setting arc-wise empirical parameter size; inconsistent sizes found "<<currentIndex<<" "<<
                       getParameterSize( )<<std::endl;
        }
    }

    int getParameterSize( ){ return singleArcParameterSize_ * ( arcStartTimeList_.size( ) - 1 ); }

    int getSingleArcParameterSize( ){ return singleArcParameterSize_; }

    std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > > getIndices( ){ return accelerationIndices_; }

    boost::shared_ptr< interpolators::LookUpScheme< double > > getArcTimeLookupScheme( )
    {
        return empiricalAccelerationInterpolator_->getLookUpScheme( );
    }

protected:

private:
    boost::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > empiricalAcceleration_;

    boost::shared_ptr< interpolators::PiecewiseConstantInterpolator< double, Eigen::Matrix3d > > empiricalAccelerationInterpolator_;

    boost::function< Eigen::Matrix3d( const double ) > empiricalAccelerationFunction_;

    std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > > accelerationIndices_;

    std::vector< Eigen::Matrix3d > empiricalAccelerationList_;

    std::vector< double > arcStartTimeList_;

    int singleArcParameterSize_;

};

}

}

#endif // EMPIRICALACCELERATIONCOEFFICIENTS_H
