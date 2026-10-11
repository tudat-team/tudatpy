#ifndef TUDAT_BARYCENTER_ORIGIN_ACCELERATION_PARTIAL_H
#define TUDAT_BARYCENTER_ORIGIN_ACCELERATION_PARTIAL_H

#include "tudat/astro/basic_astro/barycenterOriginAcceleration.h"
#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"

namespace tudat
{
namespace acceleration_partials
{

//! Exact chain rule for -a_B; prescribed offset accelerations have no estimated-state dependency.
class BarycenterOriginAccelerationPartial : public AccelerationPartial
{
public:
    BarycenterOriginAccelerationPartial( const std::string& acceleratedBody,
                                         const std::shared_ptr< basic_astrodynamics::BarycenterOriginAcceleration >& model,
                                         const std::vector< std::shared_ptr< AccelerationPartial > >& externalPartials ):
        AccelerationPartial( acceleratedBody, "", model, basic_astrodynamics::barycenter_origin_acceleration ),
        externalPartials_( externalPartials )
    {}

    std::pair< std::function< void( Eigen::Block< Eigen::MatrixXd > ) >, int > getDerivativeFunctionWrtStateOfIntegratedBody(
            const std::pair< std::string, std::string >& state,
            const propagators::IntegratedStateType type ) override
    {
        std::vector< std::function< void( Eigen::Block< Eigen::MatrixXd > ) > > functions;
        int columns = 0;
        for( const auto& partial : externalPartials_ )
        {
            auto entry = partial->getDerivativeFunctionWrtStateOfIntegratedBody( state, type );
            if( entry.second )
            {
                functions.push_back( entry.first );
                columns = entry.second;
            }
        }
        if( functions.empty( ) )
        {
            return { {}, 0 };
        }
        return { [ functions, columns ]( Eigen::Block< Eigen::MatrixXd > result ) {
                    Eigen::MatrixXd temporary = Eigen::MatrixXd::Zero( 3, columns );
                    for( const auto& function : functions )
                    {
                        temporary.setZero( );
                        function( temporary.block( 0, 0, 3, columns ) );
                        result -= temporary;
                    }
                },
                 columns };
    }

    bool isStateDerivativeDependentOnIntegratedAdditionalStateTypes( const std::pair< std::string, std::string >& state,
                                                                     const propagators::IntegratedStateType type ) override
    {
        for( const auto& partial : externalPartials_ )
        {
            if( partial->isStateDerivativeDependentOnIntegratedAdditionalStateTypes( state, type ) )
            {
                return true;
            }
        }
        return false;
    }

    void update( const double time ) override
    {
        for( const auto& partial : externalPartials_ )
        {
            partial->update( time );
        }
        currentTime_ = time;
    }
    void resetCurrentTimeOfMemberObjects( ) override
    {
        for( const auto& partial : externalPartials_ )
        {
            partial->resetCurrentTime( );
        }
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter ) override
    {
        return parameterPartial( parameter );
    }
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter ) override
    {
        return parameterPartial( parameter );
    }

    void wrtPositionOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > block,
                                       const bool add = true,
                                       const int row = 0,
                                       const int column = 0 ) override
    {
        addStateBlock( acceleratedBody_, block, add, row, column, 0 );
    }
    void wrtVelocityOfAcceleratedBody( Eigen::Block< Eigen::MatrixXd > block,
                                       const bool add = true,
                                       const int row = 0,
                                       const int column = 3 ) override
    {
        addStateBlock( acceleratedBody_, block, add, row, column, 3 );
    }
    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd >, const bool = true, const int = 0, const int = 0 ) override {}

private:
    void addStateBlock( const std::string& name,
                        Eigen::Block< Eigen::MatrixXd > block,
                        const bool add,
                        const int row,
                        const int column,
                        const int offset )
    {
        auto entry = getDerivativeFunctionWrtStateOfIntegratedBody( { name, "" }, propagators::translational_state );
        if( entry.second )
        {
            Eigen::MatrixXd temporary = Eigen::MatrixXd::Zero( 3, 6 );
            entry.first( temporary.block( 0, 0, 3, 6 ) );
            block.block( row, column, 3, 3 ) += ( add ? 1.0 : -1.0 ) * temporary.block( 0, offset, 3, 3 );
        }
    }
    template< typename Parameter >
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > parameterPartial( const std::shared_ptr< Parameter >& parameter )
    {
        std::pair< std::function< void( Eigen::MatrixXd& ) >, int > combined = { {}, 0 };
        for( const auto& partial : externalPartials_ )
        {
            combined = orbit_determination::createMergedParameterPartialFunction( combined,
                                                                                  partial->getParameterPartialFunction( parameter ) );
        }
        return combined;
    }
    std::vector< std::shared_ptr< AccelerationPartial > > externalPartials_;
};

}  // namespace acceleration_partials
}  // namespace tudat
#endif
