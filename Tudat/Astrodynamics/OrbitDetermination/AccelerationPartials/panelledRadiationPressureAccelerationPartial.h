#ifndef TUDAT_PANELLEDRADIATIONPRESSUREACCELERATIONPARTIAL_H
#define TUDAT_PANELLEDRADIATIONPRESSUREACCELERATIONPARTIAL_H

#include <boost/shared_ptr.hpp>

#include "Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/accelerationPartial.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/panelledRadiationPressure.h"

namespace tudat
{

namespace acceleration_partials
{

class PanelledRadiationPressurePartial: public AccelerationPartial
{
public:
    PanelledRadiationPressurePartial(
            const std::shared_ptr< electro_magnetism::PanelledRadiationPressureAcceleration > radiationPressureAcceleration,
            const std::shared_ptr< electro_magnetism::PanelledRadiationPressureInterface > radiationPressureInterface,
            const std::string& acceleratedBody, const std::string& acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::panelled_radiation_pressure_acceleration ),
        radiationPressureAcceleration_( radiationPressureAcceleration ), radiationPressureInterface_( radiationPressureInterface )
    { }

    ~PanelledRadiationPressurePartial( ){ }

    void wrtPositionOfAcceleratedBody(
            Eigen::Block< Eigen::MatrixXd > partialMatrix,
            const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
    }

    void wrtPositionOfAcceleratingBody( Eigen::Block< Eigen::MatrixXd > partialMatrix,
                                        const bool addContribution = 1, const int startRow = 0, const int startColumn = 0 )
    {
        if( addContribution )
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) -= currentPartialWrtPosition_;
        }
        else
        {
            partialMatrix.block( startRow, startColumn, 3, 3 ) += currentPartialWrtPosition_;
        }
    }

    void update( const double currentTime = TUDAT_NAN );

    bool isStateDerivativeDependentOnIntegratedNonTranslationalState(
            const std::pair< std::string, std::string >& stateReferencePoint,
            const propagators::IntegratedStateType integratedStateType )
    {
        if( ( stateReferencePoint.first == acceleratedBody_ )
                && ( integratedStateType == propagators::body_mass_state ||
                     integratedStateType == propagators::rotational_state ) )
        {
            throw std::runtime_error(
                        "Warning, dependency of panelled radiation pressure on body masses and rotational state not yet implemented" );
        }
        return 0;
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }


private:

    std::shared_ptr< electro_magnetism::PanelledRadiationPressureAcceleration > radiationPressureAcceleration_;

    std::shared_ptr< electro_magnetism::PanelledRadiationPressureInterface > radiationPressureInterface_;

    Eigen::Matrix3d currentPartialWrtPosition_;
};

}

}

#endif // TUDAT_PANELLEDRADIATIONPRESSUREACCELERATIONPARTIAL_H
