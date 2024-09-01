#ifndef TUDAT_PANELLEDRADIATIONPRESSUREACCELERATIONPARTIAL_H
#define TUDAT_PANELLEDRADIATIONPRESSUREACCELERATIONPARTIAL_H



#include "tudat/astro/orbit_determination/acceleration_partials/accelerationPartial.h"
#include "tudat/astro/orbit_determination/acceleration_partials/radiationPressureAccelerationPartial.h"

namespace tudat
{

namespace acceleration_partials
{

//! Class to calculate the partials of the panelled radiation pressure acceleration w.r.t. parameters and states.
class PanelledRadiationPressurePartial: public AccelerationPartial
{
public:

    //! Constructor.
    /*!
     * Constructor.
     * \param radiationPressureAcceleration Panelled radiation pressure acceleration model.
     * \param radiationPressureInterface Interface object to retrieve and/or compute the properties of
     * the panelled radiation pressure model.
     * \param acceleratedBody Name of the body undergoing acceleration.
     * \param acceleratingBody Name of the body exerting acceleration.
     */
    PanelledRadiationPressurePartial(
        const std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration > radiationPressureAcceleration,
        const std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > panelledTargetModel,
            const std::string& acceleratedBody, const std::string& acceleratingBody ):
        AccelerationPartial( acceleratedBody, acceleratingBody, basic_astrodynamics::radiation_pressure ),
        radiationPressureAcceleration_( radiationPressureAcceleration ), panelledTargetModel_( panelledTargetModel ),
        numberOfBodyFixedPanels_( panelledTargetModel_->getBodyFixedPanels( ).size( ) )
    {
        for( unsigned int i = 0; i < numberOfBodyFixedPanels_; i++ )
        {
            if( panelledTargetModel_->getBodyFixedPanels( ).at( i )->getTrackedBody( ) != "" )
            {
                trackedTargetPerPanel_[ i ] = panelledTargetModel_->getBodyFixedPanels( ).at( i )->getTrackedBody( );
                if(std::find(trackedTargetList_.begin(), trackedTargetList_.end(),
                             panelledTargetModel_->getBodyFixedPanels( ).at( i )->getTrackedBody( ) ) == trackedTargetList_.end( ) )
                {
                    trackedTargetList_.push_back( panelledTargetModel_->getBodyFixedPanels( ).at( i )->getTrackedBody( ) );
                }
            }
        }

        for( auto it : panelledTargetModel_->getSegmentFixedPanels( ) )
        {
            for( unsigned int i = 0; i < it.second.size( ); i++ )
            {
                if( it.second.at( i )->getTrackedBody( ) != "" )
                {
                    trackedTargetPerSegmentPanel_[ it.first ][ i ] = it.second.at( i )->getTrackedBody( );
                    if(std::find(trackedTargetList_.begin(), trackedTargetList_.end(),
                                 it.second.at( i )->getTrackedBody( ) ) == trackedTargetList_.end( ) )
                    {
                        trackedTargetList_.push_back( it.second.at( i )->getTrackedBody( ) );
                    }
                }
            }
        }
    }

    //! Destructor.
    ~PanelledRadiationPressurePartial( ){ }

    //! Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body undergoing acceleration
     *  and adding it to existing partial block.
     *  Update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  undergoing acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
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

    //! Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration.
    /*!
     *  Function for calculating the partial of the acceleration w.r.t. the position of body exerting acceleration and
     *  adding it to the existing partial block.
     *  The update( ) function must have been called during current time step before calling this function.
     *  \param partialMatrix Block of partial derivatives of acceleration w.r.t. Cartesian position of body
     *  exerting acceleration where current partial is to be added.
     *  \param addContribution Variable denoting whether to return the partial itself (true) or the negative partial (false).
     *  \param startRow First row in partialMatrix block where the computed partial is to be added.
     *  \param startColumn First column in partialMatrix block where the computed partial is to be added.
     */
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

    void wrtSpecularReflectivity( Eigen::MatrixXd& partial, const std::string& panelTypeId );

    void wrtDiffuseReflectivity( Eigen::MatrixXd& partial, const std::string& panelTypeId );

    //! Function for updating partial w.r.t. the bodies' positions
    /*!
     *  Function for updating common blocks of partial to current state. For the panelled radiation
     *  pressure acceleration, position partial is computed and set.
     *  \param currentTime Time at which partials are to be calculated
     */
    void update( const double currentTime = TUDAT_NAN );


    //! Function for determining if the acceleration is dependent on a non-translational integrated state.
    /*!
     *  Function for determining if the acceleration is dependent on a non-translational integrated state.
     *  \param stateReferencePoint Reference point id of propagated state
     *  \param integratedStateType Type of propagated state for which dependency is to be determined.
     *  \return True if dependency exists (non-zero partial), false otherwise.
     */
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

    //! Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a double parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency, 1 otherwise).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int >
    getParameterPartialFunction( std::shared_ptr< estimatable_parameters::EstimatableParameter< double > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        int numberOfRows = 0;

        // Check if parameter dependency exists.
        if( parameter->getParameterName( ).second.first == acceleratedBody_ )
        {
            if( parameter->getParameterName( ).second.second == acceleratingBody_ )
            {
                switch( parameter->getParameterName( ).first )
                {
                case estimatable_parameters::source_direction_radiation_pressure_scaling_factor:

                    partialFunction = std::bind( &computeRadiationPressureAccelerationWrtSourceDirectionScaling,
                                                 radiationPressureAcceleration_,
                                                 std::placeholders::_1 );
                    numberOfRows = 1;

                    break;
                case estimatable_parameters::source_perpendicular_direction_radiation_pressure_scaling_factor:

                    partialFunction = std::bind( &computeRadiationPressureAccelerationWrtSourcePerpendicularDirectionScaling,
                                                 radiationPressureAcceleration_,
                                                 std::placeholders::_1 );
                    numberOfRows = 1;

                    break;

                default:
                    break;
                }
            }
            switch( parameter->getParameterName( ).first )
                {
                case estimatable_parameters::specular_reflectivity:
                    partialFunction = std::bind( &PanelledRadiationPressurePartial::wrtSpecularReflectivity,
                                                 this, std::placeholders::_1, parameter->getParameterName( ).second.second );
                    numberOfRows = 1;
                    break;

                case estimatable_parameters::diffuse_reflectivity:
                    partialFunction = std::bind( &PanelledRadiationPressurePartial::wrtDiffuseReflectivity,
                                                 this, std::placeholders::_1, parameter->getParameterName( ).second.second );
                    numberOfRows = 1;
                    break;
                default:
                    break;
                }

        }
        return std::make_pair( partialFunction, numberOfRows );
    }

    //! Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
    /*!
     *  Function for setting up and retrieving a function returning a partial w.r.t. a vector parameter.
     *  Function returns empty function and zero size indicator for parameters with no dependency for current acceleration.
     *  \param parameter Parameter w.r.t. which partial is to be taken.
     *  \return Pair of parameter partial function and number of columns in partial (0 for no dependency).
     */
    std::pair< std::function< void( Eigen::MatrixXd& ) >, int > getParameterPartialFunction(
            std::shared_ptr< estimatable_parameters::EstimatableParameter< Eigen::VectorXd > > parameter )
    {
        std::function< void( Eigen::MatrixXd& ) > partialFunction;
        return std::make_pair( partialFunction, 0 );
    }

    Eigen::Matrix< double, 1, 3 > getCurrentCosineAnglePartial( )
    {
        return currentCosineAnglePartial_;
    }

private:

    //! Pointer to the panelled radiation pressure acceleration model.
    std::shared_ptr< electromagnetism::IsotropicPointSourceRadiationPressureAcceleration > radiationPressureAcceleration_;

    //! Pointer to the panelled radiation pressure interface.
    std::shared_ptr< electromagnetism::PaneledRadiationPressureTargetModel > panelledTargetModel_;

    unsigned int numberOfBodyFixedPanels_;

    std::map< int, std::string > trackedTargetPerPanel_;

    std::map< std::string, std::map< int, std::string > > trackedTargetPerSegmentPanel_;

    std::vector< std::string > trackedTargetList_;

    //! Current partial of acceleration w.r.t. position of body undergoing acceleration (equal to minus partial w.r.t.
    //! position of body exerting acceleration).
    Eigen::Matrix3d currentPartialWrtPosition_;

    Eigen::Matrix3d currentSourceUnitVectorPartial_;

    Eigen::Matrix< double, 1, 3 > currentRadiationPressurePositionPartial_;

    Eigen::Matrix< double, 1, 3 > currentCosineAnglePartial_;


};

}

}

#endif // TUDAT_PANELLEDRADIATIONPRESSUREACCELERATIONPARTIAL_H
