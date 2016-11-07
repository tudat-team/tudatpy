#ifndef TUDAT_POSITIONPARTIALS_H
#define TUDAT_POSITIONPARTIALS_H

#include <vector>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include <boost/function.hpp>

#include "Tudat/Mathematics/BasicMathematics/linearAlgebraTypes.h"

#include "Tudat/Astrodynamics/Ephemerides/simpleRotationalEphemeris.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/rotationMatrixPartial.h"
#include "Tudat/Astrodynamics/OrbitDetermination/ObservationPartials/observationPartial.h"
#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"

namespace tudat
{

namespace observation_partials
{

//! Function to calculate the partial of position of a point on body wrt position of that body.
/*!
 *  Function to calculate the partial of position of a point on body wrt position of that body,
 *  with both positions expressed in the same frame.
 *  \return Requested partial (3x3 identity matrix)
 */
Eigen::Matrix3d calculatePartialOfPointPositionWrtBodyPosition( );

//! Function to calculate the partial of position of a point on a body wrt its body-fixed position
/*!
 *  Function to calculate the partial of position of a point on a body (expressed in a non-corotating, non-body fixed frame)
 *  wrt its body-fixed position, e.g. the partial of a ground station's inertial
 *  position wrt its body-fixed position.
 *  \param rotationMatrixToInertialFrame Rotation matrix from body-fixed to inertial frame.
 *  \return Partial of position of a point on a body wrt its body-fixed position
 */
Eigen::Matrix3d calculatePartialOfPointPositionWrtBodyFixedPointPosition(
        const Eigen::Matrix3d& rotationMatrixToInertialFrame );


//! Base class for calculating the partial of an inertial position wrt a parameter.
/*!
 *  Base class for calculating the partial of an inertial position wrt a parameter. A derived class is implemented for
 *  each (type of) estimatable parameter. A separate instance of the class must be made for each distinct state.
 *  Note that partials wrt parameters describing a property of a rotation matrix from a local to the inertial frame is
 *  implemented in the RotationMatrixPartial class, which is
 *  then used by the PositionPartialWrtRotationMatrixParameter derived class of this class
 */
class PositionPartial
{
public:

    //! Destructor.
    virtual ~PositionPartial( ){ }

    //! Pure virtual base class function for determining partial at current time and body state.
    /*!
     *  Pure virtual base class function for determining partial at current time and body state wrt parameter of specific
     *  implemented derived. class.
     *  \param state Current inertial state of point of which partial is to be calculated by derived class
     *  \param time Current time
     *  \return Partial of point position wrt parameter (with specific parameter determined by derived class implementation).
     */
    virtual Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartial(
            const basic_mathematics::Vector6d& state, const double time ) = 0;
};

class PositionPartialWrtPosition: public PositionPartial
{
public:

    PositionPartialWrtPosition( ){ }

    ~PositionPartialWrtPosition( ){ }

    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartial(
            const basic_mathematics::Vector6d& state,
            const double time )
    {
        return calculatePartialOfPointPositionWrtBodyPosition( );
    }
};

class PositionPartialWrtRotationMatrixParameter: public PositionPartial
{
public:
    PositionPartialWrtRotationMatrixParameter(
            const boost::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject,
            const boost::function< Eigen::Vector3d( const double ) > positionFunctionInLocalFrame ):
        rotationMatrixPartialObject_( rotationMatrixPartialObject ),
        positionFunctionInLocalFrame_( positionFunctionInLocalFrame ){ }

    ~PositionPartialWrtRotationMatrixParameter( ){ }

    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartial(
            const basic_mathematics::Vector6d& state,
            const double time )
    {
        return rotationMatrixPartialObject_->calculatePartialOfRotatedVector( time, positionFunctionInLocalFrame_( time ) );
    }
private:

    boost::shared_ptr< RotationMatrixPartial > rotationMatrixPartialObject_;

    boost::function< Eigen::Vector3d( const double ) > positionFunctionInLocalFrame_;
};


class PositionPartialWrtBodyFixedPosition: public PositionPartial
{
public:

    PositionPartialWrtBodyFixedPosition( const boost::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel ):
        bodyRotationModel_( bodyRotationModel ){ }

    ~PositionPartialWrtBodyFixedPosition( ){ }

    Eigen::Matrix< double, 3, Eigen::Dynamic > calculatePartial(
            const basic_mathematics::Vector6d& state,
            const double time )
    {
        return calculatePartialOfPointPositionWrtBodyFixedPointPosition( Eigen::Matrix3d( bodyRotationModel_->getRotationToBaseFrame( time ) ) );
    }

private:

    boost::shared_ptr< ephemerides::RotationalEphemeris > bodyRotationModel_;
};

class PositionObservationScaling: public PositionPartialScaling
{
public:
    ~PositionObservationScaling( ){ }

    void update( const std::vector< basic_mathematics::Vector6d >& linkEndStates,
                 const std::vector< double >& times,
                 const observation_models::LinkEndType fixedLinkEnd )
    {

    }

    Eigen::Matrix< double, 3, 3 > getScalingFactor( const observation_models::LinkEndType linkEndType, const observation_models::LinkEndType referenceTimeLinkEnd  )
    {
        return -Eigen::Matrix3d::Identity( );
    }

private:

};

class PositionObervationPartial: public ObservationPartial< 3 >
{

public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 3, Eigen::Dynamic >, double > > PositionObservationPartialReturnType;

    PositionObervationPartial( const boost::shared_ptr< PositionObservationScaling > positionObservationScaler,
                            const std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > >& positionPartialList,
                            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier ):
        ObservationPartial< 3 >( parameterIdentifier ), positionObservationScaler_( positionObservationScaler ), positionPartialList_( positionPartialList )
    {

    }

    ~PositionObervationPartial( ) { }

    virtual PositionObservationPartialReturnType calculatePartial(
            const std::vector< basic_mathematics::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime )
    {
        PositionObservationPartialReturnType returnPartial;

        basic_mathematics::Vector6d currentState;
        double currentTime;

        for( positionPartialIterator_ = positionPartialList_.begin( ); positionPartialIterator_ != positionPartialList_.end( );
             positionPartialIterator_++ )
        {
            if( positionPartialIterator_->first == observation_models::observed_body )
            {
                currentState = states[ 0 ];
                currentTime = times[ 0 ];
            }
            else
            {
                std::cerr<<"Error when calculating position observation partial"<<std::endl;
            }

            returnPartial.push_back(
                        std::make_pair(
                            positionObservationScaler_->getScalingFactor( positionPartialIterator_->first, linkEndOfFixedTime ) *
                            ( positionPartialIterator_->second->calculatePartial(
                                  currentState, currentTime ) ), currentTime ) );
        }

        return returnPartial;
    }

protected:

    boost::shared_ptr< PositionObservationScaling > positionObservationScaler_;


    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > > positionPartialList_;

    std::map< observation_models::LinkEndType, boost::shared_ptr< PositionPartial > >::iterator positionPartialIterator_;
};

}

}


#endif // TUDAT_POSITIONPARTIALS_H
