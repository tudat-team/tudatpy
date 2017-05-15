#ifndef RELATIVISTICTIMESTATEDERIVATIVE_H
#define RELATIVISTICTIMESTATEDERIVATIVE_H

#include <map>
#include <vector>

#include <boost/function.hpp>
#include <boost/bind.hpp>

#include <string>

#include "Tudat/Basics/basicTypedefs.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/timeConversions.h"
#include "Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h"

namespace tudat
{

namespace propagators
{

enum RelativisticTimeStateDerivativeType
{
    first_order_barycentric_to_bodycentric = 0,
    first_order_bodycentric_to_proper_time = 1
};

//! Base class for relativistic time difference calculations.
/*!
 *  Base class for relativistic time difference calculations. Derived classes implement various kinds of relativistic
 *  time conversions, such as TCG <-> TCB (1st and 2nd order). The conversion between time scale t1 and time scale t2 are assumed to be of the
 *  form t2-t1 = integral(A(t))+B(t)
 */
class RelativisticTimeStateDerivative: public propagators::SingleStateTypeDerivative< double, double >
{
public:

    //! Base class constructor, sets common member variables for derived classes.
    /*!
     *  Base class constructor, sets common member variables for derived classes.
     *  \param bodyMap Map of body objects.
     *  \param centralBody Name of central body
     *  \param externalBodies List of external body names
     */
    RelativisticTimeStateDerivative(
            const std::pair< std::string, std::string >& referencePointId,
            const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType ):
        propagators::SingleStateTypeDerivative< double, double >( propagators::relativistic_time_rate ),
        relativisticStateDerivativeType_( relativisticStateDerivativeType ),
        referencePointId_( referencePointId )
    { }

    //!  Destructor
    /*!
       *   Destructor
       */
    virtual ~RelativisticTimeStateDerivative( ){ }



    Eigen::Matrix< double, Eigen::Dynamic, 1 > convertCurrentStateToGlobalRepresentation(
            const Eigen::Matrix< double, Eigen::Dynamic, 1 >& internalSolution, const double& time )
    {
        return internalSolution;
    }

    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > convertFromOutputSolution(
            const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& outputSolution, const double& time )
    {
        return outputSolution;
    }

    Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic > convertToOutputSolution(
            const Eigen::Matrix< double, Eigen::Dynamic, Eigen::Dynamic >& internalSolution, const double& time )
    {
        return internalSolution;
    }

    virtual int getStateSize( )
    {
        return 1;
    }

    std::string getCentralBody( )
    {
        return referencePointId_.first;
    }

    std::pair< std::string, std::string > getReferencePoint( )
    {
        return referencePointId_;
    }

    propagators::RelativisticTimeStateDerivativeType getRelativisticStateDerivativeType( )
    {
        return relativisticStateDerivativeType_;
    }

protected:

    propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType_;

    std::pair< std::string, std::string > referencePointId_;
};

class PostNewtonianRelativisticTimeStateDerivative: public RelativisticTimeStateDerivative
{
public:

    //! Base class constructor, sets common member variables for derived classes.
    /*!
     *  Base class constructor, sets common member variables for derived classes.
     *  \param bodyMap Map of body objects.
     *  \param centralBody Name of central body
     *  \param externalBodies List of external body names
     */
    PostNewtonianRelativisticTimeStateDerivative(
            const std::map< std::string, boost::shared_ptr< bodies::Body > >& bodyMap,
            const std::pair< std::string, std::string >& referencePoint,
            const std::vector< std::string >& externalBodies,
            const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType,
            const boost::function< double( const double ) > timeVariableConversionFunction = &basic_astrodynamics::doDummyTimeConversion< double >,
            const double distanceScalingFactor = 1.0 );

    //!  Destructor
    /*!
       *   Destructor
       */
    virtual ~PostNewtonianRelativisticTimeStateDerivative( ){ }


    boost::function< double( const double ) > getTimeVariableConversionFunction( )
    {
        return timeVariableConversionFunction_;
    }

protected:


    //! State function of central body.
    boost::function< basic_mathematics::Vector6d( ) > centralBodyStateFunction_;

    //! List of functions returning gravitational parameters of bodies influencing time conversion.
    std::vector< boost::function< double( ) > > externalBodyGravitationalParameterFunctions_;

    //! List of functions returning Cartesian states as function of baseFrameTime (see updateStateDerivativeModel) of bodies influencing time conversion.
    std::vector< boost::function< basic_mathematics::Vector6d( ) > > externalBodyStateFunctions_;

    boost::shared_ptr< propagators::EnvironmentUpdater< double, double > > environmentUpdater_;

    boost::function< double( const double ) > timeVariableConversionFunction_;

    double distanceScalingFactor_;
};

//! Class for generating time ephemerides of 1st order (i.e. up to c^{-2} terms) using IAU 2000 resolutions on relativity
class FirstOrderBarycentricToBodyCentricTimeStateDerivative: public PostNewtonianRelativisticTimeStateDerivative
{
public:

    //! Class contructor, sets required properties of bodies.
    /*!
     *  Class contructor, sets required properties of bodies, such as state functions and gravitational parameters.
     *  \param bodyMap Map of all bodies, from which required functions will be retrieved and set as member variables in this class
     *  \param centralBody Name of central body used in conversion
     *  \param externalBodies List of bodies causing gravitational variation of time conversion, from which 1st order contributions will
     *  be calculated and combined.
     */
    FirstOrderBarycentricToBodyCentricTimeStateDerivative(
            const std::map< std::string, boost::shared_ptr< bodies::Body > >& bodyMap,
            const std::string& centralBody,
            const std::vector< std::string >& externalBodies,
            const std::map< std::string, std::pair< int, int > > sphericalHarmonicGravityExpansions = ( std::map< std::string, std::pair< int, int > >( ) ),
            const boost::function< double( const double ) > timeVariableConversionFunction = &basic_astrodynamics::doDummyTimeConversion< double >,
            const double distanceScalingFactor = 1.0,
            const propagators::RelativisticTimeStateDerivativeType relativisticStateDerivativeType =
            propagators::first_order_barycentric_to_bodycentric);

    //! Destructor.
    virtual ~FirstOrderBarycentricToBodyCentricTimeStateDerivative( ){ }

    //! Function for evaluating the integrand which is required for the conversion.
    /*!
     *  Function for evaluating the integrand which is required for the conversion, i.e. the c^{-2} integral contribution
     *  from Soffel et al. Eq. (58) is calculated.
     *  \param currentGlobalTime Current time, in frame from which conversion is done (i.e. independent variable of equations of motion of bodies)
     *  \param integratedValue Value of integral up to current time.
     *  \return Value of integrand at current time.
     */
    virtual Eigen::Matrix< double, Eigen::Dynamic, 1 > calculateSystemStateDerivative(
            const double currentGlobalTime, const Eigen::Matrix< double, Eigen::Dynamic, 1 >& integratedValue );

    //! Function to update all environment variables to current time.
    /*!
     *  Function to update environment variables (central body velocity and external potential) to current time.
     *  \param baseFrameTime Current time, in frame from which conversion is done (i.e. independent variable of equations of motion of bodies)
     */
    virtual void updateStateDerivativeModel( double baseFrameTime );

    double getCurrentExternalBodyDistance( const int index )
    {
        return currentExternalBodyDistances_.at( index );
    }

    double getCurrentExternalBodyGravitationalParameter( const int index )
    {
        return currentExternalBodyGravitationalParameters_.at( index );
    }

    basic_mathematics::Vector6d getCurrentExternalBodyStates( const int index )
    {
        return currentExternalBodyStates_.at( index );
    }

    std::vector< std::string > getExternalBodies( )
    {
        return externalBodies_;
    }

    basic_mathematics::Vector6d getCurrentCentralBodyState( )
    {
        return currentCentralBodyState_;
    }


protected:



    //! Function to update all environment variables common to this class, and its (2nd order) derived class to current time.
    /*!
     *  Function to update all environment variables common to this class, and its (2nd order) derived class to current time, namely
     *  body states and gravitational parameters.
     *  \param baseFrameTime Current time, in frame from which conversion is done (i.e. independent variable of equations of motion of bodies)
     */
    void updateBaseVariables( double baseFrameTime );

    //! Current speed of central body
    double currentVelocity_;

    //! Current state of central body
    basic_mathematics::Vector6d currentCentralBodyState_;

    //! Current distances from external bodies to central body.
    std::vector< double > currentExternalBodyDistances_;

    //! Current gravitational parameters of external bodies.
    std::vector< double > currentExternalBodyGravitationalParameters_;

    //! Current states of external bodies.
    std::vector< basic_mathematics::Vector6d > currentExternalBodyStates_;

    std::vector< std::string > externalBodies_;

    //! Current scalar potential due to external bodies.
    double currentExternalPotential_;

    std::map< int, boost::function< double( Eigen::Vector3d ) > > higherOrderGravityFieldPotentialFunctions_;// only dergee > 0 term contributions, added for 1st order variation.

    std::map< int, boost::function< double( Eigen::Vector3d ) > >::iterator shPotentialIterator_;

    std::map< int, double > fullShPotentials_;
};

//! Class for generating time ephemerides of 1st order (i.e. up to c^{-2} terms) of ground station w.r.t. a PCRS
/*!
 *  Class for generating time ephemerides of 1st order (i.e. up to c^{-2} terms) of ground station w.r.t. a PCRS. Conversion is based on
 *  Eq. (22) of Turyshev et al. (2012)
 */
class FirstOrderBodyCentricToTopoCentricTimeCalculator: public PostNewtonianRelativisticTimeStateDerivative
{
public:

    //! Constructor of bodycentric<->topocentric time calculator for a single ground station on a celestial body.
    /*!
     *  Constructor of bodycentric<->topocentric time calculator for a single ground station on a celestial body. The class is
     *  basically a wrapper for automatically calculating and setting the result of Eq. (22) of Turyshev et al. (2013). The conversion
     *  to/from a barycentric frame is done by a separate instance of the FirstOrderBarycentricToBodyCentricTimeStateDerivative or
     *  SecondOrderBarycentricToBodyCentricTimeStateDerivative class.
     *  Function to add ground stations to converter, for inclusion of topocentic<->bodycentric times. For all added ground stations,
     *  automatic conversion between any of the topocentric, bodycentric and barycentric frames can be achieved.
     *  \param bodyMap List of all bodies in simulation.
     *  \param centralBody Name of body on which ground station is located.
     *  \param groundStation Name of ground station on centralBody for which topocentric frame is to be defined.
     *  \param maximumSphericalHarmonicDegree Value denoting maximum sh gravity degree of central body to use in evaluation of
     *  gravitational potential at ground station.
     *  \param externalBodies List of bodies for which the tidal contribution to the potential is to be included in the body<->topo-
     *  centric conversion.
     *  \param useAccelerationTerm Boolean denoting whether to use acceleration term in Eq. (22) of Turyshev et al.(2013)
     *  \param useTimeDependentBodyFixedPosition Boolean denoting whether use time-dependency of ground station position, or use its nominal
     *  position at all times.
     */
    FirstOrderBodyCentricToTopoCentricTimeCalculator( const std::map< std::string, boost::shared_ptr< bodies::Body > >& bodyMap,
                                                      const std::string& centralBody,
                                                      const std::vector< std::string >& externalBodies,
                                                      const std::string& groundStation,
                                                      const int maximumSphericalHarmonicDegree = 0,
                                                      const bool useAccelerationTerm = 0,
                                                      const bool useTimeDependentBodyFixedPosition = 0 );

    //! Destructor.
    ~FirstOrderBodyCentricToTopoCentricTimeCalculator( ){ }

    //! Function for evaluating the integrand which is required for the conversion.
    /*!
     *  Function for evaluating the integrand which is required for the conversion, i.e. the integral contribution
     *  from Turyshev et al. (2013) Eq. (22) is calculated.
     *  \param currentGlobalTime Current time, in frame from which conversion is done (i.e. bodycentric coordinate time)
     *  \param integratedValue Value of integral up to current time.
     *  \return Value of integrand at current time.
     */
    Eigen::Matrix< double, Eigen::Dynamic, 1 > calculateSystemStateDerivative(
            const double currentGlobalTime, const Eigen::Matrix< double, Eigen::Dynamic, 1 >& integratedValue );


    Eigen::Vector3d getGroundStationPositionInBodyCenteredInertialFrame( const double time )
    {
        return toInertialFrameTransformation_( time ) * pointPositionFunctionInPcrs_( time );
    }


    //! Function to update all environment variables to current time.
    /*!
     *  Function to update environment variables to current time (central body potential, station state and central body state).
     *  \param baseFrameTime Current time, in frame from which conversion is done (i.e. bodycentric coordinate time).
     */
    void updateStateDerivativeModel( const double baseFrameTime );

private:


    boost::function< Eigen::Vector3d( const double ) > pointPositionFunctionInPcrs_;

    boost::function< Eigen::Vector3d( const double ) > centralBodyRotationVector_;

    boost::function< double( Eigen::Vector3d ) > localCentralBodyPotentialFunction_;

    boost::function< Eigen::Vector3d( const double ) > centralBodyAccelerationFunction_;

    boost::function< Eigen::Quaterniond( const double ) > toInertialFrameTransformation_;


    std::vector< double > currentExternalBodyDistances_;

    std::vector< double > currentExternalBodyGravitationalParameters_;

    std::vector< Eigen::Vector3d > currentExternalBodyRelativePositions_;


    Eigen::Vector3d currentCentralBodyBarycentricVelocity_;

    Eigen::Vector3d currentBarycentricAccelerationOfCentralBody_;

    Eigen::Vector3d currentPointPositionInPcrs_;

    Eigen::Vector3d currentPointVelocityInPcrs_;

    double currentLocalPotential_;
};

}

}
#endif // RELATIVISTICTIMESTATEDERIVATIVE_H
