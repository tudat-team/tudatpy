#ifndef TUDAT_ESTIMATABLEPARAMETERS_H
#define TUDAT_ESTIMATABLEPARAMETERS_H

#include <iomanip>
#include <iostream>
#include <string>
#include <vector>

#include <boost/shared_ptr.hpp>
#include <boost/assign/list_of.hpp>

#include <Eigen/Geometry>


#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"
namespace tudat
{

namespace estimatable_parameters
{

//! List of parameters that can be estimated by the orbit determination code.
enum EstimatebleParametersEnum
{
    initial_body_state,
    gravitational_parameter,
    constant_drag_coefficient,
    radiation_pressure_coefficient

};

//! Function to determine whether the given parameter represents an initial dynamical state, or a static parameter.
/*!
 * Function to determine whether the given parameter represents an initial dynamical state, or a static parameter.
 * \param parameterType Parameter identifier.
 * \return True if parameter is an initial dynamical state.
 */
bool isParameterDynamicalPropertyInitialState( const EstimatebleParametersEnum parameterType );

//! Function to determine whether the given (non-dynamical) parameter is a double or vector parameter.
/*!
 * Function to determine whether the given (non-dynamical) parameter is a double or vector parameter.
 * \param parameterType Parameter identifier.
 * \return True if parameter is a double parameter.
 */
bool isDoubleParameter( const EstimatebleParametersEnum parameterType );

//! Typedef for full parameter identifier.
typedef std::pair< EstimatebleParametersEnum, std::pair< std::string, std::string > > EstimatebleParameterIdentifier;



//! Base class for a parameter that is to be estimated.
/*!
 *  Base class for a parameter that is to be estimated. A separate derived class is to be made for each type of parameter
 *  (i.e. gravitational parameter, initial translational state, etc. ).
 */
template< typename ParameterType >
class EstimatableParameter
{

public:
    //! Constructor.
    /*!
     *  Constructor taking parameter name and associated body. All parameters are identified by a these two variables.
     *  Any additional information that may be required for uniquely defining a parameter is to be defined in the derived class.
     *  \param parameterName Enum value defining the type of the parameter.
     *  \param associatedBody Name of body associated with patameters
     *  \param pointOnBodyId Reference point on body associated with parameter (empty by default).
     */
    EstimatableParameter( const EstimatebleParametersEnum parameterName,
                          const std::string& associatedBody,
                          const std::string& pointOnBodyId = ""  ):
        parameterName_( std::make_pair( parameterName, std::make_pair( associatedBody, pointOnBodyId ) ) ){ }

    //! Virtual destructor.
    virtual ~EstimatableParameter( ) { }

    //! Pure virtual function to retrieve the value of the parameter
    /*!
     *  Pure virtual function to retrieve the value of the parameter
     *  \return Current value of parameter.
     */
    virtual ParameterType getParameterValue( ) = 0;

    //! Pure virtual function to (re)set the value of the parameter.
    /*!
     *  Pure virtual function to (re)set the value of the parameter.
     *  \param parameterValue to which the parameter is to be set.
     */
    virtual void setParameterValue( ParameterType parameterValue ) = 0;

    //! Function to retrieve the type and associated body of the parameter.
    /*!
     *  Function to retrieve the type and associated body of the parameter.
     *  \return Identifier of parameter as a pair of parameter type and body of which parameter is a property.
     */
    EstimatebleParameterIdentifier getParameterName( ) { return parameterName_; }

    //! Function to retrieve the size of the parameter
    /*!
     *  Pure virtual function to retrieve the size of the parameter (i.e. 1 for double parameters)
     *  \return Size of parameter value.
     */
    virtual int getParameterSize( ) = 0;

    //! Function to return additional identifier for parameter
    /*!
     *  Function to return additional identifier for parameter, beyond information stored in parameterName_, default
     *  none.
     *  \return Additional identifier for parameter (default empty string).
     */
    virtual std::string getSecondaryIdentifier( )
    {
        return "";
    }

protected:

    //! Identifier of parameter.
    EstimatebleParameterIdentifier parameterName_;
};

//! Container class for all parameters that are to be estimated.
/*!
 *  Container class for all parameters that are to be estimated. Class is templated with the scalar type used for the
 *  estimation of any initial dynamical states that may be included
 */
template< typename InitialStateParameterType = double >
class EstimatableParameterSet
{
public:

    //! Constructor of parameter set.
    /*!
     *  Constructor of parameter set.
     *  \param estimatedDoubleParameters List of double parameters that are estimated.
     *  \param estimatedVectorParameters List of vector parameters that are estimated.
     *  \param estimateInitialStateParameters List of initial dynamical states that are to be estimated.
     */
    EstimatableParameterSet(
            const std::vector< boost::shared_ptr< EstimatableParameter< double > > >& estimatedDoubleParameters,
            const std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& estimatedVectorParameters,
            const std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix
            < InitialStateParameterType, Eigen::Dynamic, 1 > > > >& estimateInitialStateParameters =
            ( std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix
            < InitialStateParameterType, Eigen::Dynamic, 1 > > > >( ) ) ):
        estimatedDoubleParameters_( estimatedDoubleParameters ), estimatedVectorParameters_( estimatedVectorParameters ),
        estimateInitialStateParameters_( estimateInitialStateParameters )
    {
        // Initialize total number of parameters to 0.
        estimatedParameterSetSize_ = 0;
        initialDynamicalStateParameterSize_ = 0;

        // Iterate over all double parameters and add to parameter size.
        for( unsigned int i = 0; i < estimateInitialStateParameters_.size( ); i++ )
        {
            initialStateParameters_[ estimatedParameterSetSize_ ] = estimateInitialStateParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimateInitialStateParameters_[ i ]->getParameterSize( ) ) );
            estimatedParameterSetSize_ += estimateInitialStateParameters_[ i ]->getParameterSize( );
            initialDynamicalStateParameterSize_ += estimateInitialStateParameters_[ i ]->getParameterSize( );
        }

        // Iterate over all double parameters and add to parameter size and set indices in parameterIndices_
        for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
        {
            doubleParameters_[ estimatedParameterSetSize_ ] = estimatedDoubleParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_, 1 ) );
            estimatedParameterSetSize_++;
        }

        // Iterate over all vector parameter, add to total number of parameters and set indices in parameterIndices_
        for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
        {
            vectorParameters_[ estimatedParameterSetSize_ ] = estimatedVectorParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimatedVectorParameters_[ i ]->getParameterSize( ) ) );
            estimatedParameterSetSize_ += estimatedVectorParameters_[ i ]->getParameterSize( );
        }

        totalParameterSetSize_ = estimatedParameterSetSize_;
    }

    //! Function to return the total number of parameter values (including consider parameters)
    /*!
     *  Function to return the total number of parameter values (including consider parameters)
     *  \return Size of parameter vector (including consider parameters)
     */
    int getParameterSetSize( )
    {
        return totalParameterSetSize_;
    }

    //! Function to return the total number of parameter values (excluding consider parameters).
    /*!
     *  Function to return the total number of parameter values (excluding consider parameters)
     *  \return Size of parameter vector (excluding consider parameters)
     */
    int getEstimatedParameterSetSize( )
    {
        return estimatedParameterSetSize_;
    }

    //! Function to return the total number of initial state values that are estimated.
    /*!
     *  Function to return the total number of initial state values that are estimated.
     *  \return Function to return the total number of initial state values that are estimated.
     */
    int getInitialDynamicalStateParameterSize( )
    {
        return initialDynamicalStateParameterSize_;
    }

    //! Function that returns a vector containing all current parameter values
    /*!
     *  Function that returns a vector containing all current parameter values. The total vector starts with the initial
     *  state parameters, followed by the double and vector parameters, respectively.
     *  Initial state, double and vector parameter values are concatenated in the order in which they are set in the
     *  estimateInitialStateParameters_, doubleParameters_ and vectorParameters_ members.
     *  \return Vector containing all parameter values
     */
    template< typename ParameterScalar >
    Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 > getFullParameterValues( )
    {
        Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >  parameterValues =
                Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >::Zero( totalParameterSetSize_ );

        int currentStartIndex = 0;

        // Retrieve initial state parameter values.
        for( unsigned int i = 0; i < estimateInitialStateParameters_.size( ); i++ )
        {
            parameterValues.segment( currentStartIndex, estimateInitialStateParameters_[ i ]->getParameterSize( ) ) =
                    estimateInitialStateParameters_[ i ]->getParameterValue( ).template cast< ParameterScalar >( );
            currentStartIndex += estimateInitialStateParameters_[ i ]->getParameterSize( );
        }

        // Retrieve double parameter values.
        for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
        {
            parameterValues( currentStartIndex ) = static_cast< ParameterScalar >(
                        estimatedDoubleParameters_[ i ]->getParameterValue( ) );
            currentStartIndex++;
        }

        // Retrieve vector parameter values.
        for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
        {
            parameterValues.segment( currentStartIndex, estimatedVectorParameters_[ i ]->getParameterSize( ) ) =
                    estimatedVectorParameters_[ i ]->getParameterValue( ).template cast< ParameterScalar >( );
            currentStartIndex += estimatedVectorParameters_[ i ]->getParameterSize( );
        }

        return parameterValues;
    }

    //! Function to reset all parameter values.
    /*!
     *  Function to reset all parameter values.
     *  \param newParameterValues New parameter values. Order of values in vector must be same order as return vector of getFullParameterValues
     */
    template< typename ParameterScalar >
    void resetParameterValues( const Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >& newParameterValues )
    {
        // Check input consistency
        if( newParameterValues.rows( ) != totalParameterSetSize_ )
        {
            throw std::runtime_error( "Error when resetting parameters of parameter set, given vector has size " +
                                      boost::lexical_cast< std::string >( newParameterValues.rows( ) ) +
                                       ", while internal size is " + boost::lexical_cast< std::string >( totalParameterSetSize_ ) );
        }
        else
        {
            int currentStartIndex = 0;

            for( unsigned int i = 0; i < estimateInitialStateParameters_.size( ); i++ )
            {
                estimateInitialStateParameters_[ i ]->setParameterValue(
                            newParameterValues.segment( currentStartIndex, estimateInitialStateParameters_[ i ]->getParameterSize( ) ).
                            template cast< InitialStateParameterType >( ) );
                currentStartIndex += estimateInitialStateParameters_[ i ]->getParameterSize( );
            }

            // Set double parameter values.
            for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
            {
                estimatedDoubleParameters_[ i ]->setParameterValue( static_cast< double >( newParameterValues( currentStartIndex ) ) );
                currentStartIndex++;
            }

            // Set vector parameter values.

            for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
            {
                estimatedVectorParameters_[ i ]->setParameterValue(
                            newParameterValues.segment( currentStartIndex, estimatedVectorParameters_[ i ]->getParameterSize( ) ).
                            template cast< double >( ) );

                currentStartIndex += estimatedVectorParameters_[ i ]->getParameterSize( );
            }
        }
    }

    //! Function to retrieve double parameter objects.
    /*!
     *  Function to retrieve double parameter objects.
     *  \return Vector containing all double parameter objects
     */
    std::map< int, boost::shared_ptr< EstimatableParameter< double > > > getDoubleParameters( )
    {
        return doubleParameters_;
    }

    //! Function to retrieve vector parameter objects.
    /*!
     *  Function to retrieve vector parameter objects.
     *  \return Vector containing all vector parameter objects
     */
    std::map< int, boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > getVectorParameters( )
    {
        return vectorParameters_;
    }

    //! Function to get list of initial dynamical states that are to be estimated.
    //!
    /*!
     *  Function to get list of initial dynamical states that are to be estimated.
     *  \return List of initial dynamical states that are to be estimated.
     */
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    getEstimatedInitialStateParameters( )
    {
        return estimateInitialStateParameters_;
    }

    //! Function to retrieve list of start indices and sizes (map keys) of estimated parameters.
    /*!
     *  Function to retrieve list of start indices and sizes (map keys) of estimated parameters.
     *  \return List of start indices and sizes (map keys) of estimated parameters.
     */
    std::vector< std::pair< int, int > > getParametersIndices( )
    {
        return parameterIndices_;
    }


protected:

    //! Total size of all initial dynamical states that is to be estimated.
    int initialDynamicalStateParameterSize_;

    //! Total number of parameter values (including currently non yet implemented consider parameters).
    int totalParameterSetSize_;

    //! Total number of estimated parameter values (excluding currently non yet implemented consider parameters).
    int estimatedParameterSetSize_;

    //! List of start indices and sizes (map keys) of estimated parameters.
    /*!
     * List of start indices and sizes (map keys) of estimated parameters, in order of vector
     * estimateInitialStateParameters_, followed by estimatedDoubleParameters_, followed by estimatedVectorParameters_.
     */
    std::vector< std::pair< int, int > > parameterIndices_;

    //! List of double parameters that are to be estimated.
    std::vector< boost::shared_ptr< EstimatableParameter< double > > > estimatedDoubleParameters_;

    //! List of vector parameters that are to be estimated.
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedVectorParameters_;

    //! List of initial dynamical states that are to be estimated.
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    estimateInitialStateParameters_;

    //! Map of double parameters that are to be estimated, with start index in total parameter vector as key.
    std::map< int, boost::shared_ptr< EstimatableParameter< double > > > doubleParameters_;

    //! Map of vector parameters that are to be estimated, with start index in total parameter vector as key.
    std::map< int, boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters_;

    //! Map of initial dynamical states that are to be estimated, with start index in total parameter vector as key.
    std::map< int, boost::shared_ptr<
    EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialStateParameters_;

};

//! Class for providing settings for parameters to estimate
/*!
 *  Class for providing settings for parameters to estimate. This class is a functional base class for parameters that
 *  require no information in addition to their type.
 *  This class can be used for the easy setup of parameter objects (see createEstimatableParameters.h), but users may also
 *  chose to do so manually. (Derived) Class members are all public, for ease of access and modification.
 */
class EstimatableParameterSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor, takes parameter type and body of which it is a property.
     *  \param associatedBody Body of which parameter is a property.
     *  \param parameterType Type of parameter.
     *  \param pointOnBodyId Reference point on body associated with parameter (empty by default).
     */
    EstimatableParameterSettings( const std::string associatedBody ,
                                  const EstimatebleParametersEnum parameterType,
                                  const std::string pointOnBodyId = "" ):
        parameterType_( std::make_pair( parameterType, std::make_pair( associatedBody, pointOnBodyId ) ) ){ }

    //! Virtual destructor
    /*!
     *  Virtual destructor
     */
    virtual ~EstimatableParameterSettings( ){ }


    //! Identifier for parameter.
    /*!
     *  Identifier for parameter, contains type of parameter and body of which parameter is a property.
     */
    EstimatebleParameterIdentifier parameterType_;

};

//! Class to define settings for estimating an initial translational state.
template< typename InitialStateParameterType >
class InitialTranslationalStateEstimatableParameterSettings: public EstimatableParameterSettings
{
public:

    //! Constructor, sets initial value of translational state.
    /*!
     * Constructor, sets initial value of translational state.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialStateValue Current value of initial state (w.r.t. centralBody)
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, 6, 1 > initialStateValue,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, initial_body_state ), initialTime_( TUDAT_NAN ),
        initialStateValue_( initialStateValue ),
        centralBody_( centralBody ), frameOrientation_( frameOrientation ){ }

    //! Constructor, without initial value of translational state.
    /*!
     * Constructor, without initial value of translational state. Current initial state is retrieved from environment
     * (ephemeris objects) during creation of parameter object.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialTime Time at which initial state is defined.
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const double initialTime,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, initial_body_state ), initialTime_( initialTime ),
        centralBody_( centralBody ), frameOrientation_( frameOrientation ){ }

    //! Time at which initial state is defined (NaN for user-defined initial state value).
    double initialTime_;

    //! Current value of initial state (w.r.t. centralBody), set manually by used.
    Eigen::Matrix< InitialStateParameterType, 6, 1 > initialStateValue_;

    //! Body w.r.t. which the initial state is to be estimated.
    std::string centralBody_;

    //! Orientation of the frame in which the state is defined.
    std::string frameOrientation_;

};

//! Function to get the list of names of bodies for which initial translational dynamical state is estimated.
/*!
 *  Function to get the list of names of bodies for which initial translational dynamical state is estimated.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return List of names of bodies for which initial state is estimated.
 */
template< typename InitialStateParameterType >
std::vector< std::string > getListOfBodiesWithTranslationalStateToEstimate(
        const boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::string > bodiesToEstimate;

    // Retrieve initial dynamical parameters.
    std::vector< boost::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state )
        {
            bodiesToEstimate.push_back(  initialDynamicalParameters.at( i )->getParameterName( ).second.first );
        }
    }

    return bodiesToEstimate;
}

//! Function to get the complete list of initial dynamical states that are to be estimated, sorted by dynamics type.
/*!
 *  Function to get the complete list of initial dynamical states that are to be estimated, sorted by dynamics type.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return Map containing dynamics type (key) and vector of pairs: list of bodies (first in pair) with reference point
 *  identifier (second in pair; empty if not relevant) for which given dynamics type is estimated.
 */
template< typename InitialStateParameterType >
std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > >
getListOfInitialDynamicalStateParametersEstimate(
        const boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    // Retrieve initial dynamical parameters.
    std::vector< boost::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > > initialDynamicalStateParametersEstimate;
    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( initialDynamicalParameters.at( i )->getParameterName( ).first == initial_body_state )
        {
            initialDynamicalStateParametersEstimate[ propagators::transational_state ].push_back(
                        initialDynamicalParameters.at( i )->getParameterName( ).second );
        }       
    }

    return initialDynamicalStateParametersEstimate;
}

//! Function to get initial state vector of estimated dynamical states.
/*!
 *  Function to get initial state vector of estimated dynamical states (i.e. presently estimated state at propagation
 *  start time.
 *  \param estimatableParameters Object containing all parameters that are to be estimated.
 *  \return State vector of estimated dynamics at propagation start time.
 */
template< typename InitialStateParameterType = double >
Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getInitialStateVectorOfBodiesToEstimate(
        const boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    // Retrieve initial dynamical parameters.
    std::vector< boost::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );

    // Initialize state vector.
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialStateVector =
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >::Zero(
                estimatableParameters->getInitialDynamicalStateParameterSize( ), 1 );

    int vectorSize = 0;
    // Iterate over list of bodies of which the partials of the accelerations acting on them are required.
    for( unsigned int i = 0; i < initialDynamicalParameters.size( ); i++ )
    {
        if( isParameterDynamicalPropertyInitialState( initialDynamicalParameters.at( i )->getParameterName( ).first ) )
        {
            int currentParameterSize = initialDynamicalParameters.at( i )->getParameterSize( );
            initialStateVector.block( vectorSize, 0, currentParameterSize, 1 ) = initialDynamicalParameters.at( i )->getParameterValue( );

            vectorSize += currentParameterSize;
        }
    }

    return initialStateVector.block( 0, 0, vectorSize, 1 );
}

}

}

#endif // TUDAT_ESTIMATABLEPARAMETERS_H
