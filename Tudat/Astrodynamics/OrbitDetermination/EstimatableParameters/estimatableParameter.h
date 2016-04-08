#ifndef ESTIMATABLEPARAMETERS_H
#define ESTIMATABLEPARAMETERS_H

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

bool isParameterDynamicalPropertyInitialState( const EstimatebleParametersEnum parameterType );

std::string getParameterType( EstimatebleParametersEnum parameter );


typedef std::pair< EstimatebleParametersEnum, std::pair< std::string, std::string > > EstimatebleParameterIdentifier;

typedef std::vector< std::pair< std::string, EstimatebleParametersEnum > > SelectedParametersToEstimate;


//! Base class for a parameter that is to be estimated.
/*!
 *  Base class for a parameter that is to be estimated. A separate derived class is to be made for each type of parameter
 *  (i.e. gravitational parameter, h2 love number, etc. ). Body initial positions are estimated separately, not using this base class.
 *  \tparam ParameterType Data type of parameter (i.e. double or VectorXd)
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
     *  \param associatedBody Body of which parameter is a property.
     */
    EstimatableParameter( const EstimatebleParametersEnum parameterName,
                          const std::string associatedBody,
                          const std::string pointOnBodyId = ""  ):
        parameterName_( std::make_pair( parameterName, std::make_pair( associatedBody, pointOnBodyId ) ) ){ }

    //! Virtual destructor.
    /*!
     *  Virtual destructor.
     */
    virtual ~EstimatableParameter( ) { }

    //! Function to retrieve the value of the parameter
    /*!
     *  Pure virtual function to retrieve the value of the parameter
     *  \return Current value of parameter.
     */
    virtual ParameterType getParameterValue( ) = 0;

    //! Function to (re)set the value of the parameter.
    /*!
     *  Pure virtual function to (re)set the value of the parameter.
     *  \param Value to which the parameter is to be set.
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

    virtual std::string getSecondaryIdentifier( )
    {
        return "";
    }

protected:

    //! Identifier of parameter.
    /*!
     *  Identifier of parameter as a pair of parameter type and body of which parameter is a property.
     */
    EstimatebleParameterIdentifier parameterName_;
};

//! Container class for all parameters that are to be estimated.
/*!
 *  Container class for all parameters that are to be estimated, both of double and vector type.
 */
template< typename InitialStateParameterType = double >
class EstimatableParameterSet
{
public:

    //! Constructor of parameter set.
    /*!
     *  Constructor of parameter set.
     *  \param doubleParameters Vector of double parameters that are estimated.
     *  \param vectorParameters Vector of vector parameters that are estimated.
     */
    EstimatableParameterSet(
            const std::vector< boost::shared_ptr< EstimatableParameter< double > > >& estimatedDoubleParameters,
            const std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& estimatedVectorParameters,
            const std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix
            < InitialStateParameterType, Eigen::Dynamic, 1 > > > >& estimateInitialStateParameters =
            ( std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix
            < InitialStateParameterType, Eigen::Dynamic, 1 > > > >( ) ),
            const std::vector< boost::shared_ptr< EstimatableParameter< double > > >& considerDoubleParameters =
            ( std::vector< boost::shared_ptr< EstimatableParameter< double > > >( ) ),
            const std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >& considerVectorParameters =
            ( std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > >( ) ) ):
        estimatedDoubleParameters_( estimatedDoubleParameters ), estimatedVectorParameters_( estimatedVectorParameters ),
        estimateInitialStateParameters_( estimateInitialStateParameters ),
        considerDoubleParameters_( considerDoubleParameters ), considerVectorParameters_( considerVectorParameters )
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

        // Iterate over all double parameters and add to parameter size.
        for( unsigned int i = 0; i < estimatedDoubleParameters_.size( ); i++ )
        {
            doubleParameters_[ estimatedParameterSetSize_ ] = estimatedDoubleParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_, 1 ) );
            estimatedParameterSetSize_++;
        }

        // Iterate over all vector parameter, add to total number of parameters and set indices in vectorParameterIndices_
        for( unsigned int i = 0; i < estimatedVectorParameters_.size( ); i++ )
        {
            vectorParameters_[ estimatedParameterSetSize_ ] = estimatedVectorParameters_[ i ];
            parameterIndices_.push_back( std::make_pair( estimatedParameterSetSize_,
                                                         estimatedVectorParameters_[ i ]->getParameterSize( ) ) );
            estimatedParameterSetSize_ += estimatedVectorParameters_[ i ]->getParameterSize( );
        }


        // Initialize total number of consider parameters to 0.
        considerParameterSetSize_ = 0;

        // Iterate over all double parameters and add to parameter size.
        for( unsigned int i = 0; i < considerDoubleParameters.size( ); i++ )
        {
            considerParameterSetSize_ += considerDoubleParameters[ i ]->getParameterSize( );
            doubleParameters_[ i + estimatedParameterSetSize_ ] = considerDoubleParameters[ i ];
        }

        // Check input consistency.
        if( considerParameterSetSize_ != static_cast< int >( considerDoubleParameters.size( ) ) )
        {
            throw std::runtime_error( "Error when making estimatable parameter set, inconsistent parameter size of double parameter vector." );
        }

        // Iterate over all vector parameter, add to total number of parameters and set indices in vectorParameterIndices_
        int vectorParameterSize;
        for( unsigned int i = 0; i < considerVectorParameters.size( ); i++ )
        {
            vectorParameterSize = considerVectorParameters[ i ]->getParameterSize( );
            vectorParameters_[ estimatedParameterSetSize_ + considerParameterSetSize_ ] = considerVectorParameters[ i ];
            considerParameterSetSize_ += vectorParameterSize;
        }

        totalParameterSetSize_ = considerParameterSetSize_ + estimatedParameterSetSize_;
    }

    //! Function to return the total number of parameter values.
    /*!
     *  Function to return the total number of parameter values.
     *  \return Size of parameter vector
     */
    int getParameterSetSize( )
    {
        return totalParameterSetSize_;
    }

    int getEstimatedParameterSetSize( )
    {
        return estimatedParameterSetSize_;
    }

    int getConsiderParameterSetSize( )
    {
        return considerParameterSetSize_;
    }

    //! Function that returns a vector containing all current parameter values
    /*!
     *  Function that returns a vector containing all current parameter values. Double and vector parameter values are concatenated
     *  in the order in which they are in the doubleParameters_ and vectorParameters_ members.
     *  \return Vector containing all parameter values
     */
    template< typename ParameterScalar >
    Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 > getFullParameterValues( )
    {
        Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >  parameterValues =
                Eigen::Matrix< ParameterScalar, Eigen::Dynamic, 1 >::Zero( totalParameterSetSize_ );

        int currentStartIndex = 0;

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

    std::vector< boost::shared_ptr< EstimatableParameter< double > > > getEstimatedDoubleParameters( )
    {
        return estimatedDoubleParameters_;
    }

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > getEstimatedVectorParameters( )
    {
        return estimatedVectorParameters_;
    }

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    getEstimatedInitialStateParameters( )
    {
        return estimateInitialStateParameters_;
    }

    std::vector< boost::shared_ptr< EstimatableParameter< double > > > getConsiderDoubleParameters( )
    {
        return considerDoubleParameters_;
    }

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > getConsiderVectorParameters( )
    {
        return considerVectorParameters_;
    }

    std::vector< std::pair< int, int > > getParametersIndices( )
    {
        return parameterIndices_;
    }

    int getInitialDynamicalStateParameterSize( )
    {
        return initialDynamicalStateParameterSize_;
    }

protected:

    //! Total number of parameter values.
    /*!
     *  Total number of parameter values.
     */
    int totalParameterSetSize_;

    std::vector< std::pair< int, int > > parameterIndices_;

    int estimatedParameterSetSize_;

    int considerParameterSetSize_;

    //! Vector of double parameters.
    /*!
     *  Vector of double parameters.
     */
    std::map< int, boost::shared_ptr< EstimatableParameter< double > > > doubleParameters_;

    //! Vector of vector parameters.
    /*!
     *  Vector of vector parameters.
     */
    std::map< int, boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > vectorParameters_;

    std::map< int, boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    initialStateParameters_;


    std::vector< boost::shared_ptr< EstimatableParameter< double > > > estimatedDoubleParameters_;

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > estimatedVectorParameters_;

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
    estimateInitialStateParameters_;



    std::vector< boost::shared_ptr< EstimatableParameter< double > > > considerDoubleParameters_;

    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::VectorXd > > > considerVectorParameters_;

    //! List of indices denoting positions of complete vector of parameter values.
    /*!
     *  List of indices denoting positions of complete vector of parameter values, as returned by getFullParameterValues or set
     *  by resetParameterValues. Entries coincide with same index of entry in vectorParameters_. First index in each pair is the
     *  start index, the second its size.
     */
    //std::vector< std::pair< int, int > > vectorParameterIndices_; // start index; size

    int initialDynamicalStateParameterSize_;
};

//! Class for providing settings for parameters to estimate
/*!
 *  Class for providing settings for parameters to estimate. This class is a functional base class for parameters that
 *  require no information in addition to their type.
 *  This class can be used for the easy setup of parameter objects (see createEstimatableParameters.h), but users may also chose to do so manually.
 *  (Derived) Class members are all public, for ease of access and modification.
 */
class EstimatableParameterSettings
{
public:

    //! Constructor.
    /*!
     *  Constructor, takes parameter type and body of which it is a property.
     *  \param associatedBody Body of which parameter is a property.
     *  \param parameterType Type of parameter.
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

template< typename InitialStateParameterType >
class InitialTranslationalStateEstimatableParameterSettings: public EstimatableParameterSettings
{
public:
    InitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, 6, 1 > initialStateValue,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, initial_body_state ), initialStateValue_( initialStateValue ),
        centralBody_( centralBody ), frameOrientation_( frameOrientation ), isStateSet_( 1 ){ }

    InitialTranslationalStateEstimatableParameterSettings(
            const std::string& associatedBody,
            const double initialTime,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameterSettings( associatedBody, initial_body_state ), initialTime_( initialTime ), centralBody_( centralBody ),
        frameOrientation_( frameOrientation ), isStateSet_( 0 ){ }

    Eigen::Matrix< InitialStateParameterType, 6, 1 > initialStateValue_;
    double initialTime_;

    bool isStateSet_;

    std::string centralBody_;
    std::string frameOrientation_;
};

template< typename InitialStateParameterType >
std::vector< std::string > getListOfBodiesToEstimate(
        const boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< std::string > bodiesToEstimate;

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

template< typename InitialStateParameterType >
std::map< propagators::IntegratedStateType, std::vector< std::pair< std::string, std::string > > >
getListOfInitialDynamicalStateParametersEstimate(
        const boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{

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


template< typename InitialStateParameterType = double >
Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getInitialStateVectorOfBodiesToEstimate(
        const boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameters )
{
    std::vector< boost::shared_ptr< EstimatableParameter<
            Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > > initialDynamicalParameters =
            estimatableParameters->getEstimatedInitialStateParameters( );


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

#endif // ESTIMATABLEPARAMETERS_H
