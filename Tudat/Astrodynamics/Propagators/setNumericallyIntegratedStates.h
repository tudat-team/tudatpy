#ifndef SETNUMERICALLYINTEGRATEDSTATES_H
#define SETNUMERICALLYINTEGRATEDSTATES_H

#include "Tudat/SimulationSetup/body.h"
#include "Tudat/Astrodynamics/Ephemerides/frameManager.h"
#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Propagators/propagationSettings.h"


namespace tudat
{

namespace propagators
{

template< typename TimeType, typename StateScalarType >
void resetIntegratedEphemerisOfBody(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > > ephemerisInterpolator,
        const std::string bodyToIntegrate )
{
    using namespace tudat::interpolators;
    using namespace tudat::ephemerides;

    // If body has not ephemeris, give error message.
    if( bodyMap.at( bodyToIntegrate )->getEphemeris( ) == NULL )
    {
        std::cerr<<"Error, no ephemeris detected for body "<<bodyToIntegrate<<" when resetting ephemeris"<<std::endl;
    }

    // If current ephemeris is not already a tabulated ephemeris, create new ephemeris.
    else if( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                 bodyMap.at( bodyToIntegrate )->getEphemeris( ) ) == NULL )
    {
        std::cerr<<"Error A when resetting integrated ephemeris of body "<<std::endl;

    }
    // Else, update existing tabulated ephemeris
    else
    {

        boost::shared_ptr< TabulatedCartesianEphemeris< StateScalarType, TimeType > > tabulatedEphemeris =
                boost::dynamic_pointer_cast< TabulatedCartesianEphemeris<  StateScalarType, TimeType > >(
                    bodyMap.at( bodyToIntegrate )->getEphemeris( ) );
        tabulatedEphemeris->reset( ephemerisInterpolator );
    }
}

// Set state history of current body from integration result for interpolator input.
template< typename TimeType, typename StateScalarType >
std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > convertNumericalSolutionToEphemerisInput(
        const int bodyIndex,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > integrationToEphemerisFrameFunction = NULL )
{
    std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > ephemerisTable;

    if( integrationToEphemerisFrameFunction == 0 )
    {
        for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator bodyIterator =
             equationsOfMotionNumericalSolution.begin( ); bodyIterator != equationsOfMotionNumericalSolution.end( ); bodyIterator++ )
        {
            ephemerisTable[ bodyIterator->first ] = bodyIterator->second.block( 6 * bodyIndex, 0, 6, 1 );
        }
    }
    else
    {
        for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator bodyIterator =
             equationsOfMotionNumericalSolution.begin( ); bodyIterator != equationsOfMotionNumericalSolution.end( ); bodyIterator++ )
        {

            ephemerisTable[ bodyIterator->first ] = bodyIterator->second.block( 6 * bodyIndex, 0, 6, 1 ) -
                    integrationToEphemerisFrameFunction( bodyIterator->first );
        }
    }
    return ephemerisTable;
}



template< typename TimeType, typename StateScalarType >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > > createStateInterpolator(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > >& stateMap );

//! Creates the interpolator for the ephemerides from the numerical integration results.
/*!
 *  Creates the interpolator for the ephemerides from the numerical integration results using
 *  variationalEquationsSolution_ and, if applicable, variationalDerivativesSolutions_.
 *  \return Vector of body state interpolators, order determined by bodiesToIntegrate vector.
 */
template< typename TimeType, typename StateScalarType >
void createAndSetInterpolatorsForEphemerides(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& ephemerisUpdateOrder,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >& integrationToEphemerisFrameFunctions =
        std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >( ) )
{
    using namespace tudat::interpolators;

    // Iterate over all bodies that are integrated numerically and create state interpolator.
    for( unsigned int i = 0; i < ephemerisUpdateOrder.size( ); i++ )
    {
        std::vector< std::string >::const_iterator bodyFindIterator = std::find(
                    bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), ephemerisUpdateOrder.at( i ) );
        int bodyIndex = std::distance( bodiesToIntegrate.begin( ), bodyFindIterator );

        boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > integrationToEphemerisFrameFunction = NULL;

        if( integrationToEphemerisFrameFunctions.count( bodiesToIntegrate.at( bodyIndex ) ) > 0 )
        {
            integrationToEphemerisFrameFunction = integrationToEphemerisFrameFunctions.at( bodiesToIntegrate.at( bodyIndex ) );
        }

        // Create interpolator.
        boost::shared_ptr< OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > > ephemerisInterpolator =
                createStateInterpolator( convertNumericalSolutionToEphemerisInput(
                                             bodyIndex, equationsOfMotionNumericalSolution, integrationToEphemerisFrameFunction ) );

        resetIntegratedEphemerisOfBody( bodyMap, ephemerisInterpolator, bodiesToIntegrate.at( bodyIndex ) );
    }
}

//! Reset ephemerides of numerically integrated bodies.
/*!
 *  Reset ephemerides of numerically integrated bodies, i.e. use numerical integration results to create new lookup table and interpolator.
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedEphemerides(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::vector< std::string >& bodiesToIntegrate,
        std::vector< std::string > ephemerisUpdateOrder = std::vector< std::string >( ),
        const std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >& integrationToEphemerisFrameFunctions =
        std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >( ) )
{
    if( ephemerisUpdateOrder.size( ) == 0 )
    {
        ephemerisUpdateOrder = bodiesToIntegrate;
    }

    // Create interpolators from numerical integration results (states) at discrete times.
    createAndSetInterpolatorsForEphemerides( bodyMap, bodiesToIntegrate, ephemerisUpdateOrder, equationsOfMotionNumericalSolution,
                                             integrationToEphemerisFrameFunctions );

    // Having set new ephemerides, update body properties depending on ephemerides.
    for( std::map< std::string, boost::shared_ptr< simulation_setup::Body > >::const_iterator bodyIterator = bodyMap.begin( );
         bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        bodyIterator->second->updateConstantEphemerisDependentMemberQuantities( );
    }
}


std::vector< std::string > determineEphemerisUpdateorder( std::vector< std::string > integratedBodies,
                                                          std::vector< std::string > centralBodies,
                                                          std::vector< std::string > ephemerisOrigins );




template< typename TimeType, typename StateScalarType >
class IntegratedStateProcessor
{
public:
    IntegratedStateProcessor( const IntegratedStateType stateType, const std::pair< int, int > startIndexAndSize ):
        stateType_( stateType ), startIndexAndSize_( startIndexAndSize ){ }

    virtual ~IntegratedStateProcessor( ){ }

    virtual void processIntegratedStates(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& numericalSolution ) = 0;

    virtual void processIntegratedMultiArcStates(
            const std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >& numericalSolution,
            const std::vector< std::pair< double, double > >& arcStartEndTimes ) = 0;

    IntegratedStateType stateType_;

    std::pair< int, int > startIndexAndSize_;
};

template< typename TimeType, typename StateScalarType >
class TranslationalStateIntegratedStateProcessor: public IntegratedStateProcessor< TimeType, StateScalarType >
{
public:
    TranslationalStateIntegratedStateProcessor(
            const int startIndex,
            const NamedBodyMap& bodyMap,
            const std::vector< std::string >& bodiesToIntegrate,
            const std::vector< std::string >& centralBodies,
            const boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager ):
        IntegratedStateProcessor< TimeType, StateScalarType >( transational_state, std::make_pair( startIndex, 6 * bodiesToIntegrate.size( ) ) ),
        bodyMap_( bodyMap ), bodiesToIntegrate_( bodiesToIntegrate )
    {
        ephemerisUpdateOrder_ = determineEphemerisUpdateorder(
                    bodiesToIntegrate_, centralBodies, frameManager->getEphemerisOrigins( bodiesToIntegrate ) );
        integrationToEphemerisFrameFunctions_ =
                ephemerides::getTranslationFunctionsFromIntegrationFrameToEphemerisFrame< StateScalarType, TimeType >(
                    centralBodies, bodiesToIntegrate_, frameManager );
    }

    void processIntegratedStates(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& numericalSolution )
    {
        resetIntegratedEphemerides< TimeType, StateScalarType >(
                    bodyMap_, numericalSolution, bodiesToIntegrate_, ephemerisUpdateOrder_,
                    integrationToEphemerisFrameFunctions_ );
    }

    void processIntegratedMultiArcStates(
            const std::vector< std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > > >& numericalSolution,
            const std::vector< std::pair< double, double > >& arcStartEndTimes )
    {
        resetMultiArcIntegratedEphemerides< TimeType, StateScalarType >(
                    bodyMap_, numericalSolution, arcStartEndTimes,
                    bodiesToIntegrate_, ephemerisUpdateOrder_, integrationToEphemerisFrameFunctions_ );
    }

private:
    NamedBodyMap bodyMap_;

    std::vector< std::string > bodiesToIntegrate_;

    std::vector< std::string > ephemerisUpdateOrder_;

    std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > > integrationToEphemerisFrameFunctions_;
};


void checkTranslationalStatesFeasibility(
        const std::vector< std::string >& bodiesToIntegrate,
        const NamedBodyMap& bodyMap );


template< typename TimeType, typename StateScalarType >
std::map< IntegratedStateType, std::vector< boost::shared_ptr< IntegratedStateProcessor< TimeType, StateScalarType > > > >
createIntegratedStateProcessors(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const NamedBodyMap& bodyMap,
        const boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager,
        const int startIndex = 0 )
{
    std::map< IntegratedStateType, std::vector< boost::shared_ptr< IntegratedStateProcessor< TimeType, StateScalarType > > > >
            integratedStateProcessors;
    switch( propagatorSettings->stateType_ )
    {
    case hybrid:
    {
        boost::shared_ptr< MultiTypePropagatorSettings< StateScalarType > > multiTypePropagatorSettings =
                boost::dynamic_pointer_cast< MultiTypePropagatorSettings< StateScalarType > >( propagatorSettings );
        std::map< IntegratedStateType, std::vector< boost::shared_ptr< IntegratedStateProcessor< TimeType, StateScalarType > > > >
                singleTypeIntegratedStateProcessors;
        int currentStartIndex = 0;
        for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< PropagatorSettings< StateScalarType > > > >::const_iterator
             typeIterator = multiTypePropagatorSettings->propagatorSettingsMap_.begin( );
             typeIterator != multiTypePropagatorSettings->propagatorSettingsMap_.end( ); typeIterator++ )
        {
            if( typeIterator->first != hybrid )
            {
                for( unsigned int i = 0; i < typeIterator->second.size( ); i++ )
                {
                    singleTypeIntegratedStateProcessors = createIntegratedStateProcessors< TimeType, StateScalarType >(
                                typeIterator->second.at( i ), bodyMap, frameManager, currentStartIndex );
                    if( singleTypeIntegratedStateProcessors.size( ) != 1 )
                    {
                        std::cerr<<"Error when making hybrid integrated result processors, multiple types found"<<std::endl;
                    }
                    else
                    {
                        if( singleTypeIntegratedStateProcessors.begin( )->second.size( ) != 1 )
                        {
                            std::cerr<<"Error when making hybrid integrated result processors, multiple processors of single type found"<<std::endl;
                        }
                    }
                    currentStartIndex += typeIterator->second.at( i )->getSingleArcStateSize( );
                    integratedStateProcessors[ singleTypeIntegratedStateProcessors.begin( )->first ].push_back(
                                singleTypeIntegratedStateProcessors.begin( )->second.at( 0 ) );
                }
            }
            else
            {
                std::cerr<<"Error when making integrated state processors, cannot handle hybrid propagator inside hybrid propagator"<<std::endl;
            }

        }

        break;
    }
    case transational_state:
    {
        boost::shared_ptr< TranslationalStatePropagatorSettings< StateScalarType > > translationalPropagatorSettings =
                boost::dynamic_pointer_cast< TranslationalStatePropagatorSettings< StateScalarType > >( propagatorSettings );

        checkTranslationalStatesFeasibility( translationalPropagatorSettings->bodiesToIntegrate_, bodyMap );

        integratedStateProcessors[ transational_state ].push_back(
                    boost::make_shared< TranslationalStateIntegratedStateProcessor< TimeType, StateScalarType > >(
                        startIndex, bodyMap, translationalPropagatorSettings->bodiesToIntegrate_,
                        translationalPropagatorSettings->centralBodies_, frameManager ) );
        break;
    }
    default:
        std::cerr<<"Error, could not process integrated state type "<<propagatorSettings->stateType_<<std::endl;
    }
    return integratedStateProcessors;
}

template< typename TimeType, typename StateScalarType >
void resetIntegratedStates(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::map< IntegratedStateType, std::vector< boost::shared_ptr< IntegratedStateProcessor< TimeType, StateScalarType > > > >
        integratedStateProcessors )
{
    for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr< IntegratedStateProcessor< TimeType, StateScalarType > > > >::
         const_iterator updateIterator = integratedStateProcessors.begin( ); updateIterator != integratedStateProcessors.end( ); updateIterator++ )
    {
        for( unsigned int i = 0; i < updateIterator->second.size( ); i++ )
        {
            updateIterator->second.at( i )->processIntegratedStates( equationsOfMotionNumericalSolution );
        }
    }
}


}

}

#endif // SETNUMERICALLYINTEGRATEDSTATES_H
