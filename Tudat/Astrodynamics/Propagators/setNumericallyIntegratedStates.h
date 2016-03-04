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

//! Function to reset the tabulated ephemeris of a body
/*!
 * Function to reset the tabulated ephemeris of a body, this requires the requested body to possess an ephemeris
 * of type  TabulatedCartesianEphemeris< StateScalarType, TimeType >
 * \param bodyMap List of bodies used in simulations.
 * \param ephemerisInterpolator Interpolator providing the new state of the body as a function of time.
 * \param bodyToIntegrate Name of body for which the ephemeris is to be reset.
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedEphemerisOfBody(
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< interpolators::OneDimensionalInterpolator<
        TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > > ephemerisInterpolator,
        const std::string& bodyToIntegrate )
{
    using namespace tudat::interpolators;
    using namespace tudat::ephemerides;

    // If body does not have ephemeris, give error message.
    if( bodyMap.at( bodyToIntegrate )->getEphemeris( ) == NULL )
    {
        throw std::runtime_error( "Error when resetting integrated ephemeris of body " +
                                  bodyToIntegrate + "no ephemeris found" );
    }

    // If current ephemeris is not already a tabulated ephemeris, give error message.
    else if( boost::dynamic_pointer_cast< TabulatedCartesianEphemeris< StateScalarType, TimeType > >(
                 bodyMap.at( bodyToIntegrate )->getEphemeris( ) ) == NULL )
    {
        throw std::runtime_error( "Error when resetting integrated ephemeris of body " +
                                  bodyToIntegrate + "no tabulated ephemeris found" );

    }
    // Else, update existing tabulated ephemeris
    else
    {

        boost::shared_ptr< TabulatedCartesianEphemeris< StateScalarType, TimeType > > tabulatedEphemeris =
                boost::dynamic_pointer_cast< TabulatedCartesianEphemeris<  StateScalarType, TimeType > >(
                    bodyMap.at( bodyToIntegrate )->getEphemeris( ) );
        tabulatedEphemeris->resetInterpolator( ephemerisInterpolator );
    }
}

//! Function to convert output of translational motion from the numerical integrator to the required input for the ephemeris.
/*!
* Function to convert output of translational motion from the numerical integrator to the required input for the ephemeris.
* It extracts the state history of a single body from the full list of integrated states
* Additionally, it changes the origin of the reference frame in which the states are given, by using the
* integrationToEphemerisFrameFunction input variable.
* \param bodyIndex Index of integrated body for which the state is to be retrieved
* \param equationsOfMotionNumericalSolution Full numerical solution of numerical integrator, already converted to Cartesian
* states (w.r.t. the integration origin of the body of bodyIndex)
* \param integrationToEphemerisFrameFunction Function to provide the state of the ephemeris origin of the current body
* w.r.t. its integration origin.
* \return State history of body bodyIndex w.r.t. the origin with which its ephemeris is defined.
*/
template< typename TimeType, typename StateScalarType >
std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > convertNumericalSolutionToEphemerisInput(
        const int bodyIndex,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) >
        integrationToEphemerisFrameFunction = NULL )
{
    std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > ephemerisTable;

    // If no integrationToEphemerisFrameFunction is provided, origin is already correct; only extract required indices.
    if( integrationToEphemerisFrameFunction == 0 )
    {
        for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
             bodyIterator = equationsOfMotionNumericalSolution.begin( );
             bodyIterator != equationsOfMotionNumericalSolution.end( ); bodyIterator++ )
        {
            ephemerisTable[ bodyIterator->first ] = bodyIterator->second.block( 6 * bodyIndex, 0, 6, 1 );
        }
    }
    // Else, extract indices and add required translation from integrationToEphemerisFrameFunction
    else
    {
        for( typename std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >::const_iterator
             bodyIterator = equationsOfMotionNumericalSolution.begin( );
             bodyIterator != equationsOfMotionNumericalSolution.end( ); bodyIterator++ )
        {

            ephemerisTable[ bodyIterator->first ] = bodyIterator->second.block( 6 * bodyIndex, 0, 6, 1 ) -
                    integrationToEphemerisFrameFunction( bodyIterator->first );
        }
    }
    return ephemerisTable;
}

//! Function to create an interpolator for the new translational state of a body.
/*!
 * Function to create an interpolator for the new translational state of a body.
 * \param stateMap New state history, w.r.t. the required ephemeris origin.
 * \return Lagrange interpolator (order 6) that produces the required continuous state.
 */
template< typename TimeType, typename StateScalarType >
boost::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >
createStateInterpolator(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > >& stateMap );

//! Creates and resets the interpolator for the ephemerides of the integrated bodies from the numerical integration results.
/*!
 * Creates and resets the interpolator for the ephemerides of the integrated bodies from the numerical integration results.
 * \param bodyMap List of bodies used in simulations.
 * \param bodiesToIntegrate List of names of bodies which are numericall integrated (in the order in which they are
 * in the equationsOfMotionNumericalSolution map.
 * \param ephemerisUpdateOrder Order in which to update the ephemeris objects.
 * \param equationsOfMotionNumericalSolution Numerical solution of translational equations of motion, in Cartesian elements
 * w.r.t. integratation origins.
 * \param integrationToEphemerisFrameFunctions Function to provide the states of the ephemeris origins of each body
 * w.r.t. their respective integration origins.
 */
template< typename TimeType, typename StateScalarType >
void createAndSetInterpolatorsForEphemerides(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::vector< std::string >& bodiesToIntegrate,
        const std::vector< std::string >& ephemerisUpdateOrder,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >&
        integrationToEphemerisFrameFunctions =
        std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >( ) )
{
    using namespace tudat::interpolators;

    // Iterate over all bodies that are integrated numerically and create state interpolator.
    for( unsigned int i = 0; i < ephemerisUpdateOrder.size( ); i++ )
    {
        // Get index of current body to be updated in bodiesToIntegrate.
        std::vector< std::string >::const_iterator bodyFindIterator = std::find(
                    bodiesToIntegrate.begin( ), bodiesToIntegrate.end( ), ephemerisUpdateOrder.at( i ) );
        if( bodyFindIterator == bodiesToIntegrate.end( ) )
        {
            throw std::runtime_error( "Error when creating and setting ephemeris after integration, cannot find body " +
                                      ephemerisUpdateOrder.at( i ) );
        }
        int bodyIndex = std::distance( bodiesToIntegrate.begin( ), bodyFindIterator );

        // Get frame origin function if applicable
        boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > integrationToEphemerisFrameFunction =
                NULL;
        if( integrationToEphemerisFrameFunctions.count( bodiesToIntegrate.at( bodyIndex ) ) > 0 )
        {
            integrationToEphemerisFrameFunction =
                    integrationToEphemerisFrameFunctions.at( bodiesToIntegrate.at( bodyIndex ) );
        }

        // Create and reset interpolator.
        boost::shared_ptr< OneDimensionalInterpolator< TimeType, Eigen::Matrix< StateScalarType, 6, 1 > > >
                ephemerisInterpolator =
                createStateInterpolator( convertNumericalSolutionToEphemerisInput(
                                             bodyIndex, equationsOfMotionNumericalSolution,
                                             integrationToEphemerisFrameFunction ) );
        resetIntegratedEphemerisOfBody( bodyMap, ephemerisInterpolator, bodiesToIntegrate.at( bodyIndex ) );
    }
}

//! Resets the ephemerides of the integrated bodies from the numerical integration results.
/*!
 * Resets the ephemerides of the integrated bodies from the numerical integration results, and performs associated
 * computation for ephemeris-dependent environment variables.
 * \param bodyMap List of bodies used in simulations.
 * \param bodiesToIntegrate List of names of bodies which are numericall integrated (in the order in which they are
 * in the equationsOfMotionNumericalSolution map.
 * \param ephemerisUpdateOrder Order in which to update the ephemeris objects (empty if arbitrary).
 * \param equationsOfMotionNumericalSolution Numerical solution of translational equations of motion, in Cartesian elements
 * w.r.t. integratation origins.
 * \param integrationToEphemerisFrameFunctions Function to provide the states of the ephemeris origins of each body
 * w.r.t. their respective integration origins.
 */
template< typename TimeType, typename StateScalarType >
void resetIntegratedEphemerides(
        const simulation_setup::NamedBodyMap& bodyMap,
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::vector< std::string >& bodiesToIntegrate,
        std::vector< std::string > ephemerisUpdateOrder = std::vector< std::string >( ),
        const std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >&
        integrationToEphemerisFrameFunctions = std::map< std::string,
        boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >( ) )
{
    // Set update order arbitrarily if no order is provided.
    if( ephemerisUpdateOrder.size( ) == 0 )
    {
        ephemerisUpdateOrder = bodiesToIntegrate;
    }
    // Check input consistency
    else if( ephemerisUpdateOrder.size( ) != bodiesToIntegrate.size( ) )
    {
        throw std::runtime_error( "Error when resetting ephemerides, input vectors have inconsistent size" );
    }

    // Create interpolators from numerical integration results (states) at discrete times.
    createAndSetInterpolatorsForEphemerides(
                bodyMap, bodiesToIntegrate, ephemerisUpdateOrder, equationsOfMotionNumericalSolution,
                                             integrationToEphemerisFrameFunctions );

    // Having set new ephemerides, update body properties depending on ephemerides.
    for( std::map< std::string, boost::shared_ptr< simulation_setup::Body > >::const_iterator
         bodyIterator = bodyMap.begin( );
         bodyIterator != bodyMap.end( ); bodyIterator++ )
    {
        //throw std::runtime_error( "Update of ephemeris dependent quantities not yet implemented" );
        //bodyIterator->second->updateConstantEphemerisDependentMemberQuantities( );
    }
}

std::vector< std::string > determineEphemerisUpdateorder( std::vector< std::string > integratedBodies,
                                                          std::vector< std::string > centralBodies,
                                                          std::vector< std::string > ephemerisOrigins );

//! Base class for defining settings on how numerically integrated states are to be processed in the environment
/*!
 *  Base class for defining settings on how numerically integrated states are to be processed in the environment.
 *  This class is pure virtual, and a derived class must be provided for each state type.
 */
template< typename TimeType, typename StateScalarType >
class IntegratedStateProcessor
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param stateType Type of state that is to be set in environment.
     * \param startIndexAndSize Start index of current state type in integrated state vector, and its size in the vector
     * (first and second entry of pair).
     */
    IntegratedStateProcessor( const IntegratedStateType stateType, const std::pair< int, int > startIndexAndSize ):
        stateType_( stateType ), startIndexAndSize_( startIndexAndSize ){ }

    //! Virtual destructor.
    virtual ~IntegratedStateProcessor( ){ }

    //! Function that processes the entries of the stateType_ in the full numericalSolution
    /*!
     * Function that processes the entries of the stateType_ in the full numericalSolution
     * \param numericalSolution Full numerical solution, in global representation (see
     * convertToOutputSolution function in associated SingleStateTypeDerivative derived class.
     */
    virtual void processIntegratedStates(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& numericalSolution ) = 0;

    //! Type of state that is to be set in environment.
    IntegratedStateType stateType_;

    //! Start index of current state type in integrated state vector, and its size in the vector
    /*!
     *  Start index of current state type in integrated state vector, and its size in the vector
     * (first and second entry of pair).
     */
    std::pair< int, int > startIndexAndSize_;
};

//! Class used for processing numerically integrated translational states
/*!
 *  Class used for processing numerically integrated translational states, updates ephemeris object of each integrated
 *  body, performing frame translations if needed.
 */
template< typename TimeType, typename StateScalarType >
class TranslationalStateIntegratedStateProcessor: public IntegratedStateProcessor< TimeType, StateScalarType >
{
public:
    //! Constructor.
    /*!
     * Constructor.
     * \param startIndex Index in the state vector where the translational state starts.
     * \param bodyMap List of bodies used in simulations.
     * \param bodiesToIntegrate List of bodies for which the translational state is numerically integrated. Order in this
     * vector is the same as the order in state vector.
     * \param centralBodies List of origing w.r.t. which the translational states are integrated (with the same order
     * as bodiesToIntegrate).
     * \param frameManager Object to get state of one body w.r.t. another body.
     */
    TranslationalStateIntegratedStateProcessor(
            const int startIndex,
            const simulation_setup::NamedBodyMap& bodyMap,
            const std::vector< std::string >& bodiesToIntegrate,
            const std::vector< std::string >& centralBodies,
            const boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager ):
        IntegratedStateProcessor< TimeType, StateScalarType >(
            transational_state, std::make_pair( startIndex, 6 * bodiesToIntegrate.size( ) ) ),
        bodyMap_( bodyMap ), bodiesToIntegrate_( bodiesToIntegrate )
    {
        // Get update orders.
        ephemerisUpdateOrder_ = determineEphemerisUpdateorder(
                    bodiesToIntegrate_, centralBodies, frameManager->getEphemerisOrigins( bodiesToIntegrate ) );

        // Get required frame origin translations.
        integrationToEphemerisFrameFunctions_ =
                ephemerides::getTranslationFunctionsFromIntegrationFrameToEphemerisFrame< StateScalarType, TimeType >(
                    centralBodies, bodiesToIntegrate_, frameManager );
    }

    //! Function that processes the entries of the translational state in the full numericalSolution
    /*!
     * Function that processes the entries of the translational state in the full numericalSolution, extracts and converts
     * the states to the required frames, and updates the associated ephemerides.
     * \param numericalSolution Full numerical solution, in global representation (see
     * convertToOutputSolution function in NBodyStateDerivative class.
     */
    void processIntegratedStates(
            const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& numericalSolution )
    {
        resetIntegratedEphemerides< TimeType, StateScalarType >(
                    bodyMap_, numericalSolution, bodiesToIntegrate_, ephemerisUpdateOrder_,
                    integrationToEphemerisFrameFunctions_ );
    }

private:

    //! List of bodies used in simulations.
    simulation_setup::NamedBodyMap bodyMap_;

    //! List of bodies for which the translational state is numerically integrated.
    /*!
     * List of bodies for which the translational state is numerically integrated. Order in this vector is the same
     * as the order in state vector.
     */
    std::vector< std::string > bodiesToIntegrate_;

    //! Order in which to update the ephemeris objects
    std::vector< std::string > ephemerisUpdateOrder_;

    //! Function to provide the states of the ephemeris origins of each body w.r.t. their respective integration origins.
    std::map< std::string, boost::function< Eigen::Matrix< StateScalarType, 6, 1 >( const TimeType ) > >
    integrationToEphemerisFrameFunctions_;
};


void checkTranslationalStatesFeasibility(
        const std::vector< std::string >& bodiesToIntegrate,
        const simulation_setup::NamedBodyMap& bodyMap );


template< typename TimeType, typename StateScalarType >
std::map< IntegratedStateType, std::vector< boost::shared_ptr< IntegratedStateProcessor< TimeType, StateScalarType > > > >
createIntegratedStateProcessors(
        const boost::shared_ptr< PropagatorSettings< StateScalarType > > propagatorSettings,
        const simulation_setup::NamedBodyMap& bodyMap,
        const boost::shared_ptr< ephemerides::ReferenceFrameManager > frameManager,
        const int startIndex = 0 )
{
    std::map< IntegratedStateType, std::vector< boost::shared_ptr< IntegratedStateProcessor
            < TimeType, StateScalarType > > > > integratedStateProcessors;
    switch( propagatorSettings->stateType_ )
    {   
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
        throw std::runtime_error( "Error, could not process integrated state type " +
                                  boost::lexical_cast< std::string >( propagatorSettings->stateType_ ) );
    }
    return integratedStateProcessors;
}

template< typename TimeType, typename StateScalarType >
void resetIntegratedStates(
        const std::map< TimeType, Eigen::Matrix< StateScalarType, Eigen::Dynamic, 1 > >& equationsOfMotionNumericalSolution,
        const std::map< IntegratedStateType, std::vector< boost::shared_ptr
        < IntegratedStateProcessor< TimeType, StateScalarType > > > >
        integratedStateProcessors )
{
    for( typename std::map< IntegratedStateType, std::vector< boost::shared_ptr
         < IntegratedStateProcessor< TimeType, StateScalarType > > > >::
         const_iterator updateIterator = integratedStateProcessors.begin( );
         updateIterator != integratedStateProcessors.end( ); updateIterator++ )
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
