#ifndef TUDAT_INITIALTRANSLATIONALSTATE_H
#define TUDAT_INITIALTRANSLATIONALSTATE_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of an initial translational state.
template< typename InitialStateParameterType = double >
class InitialTranslationalStateParameter: public EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >
{
public:

    //! Constructor, sets initial value of translational state.
    /*!
     * Constructor, sets initial value of translational state.
     * \param associatedBody Body for which initial state is to be estimated.
     * \param initialTranslationalState Current value of initial state (w.r.t. centralBody)
     * \param centralBody Body w.r.t. which the initial state is to be estimated.
     * \param frameOrientation Orientation of the frame in which the state is defined.
     */
    InitialTranslationalStateParameter(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& initialTranslationalState,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >( initial_body_state, associatedBody ),
        initialTranslationalState_( initialTranslationalState ), centralBody_( centralBody ), frameOrientation_( frameOrientation )
    { }

    //! Function to get the current value of initial state w.r.t. centralBody.
    /*!
     * Function to get the current value of initial state w.r.t. centralBody.
     * \return The current value of initial state w.r.t. centralBody.
     */
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getParameterValue( )
    {
        return initialTranslationalState_;
    }

    //! Function to reset the current value of initial state w.r.t. centralBody.
    /*!
     * Function to reset the current value of initial state w.r.t. centralBody.
     * \param parameterValue The new value of initial state w.r.t. centralBody.
     */
    void setParameterValue( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  parameterValue )
    {
        initialTranslationalState_ = parameterValue;
    }

    //! Function to retrieve the size of the parameter (always 6).
    /*!
     *  Function to retrieve the size of the parameter (always 6).
     *  \return Size of parameter value (always 6).
     */
    int getParameterSize( )
    {
        return 6;
    }

    //! Function to get the name of the body w.r.t. which the initial state is to be estimated.
    /*!
     * Function to get the name of the body w.r.t. which the initial state is to be estimated.
     * \return Name of the body w.r.t. which the initial state is to be estimated.
     */
    std::string getCentralBody( )
    {
        return centralBody_;
    }

private:

    //! Current value of initial state (w.r.t. centralBody)
    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialTranslationalState_;

    //!  Body w.r.t. which the initial state is to be estimated.
    std::string centralBody_;

    //! Orientation of the frame in which the state is defined.
    std::string frameOrientation_;

};

template< typename InitialStateParameterType = double >
int getSingleArcParameterSetSize(
        boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameterSet )
{
    int totalParameterSetSize = estimatableParameterSet->getEstimatedParameterSetSize( );
    std::vector< boost::shared_ptr< EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > > > >
            initialStateParameters = estimatableParameterSet->getEstimatedInitialStateParameters( );

    for( unsigned int i = 0; i < initialStateParameters.size( ); i++ )
    {
        if( ( initialStateParameters.at( i )->getParameterName( ).first != initial_body_state ) )
        {
            throw std::runtime_error( "Error when getting single arc paramater vector, did not recognize initial state parameter " +
                        boost::lexical_cast< std::string >( initialStateParameters.at( i )->getParameterName( ).first ) );
        }
    }
    return totalParameterSetSize;
}

template< typename InitialStateParameterType = double >
int getSingleArcInitialDynamicalStateParameterSetSize(
        boost::shared_ptr< EstimatableParameterSet< InitialStateParameterType > > estimatableParameterSet )
{
    return getSingleArcParameterSetSize( estimatableParameterSet ) -
            ( estimatableParameterSet->getEstimatedParameterSetSize( ) - estimatableParameterSet->getInitialDynamicalStateParameterSize( ) );
}


}

}

#endif // TUDAT_INITIALTRANSLATIONALSTATE_H
