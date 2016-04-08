#ifndef INITIALTRANSLATIONALSTATE_H
#define INITIALTRANSLATIONALSTATE_H

#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

template< typename InitialStateParameterType = double >
class InitialTranslationalStateParameter: public EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >
{
public:

    InitialTranslationalStateParameter(
            const std::string& associatedBody,
            const Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >& initialTranslationalState,
            const std::string& centralBody = "SSB", const std::string& frameOrientation = "ECLIPJ2000" ):
        EstimatableParameter< Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > >( initial_body_state, associatedBody ),
        initialTranslationalState_( initialTranslationalState ), centralBody_( centralBody ), frameOrientation_( frameOrientation )
    { }


    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > getParameterValue( )
    {
        return initialTranslationalState_;
    }

    void setParameterValue( Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 >  parameterValue )
    {
        initialTranslationalState_ = parameterValue;
    }

    int getParameterSize( )
    {
        return 6;
    }

    std::string getCentralBody( )
    {
        return centralBody_;
    }


private:

    Eigen::Matrix< InitialStateParameterType, Eigen::Dynamic, 1 > initialTranslationalState_;

    std::string centralBody_;

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

#endif // INITIALTRANSLATIONALSTATE_H
