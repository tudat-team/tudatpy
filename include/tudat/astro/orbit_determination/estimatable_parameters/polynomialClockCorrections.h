#ifndef POLYNOMIALCLOCKCORRECTIONS_H
#define POLYNOMIALCLOCKCORRECTIONS_H

#include "tudat/astro/system_models/timingSystem.h"
#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"


namespace tudat
{

namespace estimatable_parameters
{

class GlobalPolynomialClockCorrections: public EstimatableParameter< Eigen::VectorXd >
{

public:

    GlobalPolynomialClockCorrections( const std::shared_ptr< system_models::TimingSystem > timingSystem,
                                      const std::vector< int > correctionPowers,
                                      const std::string& associatedBody,
                                      const std::string& associatedStation ):
        EstimatableParameter< Eigen::VectorXd  >( global_polynomial_clock_corrections, associatedBody, associatedStation ),
        timingSystem_( timingSystem ), correctionPowers_( correctionPowers ){ }

    //! Destructor
    /*!
     *  Destructor
     */
    ~GlobalPolynomialClockCorrections( ) { }

    Eigen::VectorXd  getParameterValue( )
    {
        Eigen::VectorXd currentParameterValues = Eigen::VectorXd::Zero( getParameterSize( ) );
        for( unsigned int i = 0; i < correctionPowers_.size( ); i++ )
        {
            currentParameterValues[ i ] = timingSystem_->getMeanPolynomialDriftCoefficients( correctionPowers_[ i ] );
        }
        return currentParameterValues;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        std::map< int, double > newCorrections;
        for( unsigned int i = 0; i < correctionPowers_.size( ); i++ )
        {
            newCorrections[ correctionPowers_[ i ] ] = parameterValue[ i ];
        }
        timingSystem_->setGlobalPolynomialClockCorrections( newCorrections );
    }

    int getParameterSize( )
    {
        return correctionPowers_.size( );
    }

    std::vector< int > getCorrectionPowers( )
    {
        return correctionPowers_;
    }

    std::shared_ptr< system_models::TimingSystem > getTimingSystem( )
    {
        return timingSystem_;
    }

    observation_models::LinkEndId getLinkEndId( )
    {
        return observation_models::LinkEndId( parameterName_.second );
    }



protected:

private:

    std::shared_ptr< system_models::TimingSystem > timingSystem_;

    std::vector< int > correctionPowers_;
};

class MultiArcClockCorrections: public EstimatableParameter< Eigen::VectorXd >
{

public:

    MultiArcClockCorrections( const std::shared_ptr< system_models::TimingSystem > timingSystem,
                              const std::vector< int > correctionPowers,
                              const std::vector< int > arcIndices,
                              const std::string& associatedBody,
                              const std::string& associatedStation ):
        EstimatableParameter< Eigen::VectorXd  >( arc_wise_polynomial_clock_corrections, associatedBody, associatedStation ),
        timingSystem_( timingSystem ), correctionPowers_( correctionPowers ), arcIndices_( arcIndices )
    { }

    //! Destructor
    /*!
     *  Destructor
     */
    ~MultiArcClockCorrections( ) { }

    Eigen::VectorXd  getParameterValue( )
    {
        Eigen::VectorXd currentParameterValues = Eigen::VectorXd::Zero( getParameterSize( ) );

        int counter = 0;
        for( unsigned int i = 0; i < arcIndices_.size( ); i++ )
        {
            for( unsigned int j = 0; j < correctionPowers_.size( ); j++ )
            {
                currentParameterValues[ counter ] = timingSystem_->getPolynomialDriftCoefficients( correctionPowers_[ j ], arcIndices_[ i ] );
                counter++;
            }
        }
        return currentParameterValues;
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        std::map< int, double > newCorrections;

        int counter = 0;
        for( unsigned int i = 0; i < arcIndices_.size( ); i++ )
        {
            newCorrections.clear( );
            for( unsigned int j = 0; j < correctionPowers_.size( ); j++ )
            {
                newCorrections[ correctionPowers_[ j ] ] = parameterValue[ counter ];
                counter++;
            }
            timingSystem_->setSingleArcPolynomialClockCorrections( newCorrections, arcIndices_[ i ] );
        }
    }

    int getParameterSize( )
    {
        return correctionPowers_.size( ) * arcIndices_.size( );
    }

    std::shared_ptr< system_models::TimingSystem > getTimingSystem( )
    {
        return timingSystem_;
    }

    std::vector< int > getCorrectionPowers( )
    {
        return correctionPowers_;
    }

    std::vector< int > getArcIndices( )
    {
        return arcIndices_;
    }

    observation_models::LinkEndId getLinkEndId( )
    {
        return observation_models::LinkEndId( parameterName_.second );
    }

protected:

private:

    std::shared_ptr< system_models::TimingSystem > timingSystem_;

    std::vector< int > correctionPowers_;

    std::vector< int > arcIndices_;
};

}

}

#endif // POLYNOMIALCLOCKCORRECTIONS_H
