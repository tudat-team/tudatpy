#ifndef CLOCKPARAMETERPARTIALS_H
#define CLOCKPARAMETERPARTIALS_H

#include <cmath>

#include "tudat/astro/orbit_determination/observation_partials/observationPartial.h"
#include "tudat/astro/system_models/timingSystem.h"

namespace tudat
{

namespace observation_partials
{

class TimingPartial
{
public:
    virtual ~TimingPartial( ){ }

    virtual Eigen::Matrix< double, 1, Eigen::Dynamic > getPartialOfClockErrorWrtParameter( const double time ) = 0;

protected:

private:

};

class RangePartialWrtClockParameter: public ObservationPartial< 1 >
{
public:
    typedef std::vector< std::pair< Eigen::Matrix< double, 1, Eigen::Dynamic >, double > > RangePartialReturnType;

    RangePartialWrtClockParameter(
            const estimatable_parameters::EstimatebleParameterIdentifier parameterIdentifier,
            const std::map< int, std::shared_ptr< TimingPartial > >& timingPartials,
            const std::map< int, double > multiplierMap ): ObservationPartial< 1 >( parameterIdentifier )
    {
        for( std::map< int, std::shared_ptr< TimingPartial > >::const_iterator partialIterator = timingPartials.begin( );
             partialIterator != timingPartials.end( ); partialIterator++ )
        {
            if( multiplierMap.count( partialIterator->first ) == 0 )
            {
                std::cerr<<"Error when making range partial wrt clock parameter, multiplier not found"<<std::endl;
            }
            else
            {
                timingPartials_[ partialIterator->first ] = std::make_pair(
                            partialIterator->second, multiplierMap.at( partialIterator->first ) );
            }
        }
    }

    ~RangePartialWrtClockParameter( ){ }

    RangePartialReturnType calculatePartial(
            const std::vector< Eigen::Vector6d >& states,
            const std::vector< double >& times,
            const observation_models::LinkEndType linkEndOfFixedTime = observation_models::receiver,
            const std::shared_ptr< observation_models::ObservationAncilliarySimulationSettings > ancillarySettings = nullptr,
            const Eigen::Matrix< double, 1, 1 >& currentObservation =
            Eigen::Matrix< double, 1, 1 >::Constant( TUDAT_NAN )  )
    {
        RangePartialReturnType partialSet;

        double currentTime;
        for( timingPartialIterator_ = timingPartials_.begin( ); timingPartialIterator_ != timingPartials_.end( ); timingPartialIterator_++ )
        {
            currentTime = times.at( timingPartialIterator_->first );
            partialSet.push_back( std::make_pair(
                                      physical_constants::SPEED_OF_LIGHT * timingPartialIterator_->second.second *
                                      timingPartialIterator_->second.first->getPartialOfClockErrorWrtParameter( currentTime ), currentTime ) );
        }

        return partialSet;
    }

protected:

private:
    std::map< int, std::pair< std::shared_ptr< TimingPartial >, double > > timingPartials_;

    std::map< int, std::pair< std::shared_ptr< TimingPartial >, double > >::iterator timingPartialIterator_;

};

class TimingPartialWrtGlobalPolynomialCorrectionCoefficients: public TimingPartial
{
public:
    TimingPartialWrtGlobalPolynomialCorrectionCoefficients(
            const std::shared_ptr< system_models::TimingSystem > timingSystem,
            const std::vector< int > polynomialPowers ): timingSystem_( timingSystem ), polynomialPowers_( polynomialPowers ){ }

    Eigen::Matrix< double, 1, Eigen::Dynamic > getPartialOfClockErrorWrtParameter( const double time );

private:

    std::shared_ptr< system_models::TimingSystem > timingSystem_;

    std::vector< int > polynomialPowers_;
};

class TimingPartialWrtArcwisePolynomialCorrectionCoefficients: public TimingPartial
{
public:
    TimingPartialWrtArcwisePolynomialCorrectionCoefficients(
            const std::shared_ptr< system_models::TimingSystem > timingSystem,
            const std::vector< int > polynomialPowers,
            const std::vector< int > arcIndices ): timingSystem_( timingSystem ), polynomialPowers_( polynomialPowers ), arcIndices_( arcIndices )
    {
        numberOfArcs_ = arcIndices_.size( );
    }

    ~TimingPartialWrtArcwisePolynomialCorrectionCoefficients( ){ }

    Eigen::Matrix< double, 1, Eigen::Dynamic > getPartialOfClockErrorWrtParameter( const double time );

private:

    std::shared_ptr< system_models::TimingSystem > timingSystem_;

    std::vector< int > polynomialPowers_;

    std::vector< int > arcIndices_;

    std::vector< int >::iterator arcFindIterator_;

    int numberOfArcs_;
};

}

}

#endif // CLOCKPARAMETERPARTIALS_H
