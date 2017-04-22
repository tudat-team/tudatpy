/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_OBSERVATIONBIASPARAMETER_H
#define TUDAT_OBSERVATIONBIASPARAMETER_H

#include <iostream>

#include "Tudat/Astrodynamics/ObservationModels/observationBias.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"

namespace tudat
{

namespace estimatable_parameters
{

class ConstantObservationBiasParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     * Constructor
     */
    ConstantObservationBiasParameter(
            const boost::function< Eigen::VectorXd( ) > getCurrentBias,
            const boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
            const observation_models::LinkEnds linkEnds,
            const observation_models::ObservableType observableType ):
        EstimatableParameter< Eigen::VectorXd >( constant_additive_observation_bias, linkEnds.begin( )->second. first ),
        getCurrentBias_( getCurrentBias ), resetCurrentBias_( resetCurrentBias ),
        linkEnds_( linkEnds ), observableType_( observableType ){ }

    //! Destructor
    ~ConstantObservationBiasParameter( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        return getCurrentBias_( );
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        resetCurrentBias_( parameterValue );
    }


    int getParameterSize( )
    {

        return observation_models::getObservableSize( observableType_ );
    }


    void setObservationBiasFunctions(
            const boost::function< Eigen::VectorXd( ) > getCurrentBias,
            const boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias )
    {
        if( !getCurrentBias.empty( ) || !resetCurrentBias_.empty( ) )
        {
            std::cerr<<"Warning when resetting relative observation bias in estimation object, existing contents not empty"<<std::endl;
        }

        getCurrentBias_ = getCurrentBias;
        resetCurrentBias_ = resetCurrentBias;
    }

    observation_models::LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

protected:

private:


    boost::function< Eigen::VectorXd( ) > getCurrentBias_;

    boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias_;

    observation_models::LinkEnds linkEnds_;

    observation_models::ObservableType observableType_;

};

class ConstantRelativeObservationBiasParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     * Constructor
     */
    ConstantRelativeObservationBiasParameter(
            const boost::function< Eigen::VectorXd( ) > getCurrentBias,
            const boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
            const observation_models::LinkEnds linkEnds,
            const observation_models::ObservableType observableType ):
        EstimatableParameter< Eigen::VectorXd >( constant_relative_observation_bias, linkEnds.begin( )->second.first ),
        getCurrentBias_( getCurrentBias ), resetCurrentBias_( resetCurrentBias ), linkEnds_( linkEnds ), observableType_( observableType ){ }

    //! Destructor
    ~ConstantRelativeObservationBiasParameter( ) { }

    Eigen::VectorXd getParameterValue( )
    {
        return getCurrentBias_( );
    }

    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        resetCurrentBias_( parameterValue );
    }

    int getParameterSize( )
    {
        return observation_models::getObservableSize( observableType_ );
    }

    void setObservationBiasFunctions(
            const boost::function< Eigen::VectorXd( ) > getCurrentBias,
            const boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias )
    {
        if( !getCurrentBias.empty( ) || !resetCurrentBias_.empty( ) )
        {
            std::cerr<<"Warning when resetting relative observation bias in estimation object, existing contents not empty"<<std::endl;
        }

        getCurrentBias_ = getCurrentBias;
        resetCurrentBias_ = resetCurrentBias;
    }

    observation_models::LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }

    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }


protected:

private:

    boost::function< Eigen::VectorXd( ) > getCurrentBias_;

    boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias_;

    observation_models::LinkEnds linkEnds_;

    observation_models::ObservableType observableType_;

};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_OBSERVATIONBIASPARAMETER_H
