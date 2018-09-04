/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "Tudat/JsonInterface/Estimation/orbitDetermination.h"

namespace tudat
{

namespace simulation_setup
{

void updateInverseAPrioriCovarianceFromJSON(
        const nlohmann::json& jsonObject, const int numberOfParameters, Eigen::MatrixXd& inverseAprioriCovariance )
{
    using namespace json_interface;
    using K = Keys::Estimation;

    inverseAprioriCovariance.setZero( 0, 0 );
    try
    {
        Eigen::MatrixXd data = getValue< Eigen::MatrixXd >( jsonObject, K::inverseAprioriCovariance );

        if( data.cols( ) == 1 )
        {
            inverseAprioriCovariance = data.asDiagonal( );
        }
        else
        {
            inverseAprioriCovariance = data;
        }
    }
    catch( ... ) { }

    if( inverseAprioriCovariance.rows( ) > 0 )
    {
        if( ( inverseAprioriCovariance.rows( ) != numberOfParameters ) ||
                ( inverseAprioriCovariance.cols( ) != numberOfParameters ) )
        {
            throw std::runtime_error( "Error when reading a priori covariance from file, sizes are incompatible" );
        }
    }
    else
    {
        try
        {
            std::map< int, double > data = getValue< std::map< int, double > >( jsonObject, K::inverseAprioriCovariance );
            inverseAprioriCovariance.setZero( numberOfParameters, numberOfParameters );
            for( auto entry : data )
            {
                inverseAprioriCovariance( entry.first, entry.first ) = entry.second;
            }
        }
        catch( ... )
        {
            //            try
            //            {
            //                std::map< std::pair< int, int >, double > data =
            //                        getValue< std::map< std::pair< int, int >, double > >( jsonObject, K::inverseAprioriCovariance );

            //                for( auto entry : data )
            //                {
            //                    inverseAprioriCovariance( entry.first.first, entry.first.second ) = entry.second;
            //                }
            //            }
            //            catch( ... ){ }
        }
    }



}

void updateObservationWeightsFromJSON(
        const nlohmann::json& jsonObject,
        const std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, int > > numberOfObservations,
        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > >& observableWeights )
{
    using namespace json_interface;
    using K = Keys::Estimation;

    try
    {
        double constantWeight = getValue< double >( jsonObject, K::dataWeights );

        for( auto observable : numberOfObservations )
        {
            for( auto linkEnds : observable.second )
            {
                observableWeights[ observable.first ][ linkEnds.first ] = Eigen::VectorXd::Constant(
                            linkEnds.second, constantWeight );
            }
        }

        return;
    }
    catch( ... ){ }

//    try
//    {
//        std::map< observation_models::ObservableType, double > weightPerObservable = getValue<
//                std::map< observation_models::ObservableType, double > >( jsonObject, K::dataWeights );

//        for( auto observable : numberOfObservations )
//        {
//            for( auto linkEnds : observable.second )
//            {
//                observableWeights[ observable.first ][ linkEnds.first ] = Eigen::VectorXd::Constant(
//                            linkEnds.second, weightPerObservable.at( observable.first ) );
//            }
//        }


//        return;
//    }
//    catch( ... ){ }

//    try
//    {
//        std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, double > >
//                weightPerObservableAndLinkEnds = getValue<
//                std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, double > > >(
//                    jsonObject, K::dataWeights );

//        for( auto observable : numberOfObservations )
//        {
//            for( auto linkEnds : observable.second )
//            {
//                observableWeights[ observable.first ][ linkEnds.first ] = Eigen::VectorXd::Constant(
//                            linkEnds.second, weightPerObservableAndLinkEnds.at( observable.first ).at( linkEnds.first ) );
//            }
//        }

//        return;
//    }
//    catch( ... ){ }

//    try
//    {
//        observableWeights = getValue<
//                std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > > >(
//                    jsonObject, K::dataWeights );

//        return;
//    }
//    catch( ... ){ }
}


} // namespace propagators

} // namespace tudat
