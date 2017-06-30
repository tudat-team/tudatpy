/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_PODINPUTOUTPUTTYPES_H
#define TUDAT_PODINPUTOUTPUTTYPES_H

#include <map>
#include <vector>

#include <boost/shared_ptr.hpp>

#include <Eigen/Core>

#include "Tudat/Astrodynamics/ObservationModels/linkTypeDefs.h"
#include "Tudat/Astrodynamics/ObservationModels/observableTypes.h"

namespace tudat
{

namespace simulation_setup
{

//! Data structure used to provide input to orbit determination procedure
template< typename ObservationScalarType = double, typename TimeType = double >
class PodInput
{

public:

    //! Typedef of vector of observations of a single type
    typedef Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > ObservationVectorType;

    //! Typedef of map of pairs (observation vector and pair (vector of associated times reference link end)), link ends
    //! as map key
    typedef std::map< observation_models::LinkEnds, std::pair< ObservationVectorType,
    std::pair< std::vector< TimeType >, observation_models::LinkEndType > > > SingleObservablePodInputType;

        //! List of SingleObservablePodInputType per observation type
    typedef std::map< observation_models::ObservableType, SingleObservablePodInputType > PodInputDataType;

    //! Constructor
    /*!
     * Constructor
     * \param observationsAndTimes Total data structure of observations and associated times/link ends/type
     * \param numberOfEstimatedParameters Size of vector of estimated parameters
     * \param inverseOfAprioriCovariance A priori covariance matrix (unnormalized) of estimated parameters. None (matrix of
     * size 0) by default
     * \param initialParameterDeviationEstimate Correction to estimated parameter vector to be applied on first iteration.
     * None (vector of size 0) by default
     */
    PodInput( const PodInputDataType& observationsAndTimes,
              const int numberOfEstimatedParameters,
              const Eigen::MatrixXd inverseOfAprioriCovariance = Eigen::MatrixXd::Zero( 0, 0 ),
              const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > initialParameterDeviationEstimate =
            Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( 0, 1 ) ):
        observationsAndTimes_( observationsAndTimes ), initialParameterDeviationEstimate_( initialParameterDeviationEstimate ),
        inverseOfAprioriCovariance_( inverseOfAprioriCovariance )
    {
        if( inverseOfAprioriCovariance_.rows( ) == 0 )
        {
            inverseOfAprioriCovariance_ = Eigen::MatrixXd::Zero( numberOfEstimatedParameters, numberOfEstimatedParameters );
        }

        if( ( numberOfEstimatedParameters != inverseOfAprioriCovariance_.rows( ) ) ||
                ( numberOfEstimatedParameters != inverseOfAprioriCovariance_.cols( ) ) )
        {
            throw std::runtime_error( "Error when making POD input, size of a priori covariance is inconsistent" );
        }

        if( initialParameterDeviationEstimate_.rows( ) == 0 )
        {
            initialParameterDeviationEstimate_ =
                    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >::Zero( numberOfEstimatedParameters, 1 );
        }

        if( numberOfEstimatedParameters != initialParameterDeviationEstimate_.rows( ) )
        {
            throw std::runtime_error( "Error when making POD input, size of initial parameter deviation is inconsistent" );
        }

        setConstantWeightsMatrix( 1.0 );
    }

    //! Function to set a constant values for all observation weights
    /*!
     * Function to set a constant values for all observation weights
     * \param constantWeight Constant weight that is to be set for all observations
     */
    void setConstantWeightsMatrix( const double constantWeight = 1.0 )
    {
        std::map< observation_models::ObservableType,
                    std::map< observation_models::LinkEnds, double > > weightPerObservableAndLinkEnds;
        for( typename PodInputDataType::const_iterator observablesIterator = observationsAndTimes_.begin( );
             observablesIterator != observationsAndTimes_.end( ); observablesIterator++ )
        {
            for( typename SingleObservablePodInputType::const_iterator dataIterator =
                 observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
            {
                weightPerObservableAndLinkEnds[ observablesIterator->first ][ dataIterator->first ] = constantWeight;
            }
        }
        setConstantPerObservableAndLinkEndsWeights( weightPerObservableAndLinkEnds );
    }

    //! Function to set a values for observation weights, constant per observable type
    /*!
     * Function to set a values for observation weights, constant per observable type
     * \param weightPerObservable Values for observation weights, constant per observable type
     */
    void setConstantPerObservableWeightsMatrix(
            const std::map< observation_models::ObservableType, double > weightPerObservable )
    {
        std::map< observation_models::ObservableType,
                std::map< observation_models::LinkEnds, double > > weightPerObservableAndLinkEnds;
        for( typename PodInputDataType::const_iterator observablesIterator = observationsAndTimes_.begin( );
             observablesIterator != observationsAndTimes_.end( ); observablesIterator++ )
        {
            if( weightPerObservable.count( observablesIterator->first ) != 0 )
            {
                for( typename SingleObservablePodInputType::const_iterator dataIterator =
                     observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
                {
                    weightPerObservableAndLinkEnds[ observablesIterator->first ][ dataIterator->first ] =
                            weightPerObservable.at( observablesIterator->first );
                }
            }
        }
        setConstantPerObservableAndLinkEndsWeights( weightPerObservableAndLinkEnds );
    }

    //! Function to set a values for observation weights, constant per observable type and link ends type
    /*!
     * Function to set a values for observation weights, constant per observable type and link ends type
     * \param weightPerObservableAndLinkEnds Values for observation weights, constant per observable type and link ends type
     */
    void setConstantPerObservableAndLinkEndsWeights(
            const std::map< observation_models::ObservableType,
            std::map< observation_models::LinkEnds, double > > weightPerObservableAndLinkEnds )
    {
        for( typename PodInputDataType::const_iterator observablesIterator = observationsAndTimes_.begin( );
             observablesIterator != observationsAndTimes_.end( ); observablesIterator++ )
        {
            if( weightPerObservableAndLinkEnds.count( observablesIterator->first ) == 0 )
            {
                std::string errorMessage =
                        "Error when setting  weights per observable, observable " +
                        boost::lexical_cast< std::string >( observablesIterator->first  ) + " not found";
                throw std::runtime_error( errorMessage );
            }
            else
            {
                for( typename SingleObservablePodInputType::const_iterator dataIterator =
                     observablesIterator->second.begin( ); dataIterator != observablesIterator->second.end( ); dataIterator++  )
                {
                    if( weightPerObservableAndLinkEnds.at( observablesIterator->first ).count( dataIterator->first ) == 0 )
                    {
                        std::string errorMessage =
                                "Error when setting  weights per observable, link ends not found for observable " +
                                boost::lexical_cast< std::string >( observablesIterator->first  );
                        throw std::runtime_error( errorMessage );
                    }
                    else
                    {

                        weightsMatrixDiagonals_[ observablesIterator->first ][ dataIterator->first ] =
                                Eigen::VectorXd::Constant(
                                    dataIterator->second.first.rows( ),
                                    weightPerObservableAndLinkEnds.at( observablesIterator->first ).at( dataIterator->first ) );
                    }
                }
            }
        }
    }

    //! Function to return the total data structure of observations and associated times/link ends/type (by reference)
    /*!
     * Function to return the total data structure of observations and associated times/link ends/type (by reference)
     * \return Total data structure of observations and associated times/link ends/type (by reference)
     */
    PodInputDataType& getObservationsAndTimes( )
    {
        return observationsAndTimes_;
    }

    //! Function to return the correction to estimated parameter vector to be applied on first iteration
    /*!
     * Function to return the correction to estimated parameter vector to be applied on first iteration
     * \return Correction to estimated parameter vector to be applied on first iteration
     */
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > getInitialParameterDeviationEstimate( )
    {
        return initialParameterDeviationEstimate_;
    }

    //! A priori covariance matrix (unnormalized) of estimated parameters
    //! Function to return the a priori covariance matrix (unnormalized) of estimated parameters
    /*!
     * Function to return the a priori covariance matrix (unnormalized) of estimated parameters
     * \return A priori covariance matrix (unnormalized) of estimated parameters
     */
    Eigen::MatrixXd getInverseOfAprioriCovariance( )
    {
        return inverseOfAprioriCovariance_;
    }

    //! Function to return the weight matrix diagonals, sorted by link ends and observable type (by reference)
    /*!
     * Function to return the weight matrix diagonals, sorted by link ends and observable type (by reference)
     * \return Weight matrix diagonals, sorted by link ends and observable type (by reference)
     */
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > >&
    getWeightsMatrixDiagonals( )
    {
        return weightsMatrixDiagonals_;
    }

private:
    //! Total data structure of observations and associated times/link ends/type
    PodInputDataType observationsAndTimes_;

    //! Correction to estimated parameter vector to be applied on first iteration
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > initialParameterDeviationEstimate_;

    //! A priori covariance matrix (unnormalized) of estimated parameters
    Eigen::MatrixXd inverseOfAprioriCovariance_;

    //! Weight matrix diagonals, sorted by link ends and observable type
    std::map< observation_models::ObservableType, std::map< observation_models::LinkEnds, Eigen::VectorXd > >
    weightsMatrixDiagonals_;


};

//! Data structure through which the output of the orbit determination is communicated
template< typename ObservationScalarType = double >
struct PodOutput
{

    //! Constructor
    /*!
     * Constructor
     * \param parameterEstimate Vector of estimated parameter values.
     * \param residuals Vector of postfit observation residuals
     * \param normalizedInformationMatrix Matrix of observation partials (normalixed) used in estimation
     * (may be empty if so requested)
     * \param weightsMatrixDiagonal Diagonal of weights matrix used in the estimation
     * \param informationMatrixTransformationDiagonal Vector of values by which the columns of the unnormalized information
     * matrix were divided to normalize its entries.
     * \param inverseNormalizedCovarianceMatrix Inverse of postfit normalized covariance matrix
     * \param residualStandardDeviation Standard deviation of postfit residuals vector
     */
    PodOutput( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& parameterEstimate,
               const Eigen::VectorXd& residuals,
               const Eigen::MatrixXd& normalizedInformationMatrix,
               const Eigen::VectorXd& weightsMatrixDiagonal,
               const Eigen::VectorXd& informationMatrixTransformationDiagonal,
               const Eigen::MatrixXd& inverseNormalizedCovarianceMatrix,
               const double residualStandardDeviation ):

        parameterEstimate_( parameterEstimate ), residuals_( residuals ),
        normalizedInformationMatrix_( normalizedInformationMatrix ), weightsMatrixDiagonal_( weightsMatrixDiagonal ),
        informationMatrixTransformationDiagonal_( informationMatrixTransformationDiagonal ),
        inverseNormalizedCovarianceMatrix_( inverseNormalizedCovarianceMatrix ),
        residualStandardDeviation_( residualStandardDeviation ){ }

    Eigen::MatrixXd getUnnormalizedInverseCovarianceMatrix( )
    {
        Eigen::MatrixXd inverseUnnormalizedCovarianceMatrix = inverseNormalizedCovarianceMatrix_;

        for( int i = 0; i < informationMatrixTransformationDiagonal_.rows( ); i++ )
        {
            for( int j = 0; j < informationMatrixTransformationDiagonal_.rows( ); j++ )
            {
                inverseUnnormalizedCovarianceMatrix( i, j ) *=
                        informationMatrixTransformationDiagonal_( i ) * informationMatrixTransformationDiagonal_( j );
            }

        }

        return inverseUnnormalizedCovarianceMatrix;
    }

    //! Vector of estimated parameter values.
    Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > parameterEstimate_;

    //! Vector of postfit observation residuals
    Eigen::VectorXd residuals_;

    //! Matrix of observation partials (normalixed) used in estimation (may be empty if so requested)
    Eigen::MatrixXd normalizedInformationMatrix_;

    //! Diagonal of weights matrix used in the estimation
    Eigen::VectorXd weightsMatrixDiagonal_;

    //! Vector of values by which the columns of the unnormalized information matrix were divided to normalize its entries.
    Eigen::VectorXd informationMatrixTransformationDiagonal_;

    //! Inverse of postfit normalized covariance matrix
    Eigen::MatrixXd inverseNormalizedCovarianceMatrix_;

    //! Standard deviation of postfit residuals vector
    double residualStandardDeviation_;
};

}

}

#endif // TUDAT_PODINPUTOUTPUTTYPES_H
