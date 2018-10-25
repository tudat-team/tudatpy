/*    Copyright (c) 2010-2018, Delft University of Technology
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
#include <iostream>
#include <memory>

#include <Eigen/Core>
#include <Eigen/LU>

#include "Tudat/Basics/timeType.h"
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
        inverseOfAprioriCovariance_( inverseOfAprioriCovariance ),
        reintegrateEquationsOnFirstIteration_( true ),
        reintegrateVariationalEquations_( true ),
        saveInformationMatrix_( true ),
        printOutput_( true ),
        saveResidualsAndParametersFromEachIteration_( true ),
        saveStateHistoryForEachIteration_( false )
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

    //! Destructor
    virtual ~PodInput( ){ }

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
                        std::to_string( observablesIterator->first  ) + " not found";
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
                                std::to_string( observablesIterator->first  );
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

    //! Function to define specific settings for estimation process
    /*!
     *  Function to define specific settings for estimation process
     *  \param reintegrateEquationsOnFirstIteration Boolean denoting whether the dynamics and variational equations are to
     *  be reintegrated on first iteration, or if existing values are to be used to perform first iteration.
     *  \param reintegrateVariationalEquations Boolean denoting whether the variational equations are to be reintegrated during
     *  estimation
     *  \param saveInformationMatrix Boolean denoting whether to save the partials matrix in the output
     *  \param printOutput Boolean denoting whether to print output to th terminal when running the estimation.
     *  \param saveResidualsAndParametersFromEachIteration Boolean denoting whether the residuals and parameters from the each
     *  iteration are to be saved
     *  \param saveStateHistoryForEachIteration Boolean denoting whether the state history is to be saved on each iteration
     */
    void defineEstimationSettings( const bool reintegrateEquationsOnFirstIteration = 1,
                                   const bool reintegrateVariationalEquations = 1,
                                   const bool saveInformationMatrix = 1,
                                   const bool printOutput = 1,
                                   const bool saveResidualsAndParametersFromEachIteration = 1,
                                   const bool saveStateHistoryForEachIteration = 0 )
    {
        reintegrateEquationsOnFirstIteration_ = reintegrateEquationsOnFirstIteration;
        reintegrateVariationalEquations_ = reintegrateVariationalEquations;
        saveInformationMatrix_ = saveInformationMatrix;
        printOutput_ = printOutput;
        saveResidualsAndParametersFromEachIteration_ = saveResidualsAndParametersFromEachIteration;
        saveStateHistoryForEachIteration_ = saveStateHistoryForEachIteration;
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

    //! Function to return the boolean denoting whether the dynamics and variational equations are reintegrated on first iteration
    /*!
     * Function to return the boolean denoting whether the dynamics and variational equations are to be reintegrated on first
     * iteration
     * \return Boolean denoting whether the dynamics and variational equations are to be reintegrated on first iteration
     */
    bool getReintegrateEquationsOnFirstIteration( )
    {
        return reintegrateEquationsOnFirstIteration_;
    }

    //! Function to return the boolean denoting whether the variational equations are to be reintegrated during estimation
    /*!
     * Function to return the boolean denoting whether the variational equations are to be reintegrated during estimation
     * \return Boolean denoting whether the variational equations are to be reintegrated during estimation
     */
    bool getReintegrateVariationalEquations( )
    {
        return reintegrateVariationalEquations_;
    }

    //! Function to return the boolean denoting whether to print output to th terminal when running the estimation.
    /*!
     * Function to return the boolean denoting whether to print output to th terminal when running the estimation.
     * \return Boolean denoting whether to print output to th terminal when running the estimation.
     */
    bool getSaveInformationMatrix( )
    {
        return saveInformationMatrix_;
    }

    //! Function to return the boolean denoting whether to print output to th terminal when running the estimation.
    /*!
     * Function to return the boolean denoting whether to print output to th terminal when running the estimation.
     * \return Boolean denoting whether to print output to th terminal when running the estimation.
     */
    bool getPrintOutput( )
    {
        return printOutput_;
    }

    //! Function to return the boolean denoting whether the residuals and parameters from the each iteration are to be saved
    /*!
     * Function to return the boolean denoting whether the residuals and parameters from the each iteration are to be saved
     * \return Boolean denoting whether the residuals and parameters from the each iteration are to be saved
     */
    bool getSaveResidualsAndParametersFromEachIteration( )
    {
        return saveResidualsAndParametersFromEachIteration_;
    }

    //! Function to return the boolean denoting whether the state history is to be saved on each iteration.
    /*!
     * Function to return the boolean denoting whether the state history is to be saved on each iteration.
     * \return Boolean denoting whether the state history is to be saved on each iteration.
     */
    bool getSaveStateHistoryForEachIteration( )
    {
        return saveStateHistoryForEachIteration_;
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

    //!  Boolean denoting whether the dynamics and variational equations are to be reintegrated on first iteration
    bool reintegrateEquationsOnFirstIteration_;

    //! Boolean denoting whether the variational equations are to be reintegrated during estimation
    bool reintegrateVariationalEquations_;

    //! Boolean denoting whether to print output to th terminal when running the estimation.
    bool saveInformationMatrix_;

    //! Boolean denoting whether to print output to th terminal when running the estimation.
    bool printOutput_;

    //! Boolean denoting whether the residuals and parameters from the each iteration are to be saved
    bool saveResidualsAndParametersFromEachIteration_;

    //! Boolean denoting whether the state history is to be saved on each iteration.
    bool saveStateHistoryForEachIteration_;

};

//! Class that is used during the orbit determination/parameter estimation to determine whether the estimation is converged.
class EstimationConvergenceChecker
{
public:

    //! Constructor
    /*!
     * Constructor, sets a number of values for stopping conditions. The estimation stops if one of these is met.
     * \param maximumNumberOfIterations Maximum number of allowed iterations for estimation
     * \param minimumResidualChange Minimum required change in residual between two iterations
     * \param minimumResidual Minimum value of observation residual below which estimation is converged
     * \param numberOfIterationsWithoutImprovement Number of iterations without reduction of residual
     */
    EstimationConvergenceChecker(
            const unsigned int maximumNumberOfIterations = 5,
            const double minimumResidualChange = 0.0,
            const double minimumResidual = 1.0E-20,
            const int numberOfIterationsWithoutImprovement = 2 ):
        maximumNumberOfIterations_( maximumNumberOfIterations ), minimumResidualChange_( minimumResidualChange ),
        minimumResidual_( minimumResidual ),
        numberOfIterationsWithoutImprovement_( numberOfIterationsWithoutImprovement )
    { }

    //! Function to determine whether the estimation is deemed to be converged
    /*!
     * Function to determine whether the estimation is deemed to be converged (i.e. if it should terminate)
     * \param numberOfIterations Number of iterations of estimation procedure that have been completed
     * \param rmsResidualHistory Rms residuals at current and all previous iterations
     * \return True if estimation is to be terminated
     */
    bool isEstimationConverged( const int numberOfIterations, const std::vector< double > rmsResidualHistory )
    {
        bool isConverged = 0;
        if( numberOfIterations >= maximumNumberOfIterations_ )
        {
            std::cout << "Maximum number of iterations reached" << std::endl;
            isConverged = 1;
        }
        if( rmsResidualHistory[ rmsResidualHistory.size( ) - 1 ] < minimumResidual_ )
        {
            std::cout << "Required residual level achieved" << std::endl;
            isConverged = 1;
        }
        if( ( std::distance( rmsResidualHistory.begin( ), std::max_element(
                                 rmsResidualHistory.begin( ), rmsResidualHistory.end( ) ) ) - rmsResidualHistory.size( ) ) <
                numberOfIterationsWithoutImprovement_ )
        {
            std::cout << "Too many iterations without parameter improvement" << std::endl;
            isConverged = 1;
        }
        if( rmsResidualHistory.size( ) > 1 )
        {
            if( std::fabs( rmsResidualHistory.at( rmsResidualHistory.size( )  - 1 ) -
                           rmsResidualHistory.at( rmsResidualHistory.size( )  - 2 ) ) < minimumResidualChange_ )
            {
                isConverged = 1;
            }
        }
        return isConverged;
    }
protected:

    //! Maximum number of allowed iterations for estimation
    int maximumNumberOfIterations_;

    //! Minimum required change in residual between two iterations
    double minimumResidualChange_;

    //! Minimum value of observation residual below which estimation is converged
    double minimumResidual_;

    //!  Number of iterations without reduction of residual
    unsigned int numberOfIterationsWithoutImprovement_;
};


//! Data structure through which the output of the orbit determination is communicated
template< typename ObservationScalarType = double, typename TimeType = double  >
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
     * \param residualHistory Vector of residuals per iteration
     * \param parameterHistory Vector of parameter vectors per iteration (entry 1 is pre-estimation values)
     * \param exceptionDuringInversion Boolean denoting whether an exception was caught during inversion of normal equations
     * \param exceptionDuringPropagation Boolean denoting whether an exception was caught during (re)propagation of equations of
     * motion (and variational equations).
     */
    PodOutput( const Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 >& parameterEstimate,
               const Eigen::VectorXd& residuals,
               const Eigen::MatrixXd& normalizedInformationMatrix,
               const Eigen::VectorXd& weightsMatrixDiagonal,
               const Eigen::VectorXd& informationMatrixTransformationDiagonal,
               const Eigen::MatrixXd& inverseNormalizedCovarianceMatrix,
               const double residualStandardDeviation,
               const std::vector< Eigen::VectorXd >& residualHistory = std::vector< Eigen::VectorXd >( ),
               const std::vector< Eigen::VectorXd >& parameterHistory = std::vector< Eigen::VectorXd >( ),
               const bool exceptionDuringInversion = false,
               const bool exceptionDuringPropagation = false ):
        parameterEstimate_( parameterEstimate ), residuals_( residuals ),
        normalizedInformationMatrix_( normalizedInformationMatrix ), weightsMatrixDiagonal_( weightsMatrixDiagonal ),
        informationMatrixTransformationDiagonal_( informationMatrixTransformationDiagonal ),
        inverseNormalizedCovarianceMatrix_( inverseNormalizedCovarianceMatrix ),
        residualStandardDeviation_( residualStandardDeviation ),
        residualHistory_( residualHistory ),
        parameterHistory_( parameterHistory ),
        exceptionDuringInversion_( exceptionDuringInversion ),
        exceptionDuringPropagation_( exceptionDuringPropagation)
    { }

    //! Function to retrieve the unnormalized inverse estimation covariance matrix
    /*!
     * Function to retrieve the unnormalized inverse estimation covariance matrix
     * \return Isnverse estimation covariance matrix
     */
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

    //! Function to retrieve the unnormalized estimation covariance matrix
    /*!
     * Function to retrieve the unnormalized estimation covariance matrix
     * \return Estimation covariance matrix
     */
    Eigen::MatrixXd getUnnormalizedCovarianceMatrix( )
    {
        return getUnnormalizedInverseCovarianceMatrix( ).inverse( );
    }

    //! Function to retrieve the unnormalized formal error vector of the estimation result.
    /*!
     * Function to retrieve the unnormalized formal error vector of the estimation result.
     * \return Formal error vector of the estimation result.
     */
    Eigen::VectorXd getFormalErrorVector( )
    {
        return ( getUnnormalizedCovarianceMatrix( ).diagonal( ) ).cwiseSqrt( );
    }

    //! Function to retrieve the correlation matrix of the estimation result.
    /*!
     * Function to retrieve the correlation matrix of the estimation result.
     * \return Correlation matrix of the estimation result.
     */
    Eigen::MatrixXd getCorrelationMatrix( )
    {
        return getUnnormalizedCovarianceMatrix( ).cwiseQuotient( getFormalErrorVector( ) * getFormalErrorVector( ).transpose( ) );
    }

    //! Function to get residual vectors per iteration concatenated into a matrix
    /*!
     * Function to get residual vectors per iteration concatenated into a matrix (one column per iteration).
     * \return Residual vectors per iteration concatenated into a matrix
     */
    Eigen::MatrixXd getResidualHistoryMatrix( )
    {
        if( residualHistory_.size( ) > 0 )
        {
            Eigen::MatrixXd residualHistoryMatrix = Eigen::MatrixXd( residualHistory_.at( 0 ).rows( ), residualHistory_.size( ) );
            for( unsigned int i = 0; i < residualHistory_.size( ); i++ )
            {
                residualHistoryMatrix.block( 0, i, residualHistory_.at( 0 ).rows( ), 1 ) = residualHistory_.at( i );
            }
            return residualHistoryMatrix;
        }
        else
        {
            std::cerr << "Warning, requesting residual history, but history not saved." << std::endl;
            return Eigen::MatrixXd::Zero( 0, 0 );
        }
    }

    //! Function to get parameter vectors per iteration concatenated into a matrix
    /*!
     * Function to get parameter vectors per iteration concatenated into a matrix (one column per iteration). Column 0 contains
     * pre-estimation values
     * \return Parameter vectors per iteration concatenated into a matrix
     */
    Eigen::MatrixXd getParameterHistoryMatrix( )
    {
        if( parameterHistory_.size( ) > 0 )
        {
            Eigen::MatrixXd parameterHistoryMatrix = Eigen::MatrixXd( parameterHistory_.at( 0 ).rows( ), parameterHistory_.size( ) );
            for( unsigned int i = 0; i < parameterHistory_.size( ); i++ )
            {
                parameterHistoryMatrix.block( 0, i, parameterHistory_.at( 0 ).rows( ), 1 ) = parameterHistory_.at( i );
            }
            return parameterHistoryMatrix;
        }
        else
        {
            std::cerr << "Warning, requesting parameter history, but history not saved." << std::endl;
            return Eigen::MatrixXd::Zero( 0, 0 );
        }
    }

    //! Function to set the full state histories of numerical solutions and dependent variables
    /*!
     * Function to set the full state histories of numerical solutions and dependent variables
     * \param dynamicsHistoryPerIteration List of numerical solutions of dynamics (per iteration, per arc)
     * \param dependentVariableHistoryPerIteration List of numerical solutions of dependent variables (per iteration, per arc)
     */
    void setStateHistories(
            std::vector< std::vector< std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > >
            dynamicsHistoryPerIteration,
            std::vector< std::vector< std::map< TimeType, Eigen::VectorXd > > > dependentVariableHistoryPerIteration )
    {
        dynamicsHistoryPerIteration_ = dynamicsHistoryPerIteration;
        dependentVariableHistoryPerIteration_ = dependentVariableHistoryPerIteration;
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

    //! Vector of residuals per iteration
    std::vector< Eigen::VectorXd > residualHistory_;

    //! Vector of parameter vectors per iteration (entry 0 is pre-estimation values)
    std::vector< Eigen::VectorXd > parameterHistory_;

    //! List of numerical solutions of dynamics (per iteration, per arc)
    std::vector< std::vector< std::map< TimeType, Eigen::Matrix< ObservationScalarType, Eigen::Dynamic, 1 > > > > dynamicsHistoryPerIteration_;

    //! List of numerical solutions of dependent variables (per iteration, per arc)
    std::vector< std::vector< std::map< TimeType, Eigen::VectorXd > > > dependentVariableHistoryPerIteration_;

    //! Boolean denoting whether an exception was caught during inversion of normal equations
    bool exceptionDuringInversion_;

    //! Boolean denoting whether an exception was caught during (re)propagation of equations of motion (and variational equations)
    bool exceptionDuringPropagation_;
};


extern template class PodInput< double, double >;
extern template struct PodOutput< double >;

#if( BUILD_EXTENDED_PRECISION_PROPAGATION_TOOLS )
extern template class PodInput< long double, double >;
extern template class PodInput< double, Time >;
extern template class PodInput< long double, Time >;

extern template struct PodOutput< long double >;
#endif

}

}

#endif // TUDAT_PODINPUTOUTPUTTYPES_H
