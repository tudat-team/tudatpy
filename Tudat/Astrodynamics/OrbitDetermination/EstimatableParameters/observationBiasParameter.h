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

//! Interface class for the estimation of a constant absolute observation bias.
/*!
 *  Interface class for the estimation of a constant absolute observation bias (at given link ends and observable type).
 *  Unlike most other EstimatableParameter derived
 *  classes, this class does not have direct access to the class (ConstantObservationBias) used in the simulations for the
 *  observation bias. This is due to the fact that the ConstantObservationBias class is templated by the observable size, whil
 *  this class is not.
 */
class ConstantObservationBiasParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param getCurrentBias Function to retrieve the current observation bias.
     * \param resetCurrentBias Function to reset the current observation bias
     * \param linkEnds Observation link ends for which the bias is active.
     * \param observableType Observable type for which the bias is active.
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

    //! Function to get the current value of the constant absolute observation bias that is to be estimated.
    /*!
     * Function to get the current value of the constant absolute observation bias that is to be estimated.
     * \return Current value of the constant absolute observation bias that is to be estimated.
     */
    Eigen::VectorXd getParameterValue( )
    {
        return getCurrentBias_( );
    }

    //! Function to reset the value of the constant absolute observation bias that is to be estimated.
    /*!
     * Function to reset the value of the constant absolute observation bias that is to be estimated.
     * \param parameterValue New value of the constant absolute observation bias that is to be estimated.
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        resetCurrentBias_( parameterValue );
    }

    //! Function to retrieve the size of the parameter (equal to the size of the observable).
    /*!
     *  Function to retrieve the size of the parameter (equal to the size of the observable).
     *  \return Size of parameter value (equal to the size of the observable).
     */
    int getParameterSize( )
    {
        return observation_models::getObservableSize( observableType_ );
    }

    //! Function to reset the get/set function of the observation bias
    /*!
     * Function to reset the get/set function of the observation bias. This function is needed since te observation models/biases
     * are typically created after the estimated parameter objects
     * \param getCurrentBias New function to retrieve the current observation bias.
     * \param resetCurrentBias New function to reset the current observation bias
     */
    void setObservationBiasFunctions(
            const boost::function< Eigen::VectorXd( ) > getCurrentBias,
            const boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias )
    {
        // Check if functions already exist
        if( !getCurrentBias_.empty( ) || !resetCurrentBias_.empty( ) )
        {
            std::cerr << "Warning when resetting absolute observation bias in estimation object, existing contents not empty" << std::endl;
        }

        getCurrentBias_ = getCurrentBias;
        resetCurrentBias_ = resetCurrentBias;
    }

    //! Function to retrieve the observation link ends for which the bias is active.
    /*!
     * Function to retrieve the observation link ends for which the bias is active.
     * \return Observation link ends for which the bias is active.
     */
    observation_models::LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }

    //! Function to retrieve the observable type for which the bias is active.
    /*!
     * Function to retrieve the observable type ends for which the bias is active.
     * \return Observable type for which the bias is active.
     */
    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "for observable: (" +
                observation_models::getObservableName( observableType_, linkEnds_.size( )  ) + ") and link ends: (" +
                observation_models::getLinkEndsString( linkEnds_ ) + ")";
        return parameterDescription;
    }

protected:

private:

    //! Function to retrieve the current observation bias.
    boost::function< Eigen::VectorXd( ) > getCurrentBias_;

    //! Function to reset the current observation bia
    boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias_;

    //! Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

};


class ArcWiseObservationBiasParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param getCurrentBias Function to retrieve the current observation bias.
     * \param resetCurrentBias Function to reset the current observation bias
     * \param linkEnds Observation link ends for which the bias is active.
     * \param observableType Observable type for which the bias is active.
     */
    ArcWiseObservationBiasParameter(
            const std::vector< double > arcStartTimes,
            const boost::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
            const boost::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList,
            const int linkEndIndex,
            const observation_models::LinkEnds linkEnds,
            const observation_models::ObservableType observableType ):
        EstimatableParameter< Eigen::VectorXd >( arcwise_constant_additive_observation_bias, linkEnds.begin( )->second. first ),
        arcStartTimes_( arcStartTimes ), getBiasList_( getBiasList ), resetBiasList_( resetBiasList ),
        linkEndIndex_( linkEndIndex ), linkEnds_( linkEnds ), observableType_( observableType )
    {
        observableSize_ = observation_models::getObservableSize( observableType );
        numberOfArcs_ = arcStartTimes.size( );
    }

    //! Destructor
    ~ArcWiseObservationBiasParameter( ) { }

    //! Function to get the current value of the arc-wise absolute observation bias that is to be estimated.
    /*!
     * Function to get the current value of the arc-wise absolute observation bias that is to be estimated.
     * \return Current value of the arc-wise absolute observation bias that is to be estimated.
     */
    Eigen::VectorXd getParameterValue( )
    {
        std::vector< Eigen::VectorXd > observationBiases = getBiasList_( );
        Eigen::VectorXd currentParameterSet = Eigen::VectorXd::Zero(
                    observableSize_ * observationBiases.size( ) );
        for( unsigned int i = 0; i < observationBiases.size( ); i++ )
        {
            currentParameterSet.segment( i * observableSize_, observableSize_ ) = observationBiases.at( i );
        }
        return currentParameterSet;
    }

    //! Function to reset the value of the constant absolute observation bias that is to be estimated.
    /*!
     * Function to reset the value of the constant absolute observation bias that is to be estimated.
     * \param parameterValue New value of the constant absolute observation bias that is to be estimated.
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        std::vector< Eigen::VectorXd > observationBiases;

        for( int i = 0; i < numberOfArcs_; i++ )
        {
           observationBiases.push_back( parameterValue.segment( i * observableSize_, observableSize_ ) );
        }

        resetBiasList_( observationBiases );
    }

    //! Function to retrieve the size of the parameter (equal to the size of the observable).
    /*!
     *  Function to retrieve the size of the parameter (equal to the size of the observable).
     *  \return Size of parameter value (equal to the size of the observable).
     */
    int getParameterSize( )
    {
        return observableSize_ * numberOfArcs_;
    }

    //! Function to reset the get/set function of the observation bias
    /*!
     * Function to reset the get/set function of the observation bias. This function is needed since te observation models/biases
     * are typically created after the estimated parameter objects
     * \param getCurrentBias New function to retrieve the current observation bias.
     * \param resetCurrentBias New function to reset the current observation bias
     */
    void setObservationBiasFunctions(
            const boost::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
            const boost::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList )
    {
        // Check if functions already exist
        if( !getBiasList_.empty( ) || !resetBiasList_.empty( ) )
        {
            std::cerr << "Warning when resetting arc-wise absolute observation bias in estimation object, existing contents not empty" << std::endl;
        }

        getBiasList_ = getBiasList;
        resetBiasList_ = resetBiasList;
    }

    //! Function to retrieve the observation link ends for which the bias is active.
    /*!
     * Function to retrieve the observation link ends for which the bias is active.
     * \return Observation link ends for which the bias is active.
     */
    observation_models::LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }

    //! Function to retrieve the observable type for which the bias is active.
    /*!
     * Function to retrieve the observable type ends for which the bias is active.
     * \return Observable type for which the bias is active.
     */
    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "for observable: (" +
                observation_models::getObservableName( observableType_, linkEnds_.size( )  ) + ") and link ends: (" +
                observation_models::getLinkEndsString( linkEnds_ ) + ")";
        return parameterDescription;
    }

    const std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }


    int getLinkEndIndex( )
    {
        return linkEndIndex_;
    }

    boost::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( )
    {
        return lookupScheme_;
    }

    void setLookupScheme( const boost::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme )
    {
        lookupScheme_ = lookupScheme;
    }

protected:

private:

    const std::vector< double > arcStartTimes_;

    //! Function to retrieve the current observation bias.
    boost::function< std::vector< Eigen::VectorXd >( ) > getBiasList_;

    //! Function to reset the current observation bia
    boost::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList_;

    int linkEndIndex_;

    //! Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    int observableSize_;

    int numberOfArcs_;

    boost::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;
};

//! Interface class for the estimation of a constant relative observation bias.
/*!
 *  Interface class for the estimation of a constant relative observation bias (at given link ends and observable type).
 *  Unlike most other EstimatableParameter derived
 *  classes, this class does not have direct access to the class (ConstantRelativeObservationBias) used in the simulations for the
 *  observation bias. This is due to the fact that the ConstantRelativeObservationBias class is templated by the observable size,
 *  while this class is not.
 */
class ConstantRelativeObservationBiasParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor

    //! Constructor
    /*!
     * Constructor
     * \param getCurrentBias Function to retrieve the current relative observation bias.
     * \param resetCurrentBias Function to reset the current relative observation bias
     * \param linkEnds Observation link ends for which the bias is active.
     * \param observableType Observable type for which the bias is active.
     */
    ConstantRelativeObservationBiasParameter(
            const boost::function< Eigen::VectorXd( ) > getCurrentBias,
            const boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias,
            const observation_models::LinkEnds linkEnds,
            const observation_models::ObservableType observableType ):
        EstimatableParameter< Eigen::VectorXd >( constant_relative_observation_bias, linkEnds.begin( )->second.first ),
        getCurrentBias_( getCurrentBias ), resetCurrentBias_( resetCurrentBias ),
        linkEnds_( linkEnds ), observableType_( observableType ){ }

    //! Destructor
    ~ConstantRelativeObservationBiasParameter( ) { }

    //! Function to get the current value of the constant relative observation bias that is to be estimated.
    /*!
     * Function to get the current value of the constant relative observation bias that is to be estimated.
     * \return Current value of the constant relative observation bias that is to be estimated.
     */
    Eigen::VectorXd getParameterValue( )
    {
        return getCurrentBias_( );
    }

    //! Function to reset the value of the constant relative observation bias that is to be estimated.
    /*!
     * Function to reset the value of the constant relative observation bias that is to be estimated.
     * \param parameterValue New value of the constant relative observation bias that is to be estimated.
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        resetCurrentBias_( parameterValue );
    }

    //! Function to retrieve the size of the parameter (equal to the size of the observable).
    /*!
     *  Function to retrieve the size of the parameter (equal to the size of the observable).
     *  \return Size of parameter value (equal to the size of the observable).
     */
    int getParameterSize( )
    {
        return observation_models::getObservableSize( observableType_ );
    }

    //! Function to reset the get/set function of the observation bias
    /*!
     * Function to reset the get/set function of the observation bias. This function is needed since te observation models/biases
     * are typically created after the estimated parameter objects
     * \param getCurrentBias New function to retrieve the current observation bias.
     * \param resetCurrentBias New function to reset the current observation bias
     */
    void setObservationBiasFunctions(
            const boost::function< Eigen::VectorXd( ) > getCurrentBias,
            const boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias )
    {
        // Check if functions already exist
        if( !getCurrentBias_.empty( ) || !resetCurrentBias_.empty( ) )
        {
            std::cerr << "Warning when resetting relative observation bias in estimation object, existing contents not empty" << std::endl;
        }

        getCurrentBias_ = getCurrentBias;
        resetCurrentBias_ = resetCurrentBias;
    }

    //! Function to retrieve the observation link ends for which the bias is active.
    /*!
     * Function to retrieve the observation link ends for which the bias is active.
     * \return Observation link ends for which the bias is active.
     */
    observation_models::LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }

    //! Function to retrieve the observable type for which the bias is active.
    /*!
     * Function to retrieve the observable type ends for which the bias is active.
     * \return Observable type for which the bias is active.
     */
    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "for observable: (" +
                observation_models::getObservableName( observableType_, linkEnds_.size( ) ) + ") and link ends: (" +
                observation_models::getLinkEndsString( linkEnds_ ) + ")";
        return parameterDescription;
    }


protected:

private:

    //! Function to retrieve the current observation bias.
    boost::function< Eigen::VectorXd( ) > getCurrentBias_;

    //! Function to reset the current observation bia
    boost::function< void( const Eigen::VectorXd& ) > resetCurrentBias_;

    //! Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;
};

class ArcWiseRelativeObservationBiasParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param getCurrentBias Function to retrieve the current observation bias.
     * \param resetCurrentBias Function to reset the current observation bias
     * \param linkEnds Observation link ends for which the bias is active.
     * \param observableType Observable type for which the bias is active.
     */
    ArcWiseRelativeObservationBiasParameter(
            const std::vector< double > arcStartTimes,
            const boost::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
            const boost::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList,
            const int linkEndIndex,
            const observation_models::LinkEnds linkEnds,
            const observation_models::ObservableType observableType ):
        EstimatableParameter< Eigen::VectorXd >( arcwise_constant_relative_observation_bias, linkEnds.begin( )->second. first ),
        arcStartTimes_( arcStartTimes ), getBiasList_( getBiasList ), resetBiasList_( resetBiasList ),
        linkEndIndex_( linkEndIndex ), linkEnds_( linkEnds ), observableType_( observableType )
    {
        observableSize_ = observation_models::getObservableSize( observableType );
        numberOfArcs_ = arcStartTimes.size( );
    }

    //! Destructor
    ~ArcWiseRelativeObservationBiasParameter( ) { }

    //! Function to get the current value of the arc-wise absolute observation bias that is to be estimated.
    /*!
     * Function to get the current value of the arc-wise absolute observation bias that is to be estimated.
     * \return Current value of the arc-wise absolute observation bias that is to be estimated.
     */
    Eigen::VectorXd getParameterValue( )
    {
        std::vector< Eigen::VectorXd > observationBiases = getBiasList_( );
        Eigen::VectorXd currentParameterSet = Eigen::VectorXd::Zero(
                    observableSize_ * observationBiases.size( ) );
        for( unsigned int i = 0; i < observationBiases.size( ); i++ )
        {
            currentParameterSet.segment( i * observableSize_, observableSize_ ) = observationBiases.at( i );
        }
        return currentParameterSet;
    }

    //! Function to reset the value of the constant absolute observation bias that is to be estimated.
    /*!
     * Function to reset the value of the constant absolute observation bias that is to be estimated.
     * \param parameterValue New value of the constant absolute observation bias that is to be estimated.
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        std::vector< Eigen::VectorXd > observationBiases;

        for( int i = 0; i < numberOfArcs_; i++ )
        {
           observationBiases.push_back( parameterValue.segment( i * observableSize_, observableSize_ ) );
        }

        resetBiasList_( observationBiases );
    }

    //! Function to retrieve the size of the parameter (equal to the size of the observable).
    /*!
     *  Function to retrieve the size of the parameter (equal to the size of the observable).
     *  \return Size of parameter value (equal to the size of the observable).
     */
    int getParameterSize( )
    {
        return observableSize_ * numberOfArcs_;
    }

    //! Function to reset the get/set function of the observation bias
    /*!
     * Function to reset the get/set function of the observation bias. This function is needed since te observation models/biases
     * are typically created after the estimated parameter objects
     * \param getCurrentBias New function to retrieve the current observation bias.
     * \param resetCurrentBias New function to reset the current observation bias
     */
    void setObservationBiasFunctions(
            const boost::function< std::vector< Eigen::VectorXd >( ) > getBiasList,
            const boost::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList )
    {
        // Check if functions already exist
        if( !getBiasList_.empty( ) || !resetBiasList_.empty( ) )
        {
            std::cerr << "Warning when resetting arc-wise absolute observation bias in estimation object, existing contents not empty" << std::endl;
        }

        getBiasList_ = getBiasList;
        resetBiasList_ = resetBiasList;
    }

    //! Function to retrieve the observation link ends for which the bias is active.
    /*!
     * Function to retrieve the observation link ends for which the bias is active.
     * \return Observation link ends for which the bias is active.
     */
    observation_models::LinkEnds getLinkEnds( )
    {
        return linkEnds_;
    }

    //! Function to retrieve the observable type for which the bias is active.
    /*!
     * Function to retrieve the observable type ends for which the bias is active.
     * \return Observable type for which the bias is active.
     */
    observation_models::ObservableType getObservableType( )
    {
        return observableType_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription = getParameterTypeString( parameterName_.first ) + "for observable: (" +
                observation_models::getObservableName( observableType_, linkEnds_.size( )  ) + ") and link ends: (" +
                observation_models::getLinkEndsString( linkEnds_ ) + ")";
        return parameterDescription;
    }

    const std::vector< double > getArcStartTimes( )
    {
        return arcStartTimes_;
    }


    int getLinkEndIndex( )
    {
        return linkEndIndex_;
    }

    boost::shared_ptr< interpolators::LookUpScheme< double > > getLookupScheme( )
    {
        return lookupScheme_;
    }

    void setLookupScheme( const boost::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme )
    {
        lookupScheme_ = lookupScheme;
    }

protected:

private:

    const std::vector< double > arcStartTimes_;

    //! Function to retrieve the current observation bias.
    boost::function< std::vector< Eigen::VectorXd >( ) > getBiasList_;

    //! Function to reset the current observation bia
    boost::function< void( const std::vector< Eigen::VectorXd >& ) > resetBiasList_;

    int linkEndIndex_;

    //! Observation link ends for which the bias is active.
    observation_models::LinkEnds linkEnds_;

    //! Observable type for which the bias is active.
    observation_models::ObservableType observableType_;

    int observableSize_;

    int numberOfArcs_;

    boost::shared_ptr< interpolators::LookUpScheme< double > > lookupScheme_;
};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_OBSERVATIONBIASPARAMETER_H
