/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_DEPENDENTVARIABLESINTERFACE_H
#define TUDAT_DEPENDENTVARIABLESINTERFACE_H

#include <iostream>
#include <vector>

#include <memory>
#include <Eigen/Core>

#include "tudat/math/interpolators/oneDimensionalInterpolator.h"
#include "tudat/simulation/propagation_setup/propagationOutputSettings.h"
#include "tudat/simulation/propagation_setup/propagationOutput.h"

namespace tudat
{

namespace propagators
{

//! Base class for interface object of interpolation of numerically propagated dependent variables.
/*!
*  Base class for interface object of interpolation of numerically propagated dependent variables.
*  Derived classes implement the case of single-arc/multi-arc/hybrid combined dynamics.
*/
template< typename TimeType = double >
class DependentVariablesInterface
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dependentVariablesSettings Vector of single dependent variable settings
     */
    DependentVariablesInterface( )
    {

    }

    //! Destructor.
    virtual ~DependentVariablesInterface( ){ }

    //! Function to get the dependent variables at a given time.
    /*!
     *  Function to get the dependent variabless at a given time.
     *  \param evaluationTime Time at which to evaluate dependent variables interpolator
     *  \return Concatenated dependent variables.
     */
    virtual Eigen::VectorXd getDependentVariables( const TimeType evaluationTime ) = 0;

    virtual Eigen::VectorXd getSingleDependentVariable(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const TimeType evaluationTime ) = 0;
//
//    //! Function to get the value of a single dependent variable at a given time, from the dependent variable ID.
//    virtual Eigen::VectorXd getSingleDependentVariable(
//        const std::string dependentVariableId,
//        const int dependentVariableSize,
//        const TimeType evaluationTime ) = 0;

protected:


};

//! Interface object of interpolation of numerically propagated dependent variables for single-arc propagation/estimation.
template< typename TimeType = double >
class SingleArcDependentVariablesInterface : public DependentVariablesInterface< TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param dependentVariablesInterpolator Interpolator returning the dependent variables values as a function of time.
     * \param dependentVariablesSettings Vector of single dependent variables settings
     */
    SingleArcDependentVariablesInterface(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > dependentVariablesInterpolator,
            const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesSettings,
            const std::map< std::pair< int, int >, std::string > dependentVariableIds,
            const std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > > orderedDependentVariableSettings ):
            dependentVariablesSettings_( dependentVariablesSettings ),
            dependentVariablesInterpolator_( dependentVariablesInterpolator ),
            dependentVariableIds_( dependentVariableIds ),
            orderedDependentVariableSettings_( orderedDependentVariableSettings )
    {
        dependentVariablesSize_ = 0;
        dependentVariablesIdsAndIndices_.clear( );

        for ( unsigned int i= 0 ; i < dependentVariablesSettings_.size( ) ; i++ )
        {
            dependentVariablesTypes_.push_back( dependentVariablesSettings_[ i ]->dependentVariableType_ );
            dependentVariablesIdsAndIndices_[ getDependentVariableId( dependentVariablesSettings_[ i ] ) ] = dependentVariablesSize_;
            dependentVariablesSize_ += getDependentVariableSaveSize( dependentVariablesSettings_[ i ] );
        }
        dependentVariables_ = Eigen::VectorXd::Zero( dependentVariablesSize_ );

        if( dependentVariablesInterpolator_ != nullptr )
        {
            dependentVariablesInterpolator_->resetBoundaryHandling( interpolators::throw_exception_at_boundary );
        }
    }

    //! Destructor.
    ~SingleArcDependentVariablesInterface( ){ }

    //! Function to reset the dependent variables interpolator
    /*!
     * Function to reset the dependent variables interpolator
     * \param dependentVariablesInterpolator New interpolator returning the dependent variables values as a function of time.
     */
    void updateDependentVariablesInterpolator(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > dependentVariablesInterpolator )
    {
        dependentVariablesInterpolator_ = dependentVariablesInterpolator;
        if( dependentVariablesInterpolator_ != nullptr )
        {
            dependentVariablesInterpolator_->resetBoundaryHandling( interpolators::throw_exception_at_boundary );
        }
    }

    //! Function to get the interpolator returning the dependent variables as a function of time.
    /*!
     * Function to get the interpolator returning the dependent variables as a function of time.
     * \return Interpolator returning the dependent variable as a function of time.
     */
    std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > getDependentVariablesInterpolator( )
    {
        return dependentVariablesInterpolator_;
    }

    std::map< std::pair< int, int >, std::string > getDependentVariableIds( ) const
    {
        return dependentVariableIds_;
    }

    std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > > getOrderedDependentVariableSettings( ) const
    {
        return orderedDependentVariableSettings_;
    }

    //! Function to get the concatenated dependent variables values at a given time.
    /*!
     *  Function to get the concatenated dependent variables values at a given time.
     *  \param evaluationTime Time at which to evaluate dependent variables interpolator
     *  \return Dependent variables values
     */
    Eigen::VectorXd getDependentVariables( const TimeType evaluationTime )
    {
        dependentVariables_.setZero( );

        // Set dependent variable.
        if( dependentVariablesInterpolator_ != nullptr )
        {
            dependentVariables_ = dependentVariablesInterpolator_->interpolate( evaluationTime );
        }
        return dependentVariables_;
    }

    //! Function to get the value of a single dependent variable at a given time.
    Eigen::VectorXd getSingleDependentVariable(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const TimeType evaluationTime )
    {
        Eigen::VectorXd dependentVariable = Eigen::VectorXd( 0 );

        // Retrieve ID and index of dependent variable of interest.
        std::string dependentVariableId = getDependentVariableId( dependentVariableSettings );
        if( dependentVariablesIdsAndIndices_.count( dependentVariableId ) > 0 )
        {
            int dependentVariableIndex = dependentVariablesIdsAndIndices_.find( dependentVariableId )->second;

            // Retrieve full vector of dependent variables at a given time.
            Eigen::VectorXd fullDependentVariablesVector = getDependentVariables( evaluationTime );

            dependentVariable = fullDependentVariablesVector.segment( dependentVariableIndex, getDependentVariableSaveSize( dependentVariableSettings ) );
        }
        else
        {
            throw std::runtime_error( "Error, dependent variable " + dependentVariableId + " not found when retrieving parameter." );
        }
        return dependentVariable;
    }

//    //! Function to get the value of a single dependent variable at a given time, from the dependent variable ID.
//    Eigen::VectorXd getSingleDependentVariable(
//        const std::string dependentVariableId,
//        const int dependentVariableSize,
//        const TimeType evaluationTime )
//    {
//        Eigen::VectorXd dependentVariable = Eigen::VectorXd( dependentVariableSize );
//
//        // Retrieve index of dependent variable of interest.
//        int dependentVariableIndex = dependentVariablesIdsAndIndices_.find( dependentVariableId )->second;
//
//        // Retrieve full vector of dependent variables at a given time.
//        Eigen::VectorXd fullDependentVariablesVector = getDependentVariables( evaluationTime );
//
//        dependentVariable =
//            fullDependentVariablesVector.segment( dependentVariableIndex, dependentVariableSize );
//
//        return dependentVariable;
//    }



    //! Function to get the size of the dependent variables
    /*!
     * Function to get the size of the dependent variables
     * \return Size of dependent variables
     */
    int getDependentVariablesize( )
    {
        return dependentVariablesSize_;
    }

    //! Function to retrieve vector of single dependent variables settings.
    /*!
     * Function to retrieve the dependent variables settings object.
     * \return dependent variables settings
     */
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > getDependentVariablesSettings( ) const
    {
        return dependentVariablesSettings_;
    }


    //! Function to retrieve the map containing the dependent variables Ids and indices.
    std::map< std::string, int > getDependentVariablesIdsAndIndices( )
    {
        return dependentVariablesIdsAndIndices_;
    }

private:


    //! Vector of single dependent variable settings objects
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesSettings_;

    //! Interpolator returning the dependent variables as a function of time.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > dependentVariablesInterpolator_;

    std::map< std::pair< int, int >, std::string > dependentVariableIds_;

    std::map< std::pair< int, int >, std::shared_ptr< SingleDependentVariableSaveSettings > > orderedDependentVariableSettings_;

    //! Predefined vector to use as return value when calling getDependentVariables.
    Eigen::VectorXd dependentVariables_;

    //! Type of the dependent variables of interest
    std::vector< PropagationDependentVariables > dependentVariablesTypes_;

    //! Size of dependent variable vector
    int dependentVariablesSize_;

    //! Map containing the dependent variables Ids and indices
    std::map< std::string, int > dependentVariablesIdsAndIndices_;
};


//! Interface object of interpolation of numerically propagated dependent variables for multi-arc
//! estimation.
template< typename TimeType = double >
class MultiArcDependentVariablesInterface: public DependentVariablesInterface< TimeType >
{
public:


    //! Constructor
    /*!
     * Constructor
     * \param dependentVariablesInterpolators Vector of interpolators returning the dependent variables values as a function of time.
     * \param dependentVariablesSettings Vector of single dependent variables settings.
     * \param arcStartTimes Times at which the multiple arcs start
     * \param arcEndTimes Times at which the multiple arcs end
     */
    MultiArcDependentVariablesInterface(
            const std::vector< std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > > singleArcInterfaces,
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& arcEndTimes ):
            DependentVariablesInterface< TimeType >( ),
            singleArcInterfaces_( singleArcInterfaces ),
            arcStartTimes_( arcStartTimes ),
            arcEndTimes_( arcEndTimes )
    {

    }

    //! Destructor
    ~MultiArcDependentVariablesInterface( ){ }

    //! Function to reset the dependent variables interpolators
    /*!
     * Function to reset the dependent variables interpolators
     * \param dependentVariableInterpolators New vector of interpolators returning the dependent variables as a function of time.
     * \param arcStartTimes Times at which the multiple arcs start.
     * \param arcEndTimes Times at which the multiple arcs end
     */
    void updateDependentVariablesInterpolators(
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > > dependentVariablesInterpolators,
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& arcEndTimes  )
    {
        if( arcStartTimes.size( ) != arcEndTimes.size( ) )
        {
            throw std::runtime_error( "Error when updating MultiArcDependentVariablesInterface, incompatible time lists" );
        }

        if( dependentVariablesInterpolators.size( ) != arcEndTimes.size( ) )
        {
            throw std::runtime_error( "Error when updating MultiArcDependentVariablesInterface, incompatible interpolator lists" );
        }

        if( dependentVariablesInterpolators.size( ) != singleArcInterfaces_.size( ) )
        {
            throw std::runtime_error( "Error when updating MultiArcDependentVariablesInterface, size of interpolator lists is incompatible with number of arcs" );

        }

        arcStartTimes_ = arcStartTimes;
        arcEndTimes_ = arcEndTimes;

        for( unsigned int i = 0; i < dependentVariablesInterpolators.size( ); i++ )
        {
            singleArcInterfaces_.at( i )->updateDependentVariablesInterpolator( dependentVariablesInterpolators.at( i ) );
        }
        std::vector< double > arcSplitTimes = arcStartTimes_;
        arcSplitTimes.push_back( std::numeric_limits< double >::max( ) );

        lookUpscheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( arcSplitTimes );
    }

    //! Function to get the vector of interpolators returning the dependent variables as a function of time.
    /*!
     * Function to get the vector of interpolators returning the dependent variables as a function of time.
     * \return Vector of interpolators returning the dependent variables as a function of time.
     */
    std::vector< std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > > getSingleArcInterfaces( )
    {
        return singleArcInterfaces_;
    }

    //! Function to get the concatenated dependent variables values at a given time.
    /*!
     * Function to get the concatenated dependent variables values at a given time.
     *  \param evaluationTime Time at which to evaluate dependent variables interpolators
     *  \return Concatenated dependent variables values.
     */
    Eigen::VectorXd getDependentVariables( const TimeType evaluationTime )
    {
        Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( 0 );
        int currentArc = getCurrentArc( evaluationTime ).first;

        // Set dependent variables vector.
        if( currentArc >= 0 )
        {
            dependentVariables = singleArcInterfaces_.at( currentArc )->getDependentVariables( evaluationTime );
        }

        return dependentVariables;
    }

    //! Function to get the value of a single dependent variable at a given time.
    Eigen::VectorXd getSingleDependentVariable(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const TimeType evaluationTime )
    {
        Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( 0 );
        int currentArc = getCurrentArc( evaluationTime ).first;

        // Set dependent variables vector.
        if( currentArc >= 0 )
        {
            dependentVariables = singleArcInterfaces_.at( currentArc )->getSingleDependentVariable( dependentVariableSettings, evaluationTime );
        }

        return dependentVariables;
    }

    //! Function to retrieve the current arc for a given time
    /*!
     * Function to retrieve the current arc for a given time
     * \param evaluationTime Time at which current arc is to be determined
     * \return Pair with current arc index and associated arc initial time.
     */
    std::pair< int, double > getCurrentArc( const TimeType evaluationTime )
    {
        if( lookUpscheme_ == nullptr )
        {
            throw std::runtime_error( "Error when accessing multi-arc dependent variable interface; interface not yet set" );
        }

        int currentArc =  lookUpscheme_->findNearestLowerNeighbour( evaluationTime );
        if( evaluationTime <= arcEndTimes_.at( currentArc ) && evaluationTime >= arcStartTimes_.at( currentArc ) )
        {
            return std::make_pair( currentArc, arcStartTimes_.at( currentArc ) );
        }
        else
        {
            return std::make_pair( -1, TUDAT_NAN );
        }
    }

    //! Function to retrieve the number of arcs in dynamics
    /*!
     * Function to retrieve the number of arcs in dynamics
     * \return Number of arcs in dynamics
     */
    int getNumberOfArcs( )
    {
        return numberArcs_;
    }


private:

    std::vector< std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > > singleArcInterfaces_;

    //! Times at which the multiple arcs start
    std::vector< double > arcStartTimes_;

    //! Times at which the multiple arcs end
    std::vector< double > arcEndTimes_;

    //! Number of arcs.
    int numberArcs_;

    //! Look-up algorithm to determine the arc of a given time.
    std::shared_ptr< interpolators::HuntingAlgorithmLookupScheme< double > > lookUpscheme_;

};


//! Interface object of interpolation of numerically propagated dependent variables for a hybrid of
//! single-and multi-arc estimation (single order is put first in concatenation)
/*!
*  Interface object of interpolation of numerically propagated dependent variables for a hybrid of
*  single-and multi-arc estimation (single order is put first in concatenation).
*/
template< typename TimeType = double >
class HybridArcDependentVariablesInterface: public DependentVariablesInterface< TimeType >
{
public:

    //! Constructor
    /*!
     * Constructor
     * \param singleArcInterface Object to retrieve the dependent variable for single arc component
     * \param multiArcInterface Object to retrieve the dependent variable for multi arc component
     */
    HybridArcDependentVariablesInterface(
            const std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > singleArcInterface,
            const std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > multiArcInterface ):
            DependentVariablesInterface< TimeType >(  ),
            singleArcInterface_( singleArcInterface ), multiArcInterface_( multiArcInterface )
    {
    }

    //! Destructor
    ~HybridArcDependentVariablesInterface( ){ }

    //! Function to get the dependent variable at a given time.
    /*!
     *  Function to get the dependent variable matrix at a given time.
     *  \param evaluationTime Time at which to evaluate matrix interpolators
     *  \return Dependent variable.
     */
    Eigen::VectorXd getDependentVariables( const TimeType evaluationTime )
    {
        throw std::runtime_error( "Error when retrieving interpolated dependent variables from hybrid-arc interface. This functionality is not supported as the definition is ambiguous. Retreieve the single- or multi-arc interface, and retrieve the dependent variable from there." );
        return Eigen::VectorXd::Zero( 0 );
    }

    Eigen::VectorXd getSingleDependentVariable(
        const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
        const TimeType evaluationTime )
    {
        Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( 0 );

        // Get single-arc dependent variables.
        Eigen::VectorXd singleArcDependentVariables = singleArcInterface_->getSingleDependentVariable( dependentVariableSettings, evaluationTime );

        // Get multi-arc dependent variables.
        Eigen::VectorXd multiArcDependentVariables = multiArcInterface_->getSingleDependentVariable( dependentVariableSettings, evaluationTime );

        // Single-arc result takes preference over multi-arc result
        if( singleArcDependentVariables.rows( ) > 0 )
        {
            dependentVariables = singleArcDependentVariables;
        }
        else if( multiArcDependentVariables.rows( ) > 0 )
        {
            dependentVariables = multiArcDependentVariables;
        }

        return dependentVariables;
    }
private:

    //! Object to retrieve dependent variable for single arc component
    std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > singleArcInterface_;

    //! Object to retrieve dependent variable for multi arc component
    std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > multiArcInterface_;
};






    } // namespace propagators

} // namespace tudat

#endif // TUDAT_DEPENDENTVARIABLESINTERFACE_H
