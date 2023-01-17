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
    DependentVariablesInterface(
            const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesSettings ):
            dependentVariablesSettings_( dependentVariablesSettings )
    {
        dependentVariablesSize_ = 0;

        dependentVariablesIdsAndIndices_.clear( );

        for ( unsigned int i= 0 ; i < dependentVariablesSettings_.size( ) ; i++ )
        {
            dependentVariablesTypes_.push_back( dependentVariablesSettings_[ i ]->dependentVariableType_ );

            dependentVariablesIdsAndIndices_[ getDependentVariableId( dependentVariablesSettings_[ i ] ) ] = dependentVariablesSize_;

            dependentVariablesSize_ += getDependentVariableSaveSize( dependentVariablesSettings_[ i ] );
        }
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

    //! Function to get the value of a single dependent variable at a given time.
    Eigen::VectorXd getSingleDependentVariable(
            const std::shared_ptr< SingleDependentVariableSaveSettings > dependentVariableSettings,
            const TimeType evaluationTime )
    {
        Eigen::VectorXd dependentVariable = Eigen::VectorXd( getDependentVariableSaveSize( dependentVariableSettings ) );

        // Retrieve ID and index of dependent variable of interest.
        std::string dependentVariableId = getDependentVariableId( dependentVariableSettings );
        int dependentVariableIndex = dependentVariablesIdsAndIndices_.find( dependentVariableId )->second;

        // Retrieve full vector of dependent variables at a given time.
        Eigen::VectorXd fullDependentVariablesVector = getDependentVariables( evaluationTime );

        dependentVariable =
                fullDependentVariablesVector.segment( dependentVariableIndex, getDependentVariableSaveSize( dependentVariableSettings ) );

        return dependentVariable;
    }

    //! Function to get the value of a single dependent variable at a given time, from the dependent variable ID.
    Eigen::VectorXd getSingleDependentVariable(
            const std::string dependentVariableId,
            const int dependentVariableSize,
            const TimeType evaluationTime )
    {
        Eigen::VectorXd dependentVariable = Eigen::VectorXd( dependentVariableSize );

        // Retrieve index of dependent variable of interest.
        int dependentVariableIndex = dependentVariablesIdsAndIndices_.find( dependentVariableId )->second;

        // Retrieve full vector of dependent variables at a given time.
        Eigen::VectorXd fullDependentVariablesVector = getDependentVariables( evaluationTime );

        dependentVariable =
                fullDependentVariablesVector.segment( dependentVariableIndex, dependentVariableSize );

        return dependentVariable;
    }



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
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > getDependentVariablesSettings( )
    {
        return dependentVariablesSettings_;
    }


    //! Function to retrieve the map containing the dependent variables Ids and indices.
    std::map< std::string, int > getDependentVariablesIdsAndIndices( )
    {
        return dependentVariablesIdsAndIndices_;
    }


protected:

    //! Vector of single dependent variable settings objects
    std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesSettings_;

    //! Type of the dependent variables of interest
    std::vector< PropagationDependentVariables > dependentVariablesTypes_;

    //! Size of dependent variable vector
    int dependentVariablesSize_;

    //! Map containing the dependent variables Ids and indices
    std::map< std::string, int > dependentVariablesIdsAndIndices_;

};

//! Interface object of interpolation of numerically propagated dependent variables for single-arc propagation/estimation.
template< typename TimeType = double >
class SingleArcDependentVariablesInterface : public DependentVariablesInterface< TimeType >
{
public:

    using DependentVariablesInterface< TimeType >::dependentVariablesSize_;

    //! Constructor
    /*!
     * Constructor
     * \param dependentVariablesInterpolator Interpolator returning the dependent variables values as a function of time.
     * \param dependentVariablesSettings Vector of single dependent variables settings
     */
    SingleArcDependentVariablesInterface(
            const std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > dependentVariablesInterpolator,
            const std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > dependentVariablesSettings ):
            DependentVariablesInterface< TimeType >( dependentVariablesSettings ),
            dependentVariablesInterpolator_( dependentVariablesInterpolator )
    {
        dependentVariables_ = Eigen::VectorXd::Zero( dependentVariablesSize_ );
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
        dependentVariables_ = dependentVariablesInterpolator_->interpolate( evaluationTime );
        return dependentVariables_;
    }


private:

    //! Predefined vector to use as return value when calling getDependentVariables.
    Eigen::VectorXd dependentVariables_;

    //! Interpolator returning the dependent variables as a function of time.
    std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > dependentVariablesInterpolator_;

};


//! Interface object of interpolation of numerically propagated dependent variables for multi-arc
//! estimation.
template< typename TimeType = double >
class MultiArcDependentVariablesInterface: public DependentVariablesInterface< TimeType >
{
public:

    using DependentVariablesInterface< TimeType >::dependentVariablesSize_;

    //! Constructor
    /*!
     * Constructor
     * \param dependentVariablesInterpolators Vector of interpolators returning the dependent variables values as a function of time.
     * \param dependentVariablesSettings Vector of single dependent variables settings.
     * \param arcStartTimes Times at which the multiple arcs start
     * \param arcEndTimes Times at which the multiple arcs end
     */
    MultiArcDependentVariablesInterface(
            const std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > > dependentVariablesInterpolators,
            const std::vector< std::vector< std::shared_ptr< SingleDependentVariableSaveSettings > > > dependentVariablesSettings,
            const std::vector< double >& arcStartTimes,
            const std::vector< double >& arcEndTimes ):
            DependentVariablesInterface< TimeType >( dependentVariablesSettings.at( 0 ) ),
            dependentVariablesInterpolators_( dependentVariablesInterpolators ),
            arcStartTimes_( arcStartTimes ),
            arcEndTimes_( arcEndTimes )
    {
        if( arcStartTimes_.size( ) != arcEndTimes_.size( ) )
        {
            throw std::runtime_error( "Error when making MultiArcDependentVariablesInterface, incompatible time lists" );
        }
        numberArcs_ = arcStartTimes_.size( );

        std::vector< double > arcSplitTimes = arcStartTimes_;
        arcSplitTimes.push_back(  std::numeric_limits< double >::max( ));
        lookUpscheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( arcSplitTimes );
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
        dependentVariablesInterpolators_ = dependentVariablesInterpolators;
        arcStartTimes_ = arcStartTimes;
        arcEndTimes_ = arcEndTimes;

        std::vector< double > arcSplitTimes = arcStartTimes_;
        arcSplitTimes.push_back( std::numeric_limits< double >::max( ) );

        lookUpscheme_ = std::make_shared< interpolators::HuntingAlgorithmLookupScheme< double > >( arcSplitTimes );
    }

    //! Function to get the vector of interpolators returning the dependent variables as a function of time.
    /*!
     * Function to get the vector of interpolators returning the dependent variables as a function of time.
     * \return Vector of interpolators returning the dependent variables as a function of time.
     */
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > > getDependentVariablesInterpolators( )
    {
        return dependentVariablesInterpolators_;
    }

    //! Function to get the concatenated dependent variables values at a given time.
    /*!
     * Function to get the concatenated dependent variables values at a given time.
     *  \param evaluationTime Time at which to evaluate dependent variables interpolators
     *  \return Concatenated dependent variables values.
     */
    Eigen::VectorXd getDependentVariables( const TimeType evaluationTime )
    {
        Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( dependentVariablesSize_ );

        int currentArc = getCurrentArc( evaluationTime ).first;

        // Set dependent variables vector.
        if( currentArc >= 0 )
        {
            dependentVariables.segment( 0, dependentVariablesSize_) = dependentVariablesInterpolators_.at( currentArc )->interpolate( evaluationTime );
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

    //! List of interpolators returning the dependent variables as a function of time.
    std::vector< std::shared_ptr< interpolators::OneDimensionalInterpolator< TimeType, Eigen::VectorXd > > > dependentVariablesInterpolators_;

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
            DependentVariablesInterface< TimeType >( multiArcInterface->getDependentVariablesSettings( ) ),
            singleArcInterface_( singleArcInterface ), multiArcInterface_( multiArcInterface )
    {
        singleArcDependentVariablesSize_ = 0;
        multiArcDependentVariablesSize_ = 0;
        numberArcs_ = 0;

        if ( singleArcInterface_ != nullptr )
        {
            singleArcDependentVariablesIdsAndIndices_ = singleArcInterface_->getDependentVariablesIdsAndIndices( );
            singleArcDependentVariablesSize_ = singleArcInterface_->getDependentVariablesize( );
        }
        std::map< std::string, int > multiArcDependentVariablesIdsAndIndices;
        if ( multiArcInterface_ != nullptr )
        {
            multiArcDependentVariablesIdsAndIndices = multiArcInterface_->getDependentVariablesIdsAndIndices( );
            multiArcDependentVariablesSize_ = multiArcInterface_->getDependentVariablesize( ) /*- singleArcDependentVariablesSize_*/;
            numberArcs_ = multiArcInterface->getNumberOfArcs( );
        }

        // Check input consistency: verify that the same dependent variable is not saved twice (in both the single and multi-arc
        // dependent variables interface)
        if ( ( singleArcInterface != nullptr ) && ( multiArcInterface != nullptr ) )
        {
            for ( std::map< std::string, int >::iterator itr = multiArcDependentVariablesIdsAndIndices.begin( ) ;
                  itr != multiArcDependentVariablesIdsAndIndices.end( ) ; itr++ )
            {
                if ( singleArcDependentVariablesIdsAndIndices_.count( itr->first ) != 0 )
                {
                    // Remove dependent variable from multi-arc dependent variables interface.
                    idsAndIndicesMultiArcDependentVariablesToBeRemoved_[ itr->first ] = itr->second;

                    std::cerr << "Warning when making hybrid dependent variables interface, dependent variable "
                              << itr->first << " required by the multi-arc interface is already included in the single arc interface." << std::endl;
                }
            }
        }

        int sizeRemovedDependentVariables = 0;
        int counterDependentVariable = 0;
        for ( std::map< std::string, int >::iterator itr = multiArcDependentVariablesIdsAndIndices.begin( ) ;
              itr != multiArcDependentVariablesIdsAndIndices.end( ) ; itr++ )
        {
            if ( idsAndIndicesMultiArcDependentVariablesToBeRemoved_.count( itr->first ) == 0 )
            {
                multiArcDependentVariablesIdsAndIndices_[ itr->first ] = itr->second - sizeRemovedDependentVariables;
            }
            else
            {
                sizeRemovedDependentVariables += getDependentVariableSaveSize(
                        multiArcInterface_->getDependentVariablesSettings( )[ counterDependentVariable ] );
                multiArcDependentVariablesSize_ -= getDependentVariableSaveSize(
                        multiArcInterface_->getDependentVariablesSettings( )[ counterDependentVariable ] );
            }
            counterDependentVariable += 1;
        }

        hybridArcDependentVariablesSize_ = singleArcDependentVariablesSize_ + multiArcDependentVariablesSize_;

//        // Reset dependent variables IDs and indices map.
//        dependentVariablesIdsAndIndices_.clear( );
//        dependentVariablesIdsAndIndices_ = singleArcDependentVariablesIdsAndIndices_;
//        for ( std::map< std::string, int >::iterator itr = multiArcDependentVariablesIdsAndIndices_.begin( ) ;
//              itr != multiArcDependentVariablesIdsAndIndices_.end( ) ; itr++ )
//        {
//            dependentVariablesIdsAndIndices_[ itr->first ] = itr->second + singleArcDependentVariablesSize_;
//        }
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
        Eigen::VectorXd dependentVariables = Eigen::VectorXd::Zero( hybridArcDependentVariablesSize_ );

        // Get single-arc dependent variables.
        Eigen::VectorXd singleArcDependentVariables;
        if ( singleArcInterface_ != nullptr )
        {
            singleArcDependentVariables = singleArcInterface_->getDependentVariables( evaluationTime );
        }
        else
        {
            singleArcDependentVariables = Eigen::VectorXd::Zero( 0 );
        }

        // Get multi-arc dependent variables.
        Eigen::VectorXd multiArcDependentVariables = multiArcInterface_->getDependentVariables( evaluationTime );
        std::pair< int, double >  currentArc = multiArcInterface_->getCurrentArc( evaluationTime );


//    if ( ( multiArcDependentVariables.size( ) ) != hybridArcDependentVariablesSize_ )
//    {
//        throw std::runtime_error( "Error when getting dependent variables from hybrid arc interfaces, size inconsistent "
//                                  "with size of single- and multi-arc dependent variables." );
//    }

        dependentVariables.segment( 0, singleArcDependentVariablesSize_ ) = singleArcDependentVariables;

        if( !( currentArc.first < 0 ) )
        {
            dependentVariables.segment( singleArcDependentVariablesSize_, multiArcDependentVariablesSize_ ) = multiArcDependentVariables;

//        int currentMultiArcDependentVariableIndex = singleArcDependentVariablesSize_;
//
//        for ( std::map< std::string, int >::iterator itr = multiArcDependentVariablesIdsAndIndices_.begin( ) ;
//              itr != multiArcDependentVariablesIdsAndIndices_.end( ) ; itr++ )
//        {
//            int sizeCurrentDependentVariable = 0;
//            for ( unsigned int i = 0 ;
//                  i < multiArcInterface_->getDependentVariablesSettings( ).size( ) ; i++ )
//            {
//                if ( itr->first == getDependentVariableId( multiArcInterface_->getDependentVariablesSettings( )[ i ] ) )
//                {
//                    sizeCurrentDependentVariable = getDependentVariableSaveSize( multiArcInterface_->getDependentVariablesSettings( )[ i ] );
//                }
//            }
//
//            dependentVariables.segment( currentMultiArcDependentVariableIndex, sizeCurrentDependentVariable ) =
//                    multiArcDependentVariables.segment( itr->second, sizeCurrentDependentVariable );
//
//            currentMultiArcDependentVariableIndex += sizeCurrentDependentVariable;
//        }
        }

        return dependentVariables;
    }

private:

    //! Object to retrieve dependent variable for single arc component
    std::shared_ptr< SingleArcDependentVariablesInterface< TimeType > > singleArcInterface_;

    //! Object to retrieve dependent variable for multi arc component
    std::shared_ptr< MultiArcDependentVariablesInterface< TimeType > > multiArcInterface_;

    //! Dependent variables IDs and indices for single-arc interface.
    std::map< std::string, int > singleArcDependentVariablesIdsAndIndices_;

    //! Dependent variables IDs and indices for multi-arc interface.
    std::map< std::string, int > multiArcDependentVariablesIdsAndIndices_;

    //! IDs and indices of dependent variables to be removed from multi-arc interface (in case already included in single arc
    //! interface)
    std::map< std::string, int > idsAndIndicesMultiArcDependentVariablesToBeRemoved_;

    //! Size of single-arc dependent variables
    int singleArcDependentVariablesSize_;

    //! Size of multi-arc dependent variables.
    int multiArcDependentVariablesSize_;

    //! Size of hybrid arc dependent variables.
    int hybridArcDependentVariablesSize_;

    //! Number of arcs in multi-arc model.
    int numberArcs_;
};






    } // namespace propagators

} // namespace tudat

#endif // TUDAT_DEPENDENTVARIABLESINTERFACE_H
