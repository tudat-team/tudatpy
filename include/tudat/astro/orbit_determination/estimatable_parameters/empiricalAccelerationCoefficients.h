/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_EMPIRICALACCELERATIONCOEFFICIENTS_H
#define TUDAT_EMPIRICALACCELERATIONCOEFFICIENTS_H

#include "tudat/astro/orbit_determination/estimatable_parameters/estimatableParameter.h"
#include "tudat/astro/basic_astro/empiricalAcceleration.h"
#include "tudat/math/interpolators/piecewiseConstantInterpolator.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for estimation of a body's time-independent empirical accelerations
/*!
 * Interface class for estimation of a body's time-independent empirical accelerations. Interfaces the estimation with the
 * acceleration components in the EmpiricalAcceleration class
 */
class EmpiricalAccelerationCoefficientsParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param empiricalAcceleration Class defining properties of empirical acceleration used in propagation.
     *  \param associatedBody Body for which empirical accelerations are estimated
     *  \param centralBody Central body used for empirical accelerations
     *  \param componentsToEstimate List of empirical acceleration components that are to be estimated. Map key is functional
     *  shape of acceleration, value is list of components.
     */
    EmpiricalAccelerationCoefficientsParameter(
            const std::vector< std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > > empiricalAcceleration,
            const std::string& associatedBody,
            const std::string& centralBody,
            const std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate ):
        EstimatableParameter< Eigen::VectorXd >( empirical_acceleration_coefficients, associatedBody, centralBody ),
        empiricalAcceleration_( empiricalAcceleration )
    {
        parameterSize_ = 0;

        // Set components of empirical accelerations that are to be estimated
        for( std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
             std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > >::
             const_iterator inputComponentIterator =  componentsToEstimate.begin( );
             inputComponentIterator != componentsToEstimate.end( );
             inputComponentIterator++ )
        {
            int currentComponentIndex = static_cast< int >( inputComponentIterator->first );

            // Iterate over all components
            for( unsigned int i = 0; i < inputComponentIterator->second.size( ); i++ )
            {
                // Check feasibility of input
                if( accelerationIndices_.count( inputComponentIterator->second.at( i ) ) > 0 )
                {
                    std::vector< int > currentIndices = accelerationIndices_.at( inputComponentIterator->second.at( i ) );
                    if( std::find( currentIndices.begin( ), currentIndices.end( ), currentComponentIndex ) !=
                            currentIndices.end( ) )
                    {
                        throw std::runtime_error(
                                    "Error when creating empirical acceleration parameter object, found duplicate of component." );
                    }
                }
                accelerationIndices_[ inputComponentIterator->second.at( i ) ].push_back( currentComponentIndex );
                parameterSize_++;
            }
        }

    }

    //! Destructor
    ~EmpiricalAccelerationCoefficientsParameter( ) { }

    //! Get value of empirical acceleration components
    /*!
     *  Get value of empirical acceleration components
     *  \return Value of empirical acceleration components
     */
    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd parameter = Eigen::VectorXd::Zero( parameterSize_ );

        // Iterate over all acceleration models
        for( unsigned int i = 0; i < empiricalAcceleration_.size( ); i++ )
        {

            Eigen::VectorXd currentAccelerationParameter = Eigen::VectorXd::Zero( parameterSize_ );

            // Iterate over all components and retrieve required values
            Eigen::Matrix3d accelerationComponents = empiricalAcceleration_.at( i )->getAccelerationComponents( );
            int currentIndex = 0;
            for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator
                 componentIterator = accelerationIndices_.begin( );
                 componentIterator != accelerationIndices_.end( ); componentIterator++ )
            {
                Eigen::Vector3d accelerations =
                        accelerationComponents.block( 0, static_cast< int >( componentIterator->first ), 3, 1 );

                for( unsigned int i = 0; i < componentIterator->second.size( ); i++ )
                {
                    currentAccelerationParameter( currentIndex ) = accelerations( componentIterator->second.at( i ) );
                    currentIndex++;
                }
            }

            // Check consistency
            if( currentIndex != parameterSize_ )
            {
                throw std::runtime_error( "Error when getting empirical parameter; inconsistent sizes found." );
            }
            if( i == 0 )
            {
                parameter = currentAccelerationParameter;
            }
            else if( parameter != currentAccelerationParameter )
            {
                std::cerr<<"Warning when retrieving emp. acc. parameter from list. List entries are incompatible. First entry is: "<<std::endl<<
                           parameter<<std::endl<<" Current value is: "<<std::endl<<currentAccelerationParameter<<std::endl<<
                           "Difference is: "<<parameter - currentAccelerationParameter<<std::endl;
            }
        }

        return parameter;
    }

    //! Reset value of empirical acceleration components
    /*!
     *  Reset value of empirical acceleration components
     *  \param parameterValue New value of empirical acceleration components
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        for( unsigned int i = 0; i < empiricalAcceleration_.size( ); i++ )
        {
            int currentIndex = 0;

            // Iterate over all components and set required values
            Eigen::Matrix3d accelerationComponents = empiricalAcceleration_.at( i )->getAccelerationComponents( );
            for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator
                 componentIterator = accelerationIndices_.begin( );
                 componentIterator != accelerationIndices_.end( ); componentIterator++ )
            {
                Eigen::Vector3d accelerations =
                        accelerationComponents.block( 0, static_cast< int >( componentIterator->first ), 3, 1 );
                for( unsigned int i = 0; i < componentIterator->second.size( ); i++ )
                {
                    accelerations( componentIterator->second.at( i ) ) = parameterValue( currentIndex );
                    currentIndex++;
                }
                accelerationComponents.block( 0, static_cast< int >( componentIterator->first ), 3, 1 ) = accelerations;
            }

            // Reset components in acceleration model
            empiricalAcceleration_.at( i )->resetAccelerationComponents( accelerationComponents );

            if( currentIndex != parameterSize_ )
            {
                throw std::runtime_error(
                            "Error when getting empirical parameter size; inconsistent sizes found." );
            }
        }
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value.
     */
    int getParameterSize( )
    {
        return parameterSize_;
    }

    std::string getParameterDescription( )
    {
        std::string parameterDescription = ", components in RSW frame, functional shapes; ";

        for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator
             indexIterator = accelerationIndices_.begin( ); indexIterator != accelerationIndices_.end( ); indexIterator++ )
        {
            parameterDescription += "("+ basic_astrodynamics::getEmpiricalAccelerationFunctionalShapeString( indexIterator->first ) + ": index ";
            for( unsigned int i = 0; i < indexIterator->second.size( ); i++ )
            {
                parameterDescription += std::to_string( indexIterator->second.at( i ) );
                if( i != indexIterator->second.size( ) - 1 )
                {
                    parameterDescription += ", ";
                }
            }
            parameterDescription += ")";
        }

        return parameterDescription;
    }

    //! Function to retrieve list of components in empirical accelerations that are to be estimated.
    /*!
     * Function to retrieve list of components in empirical accelerations that are to be estimated.
     * \return List of components in empirical accelerations that are to be estimated.
     */
    std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > > getIndices( )
    {
        return accelerationIndices_;
    }

protected:

private:

    //! Class defining properties of empirical acceleration used in propagation.
    std::vector< std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > > empiricalAcceleration_;

    //! Number of empirical acceleration components that are to be estimated
    int parameterSize_;

    //! List of components in empirical accelerations that are to be estimated.
    std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > > accelerationIndices_;
};

//! Interface class for estimation of a body's time-dependent (arcwise constant) empirical accelerations
/*!
 * Interface class for estimation of a body's time-dependent (arcwise constant) empirical accelerations.
 * Interfaces the estimation with the acceleration components in the EmpiricalAcceleration class
 */
class ArcWiseEmpiricalAccelerationCoefficientsParameter: public EstimatableParameter< Eigen::VectorXd >
{

public:

    //! Constructor
    /*!
     *  Constructor
     *  \param empiricalAcceleration Class defining properties of empirical acceleration used in propagation.
     *  \param associatedBody Body for which empirical accelerations are estimated
     *  \param centralBody Central body used for empirical accelerations
     *  \param componentsToEstimate List of empirical acceleration components that are to be estimated. Map key is functional
     *  shape of acceleration, value is list of components.
     *  \param arcStartTimeList List of times at which the arcs over which empirical accelerations are constant are to be
     *  estimated
     */
    ArcWiseEmpiricalAccelerationCoefficientsParameter(
            const std::vector< std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > > empiricalAcceleration,
            const std::string& associatedBody,
            const std::string& centralBody,
            const std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
            std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > > componentsToEstimate,
            const std::vector< double > arcStartTimeList ):
        EstimatableParameter< Eigen::VectorXd >( arc_wise_empirical_acceleration_coefficients, associatedBody, centralBody ),
        empiricalAcceleration_( empiricalAcceleration ),
        arcStartTimeList_( arcStartTimeList )
    {
        singleArcParameterSize_ = 0;

        // Set components of empirical accelerations that are to be estimated
        for( std::map< basic_astrodynamics::EmpiricalAccelerationComponents,
             std::vector< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes > >::
             const_iterator inputComponentIterator =  componentsToEstimate.begin( );
             inputComponentIterator != componentsToEstimate.end( );
             inputComponentIterator++ )
        {
            int currentComponentIndex = static_cast< int >( inputComponentIterator->first );

            // Iterate over all components
            for( unsigned int i = 0; i < inputComponentIterator->second.size( ); i++ )
            {
                // Check feasibility of input
                if( accelerationIndices_.count( inputComponentIterator->second.at( i ) ) > 0 )
                {
                    std::vector< int > currentIndices = accelerationIndices_.at( inputComponentIterator->second.at( i ) );
                    if( std::find( currentIndices.begin( ), currentIndices.end( ), currentComponentIndex ) !=
                            currentIndices.end( ) )
                    {
                        throw std::runtime_error(
                                    "Error when creating arcwise empirical acceleration parameter object, found duplicate of component." );
                    }
                }
                accelerationIndices_[ inputComponentIterator->second.at( i ) ].push_back( currentComponentIndex );
                singleArcParameterSize_++;
            }
        }

        // Add maximum time to end of arc list.
        arcStartTimeList_.push_back( 1.0E300 );

        // Retrieve current empirical accelerations (set in each arc)
        Eigen::Matrix3d currentTimeInvariantEmpiricalAccelerations = empiricalAcceleration.at( 0 )->getAccelerationComponents( );
        for( unsigned int i = 1; i < empiricalAcceleration.size( ); i++ )
        {
            Eigen::Matrix3d comparisonAccelerations = empiricalAcceleration.at( i )->getAccelerationComponents( );
            if( comparisonAccelerations != currentTimeInvariantEmpiricalAccelerations )
            {
                std::cerr<<"Warning when initializing arc-wise empirical acceleration parameter, list of input acceleration models do not have identical accelerations."<<std::endl;
            }

        }
        for( unsigned int i = 0; i < arcStartTimeList_.size( ); i++ )
        {
            empiricalAccelerationList_.push_back( currentTimeInvariantEmpiricalAccelerations );
        }

        // Create piecewise constant interpolator for empirical acceleration
        empiricalAccelerationInterpolator_ = std::make_shared<
                interpolators::PiecewiseConstantInterpolator< double, Eigen::Matrix3d > >(
                    arcStartTimeList_, empiricalAccelerationList_ );
        typedef interpolators::OneDimensionalInterpolator< double, Eigen::Matrix3d > LocalInterpolator;

        empiricalAccelerationFunction_ =
                std::bind( static_cast< Eigen::Matrix3d( LocalInterpolator::* )( const double ) >
                           ( &LocalInterpolator::interpolate ), empiricalAccelerationInterpolator_, std::placeholders::_1 );

    }

    //! Destructor
    ~ArcWiseEmpiricalAccelerationCoefficientsParameter( ) { }

    //! Get value of arcwise empirical acceleration components
    /*!
     *  Get value of arcwise empirical acceleration components
     *  \return Value of arcwise empirical acceleration components
     */
    Eigen::VectorXd getParameterValue( )
    {
        Eigen::VectorXd parameter = Eigen::VectorXd::Zero( getParameterSize( ) );
        int currentIndex = 0;

        // Iterate over all arcs and retrieve required values
        for( unsigned int i = 0; i < arcStartTimeList_.size( ) - 1; i++ )
        {
            // Iterate over all components and retrieve required values
            Eigen::Matrix3d currentArcAccelerations = empiricalAccelerationList_.at( i );
            for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator
                 componentIterator = accelerationIndices_.begin( );
                 componentIterator != accelerationIndices_.end( ); componentIterator++ )
            {
                int currentBlockColumn = static_cast< int >( componentIterator->first );
                for( unsigned int j = 0; j < componentIterator->second.size( ); j++ )
                {
                    parameter( currentIndex ) = currentArcAccelerations( componentIterator->second.at( j ), currentBlockColumn );
                    currentIndex++;
                }
            }
        }

        // Check consistency
        if( currentIndex != getParameterSize( ) )
        {
            throw std::runtime_error( "Error when getting arcwise empirical parameter; inconsistent sizes found." );

        }

        return parameter;
    }

    //! Reset value of arcwise empirical acceleration components
    /*!
     *  Reset value of arcwise empirical acceleration components
     *  \param parameterValue New value of empirical acceleration components
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        for( unsigned int i = 0; i < empiricalAcceleration_.size( ); i++ )
        {
            int currentIndex = 0;

            // Iterate over all arcs and set required values
            for( unsigned int i = 0; i < arcStartTimeList_.size( ) - 1; i++ )
            {
                // Iterate over all components and set required values
                Eigen::Matrix3d currentArcAccelerations = empiricalAccelerationList_.at( i );
                for( std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > >::const_iterator
                     componentIterator = accelerationIndices_.begin( );
                     componentIterator != accelerationIndices_.end( ); componentIterator++ )
                {
                    int currentBlockColumn = static_cast< int >( componentIterator->first );
                    for( unsigned int j = 0; j < componentIterator->second.size( ); j++ )
                    {
                        currentArcAccelerations( componentIterator->second.at( j ), currentBlockColumn ) =
                                parameterValue( currentIndex );
                        currentIndex++;
                    }
                }
                empiricalAccelerationList_[ i ] = currentArcAccelerations;
            }

            // Set list of arc times and interpolation function
            empiricalAccelerationInterpolator_->resetDependentValues( empiricalAccelerationList_ );
            empiricalAcceleration_.at( i )->resetAccelerationComponentsFunction( empiricalAccelerationFunction_ );

            if( currentIndex != getParameterSize( ) )
            {
                throw std::runtime_error(
                            "Error when setting arc-wise empirical parameter size; inconsistent sizes found." );
            }
        }
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value.
     */
    int getParameterSize( )
    {
        return singleArcParameterSize_ * ( arcStartTimeList_.size( ) - 1 );
    }

    //! Function to retrieve the size of the parameter for a single arc
    /*!
     *  Function to retrieve the size of the parameter for a single arc
     *  \return Size of parameter value for a single arc
     */
    int getSingleArcParameterSize( )
    {
        return singleArcParameterSize_;
    }

    std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > > getIndices( )
    {
        return accelerationIndices_;
    }

    std::shared_ptr< interpolators::LookUpScheme< double > > getArcTimeLookupScheme( )
    {
        return empiricalAccelerationInterpolator_->getLookUpScheme( );
    }

protected:

private:

    //! Class defining properties of empirical acceleration used in propagation.
    std::vector< std::shared_ptr< basic_astrodynamics::EmpiricalAcceleration > > empiricalAcceleration_;

    //! Interpolator to compute the empirical acceleration as a function of time
    std::shared_ptr< interpolators::PiecewiseConstantInterpolator< double, Eigen::Matrix3d > >
    empiricalAccelerationInterpolator_;

    //! Function returning the empirical acceleration as a function of time (linked to empiricalAccelerationInterpolator_)
    std::function< Eigen::Matrix3d( const double ) > empiricalAccelerationFunction_;

    //! List of components in empirical accelerations that are to be estimated for every arc.
    std::map< basic_astrodynamics::EmpiricalAccelerationFunctionalShapes, std::vector< int > > accelerationIndices_;

    //! List of current empirical accelerations per arc
    std::vector< Eigen::Matrix3d > empiricalAccelerationList_;

    //! List of times at which arcs start (and final time of last of as last entry).
    std::vector< double > arcStartTimeList_;

    //! Number of empirical acceleration components that are to be estimated for a single arc
    int singleArcParameterSize_;

};

}

}

#endif // TUDAT_EMPIRICALACCELERATIONCOEFFICIENTS_H
