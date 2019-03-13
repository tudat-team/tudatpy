/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CONSTANTDRAGCOEFFICIENT_H
#define TUDAT_CONSTANTDRAGCOEFFICIENT_H

#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/Astrodynamics/OrbitDetermination/EstimatableParameters/estimatableParameter.h"
#include "Tudat/Mathematics/Interpolators/piecewiseConstantInterpolator.h"

namespace tudat
{

namespace estimatable_parameters
{

//! Interface class for the estimation of a constant body drag coefficient
class ConstantDragCoefficient: public EstimatableParameter< double >
{

public:

    //! Constructor
    /*!
     * Constructor
     * \param coefficientInterface Object that contains the aerodynamic coefficients. Constructor checks whether it is
     * consistent with this class, e.g. if the existing aerodynamic coefficients are constant
     * \param associatedBody Body for which the drag coefficient is considered.
     */
    ConstantDragCoefficient(
            const std::shared_ptr< aerodynamics::CustomAerodynamicCoefficientInterface > coefficientInterface,
            const std::string& associatedBody ):
        EstimatableParameter< double >( constant_drag_coefficient, associatedBody ),
        coefficientInterface_( coefficientInterface )
    {
        if( coefficientInterface->getNumberOfIndependentVariables( ) != 0 )
        {
            throw std::runtime_error( "Error when making ConstantDragCoefficient, coefficient interface is inconsistent" );
        }
    }

    //! Destructor
    ~ConstantDragCoefficient( ) { }

    //! Function to get the current value of the constant drag coefficient that is to be estimated.
    /*!
     * Function to get the current value of the constant drag coefficient that is to be estimated.
     * \return Current value of the constant drag coefficient that is to be estimated.
     */
    double getParameterValue( )
    {
        return coefficientInterface_->getConstantCoefficients( )( 0 );
    }

    //! Function to reset the value of the constant drag coefficient that is to be estimated.
    /*!
     * Function to reset the value of the constant drag coefficient that is to be estimated.
     * \param parameterValue New value of the constant drag coefficient that is to be estimated.
     */
    void setParameterValue( double parameterValue )
    {
        Eigen::Vector6d currentCoefficientSet =
                coefficientInterface_->getCurrentAerodynamicCoefficients( );
        currentCoefficientSet( 0 ) = parameterValue;
        coefficientInterface_->resetConstantCoefficients( currentCoefficientSet );
    }

    //! Function to retrieve the size of the parameter (always 1).
    /*!
     *  Function to retrieve the size of the parameter (always 1).
     *  \return Size of parameter value (always 1).
     */
    int getParameterSize( )
    {
        return 1;
    }

protected:

private:

    //! Object that contains the aerodynamic coefficients
    std::shared_ptr< aerodynamics::CustomAerodynamicCoefficientInterface > coefficientInterface_;

};

//! Interface class for the estimation of an arc-wise (piecewise constant) drag coefficient
class ArcWiseConstantDragCoefficient: public EstimatableParameter< Eigen::VectorXd >
{

public:
    //! Constructor.
    /*!
     * Constructor
     * \param coefficientInterface Object containing the drag coefficient to be estimated.
     * \param timeLimits Times at which the arcs are to start.
     * \param associatedBody Name of body containing the coefficientInterface object
     */
    ArcWiseConstantDragCoefficient(
            const std::shared_ptr< aerodynamics::CustomAerodynamicCoefficientInterface > coefficientInterface,
            const std::vector< double > timeLimits,
            const std::string& associatedBody ):
        EstimatableParameter< Eigen::VectorXd >( arc_wise_constant_drag_coefficient, associatedBody ),
        coefficientInterface_( coefficientInterface ), timeLimits_( timeLimits )
    {
        if( coefficientInterface->getNumberOfIndependentVariables( ) != 0 )
        {
            throw std::runtime_error( "Error when making ArcWiseConstantDragCoefficient, coefficient interface is inconsistent" );
        }

        Eigen::Vector6d aerodynamicCoefficients = coefficientInterface->getConstantCoefficients( );
        for( unsigned int i = 0; i < timeLimits.size( ); i++ )
        {
            fullAerodynamicCoefficients_.push_back( aerodynamicCoefficients );
            dragCoefficients_.push_back( aerodynamicCoefficients( 0 ) );
        }

        timeLimits_.push_back( std::numeric_limits< double >::max( ) );
        fullAerodynamicCoefficients_.push_back( aerodynamicCoefficients );

        coefficientInterpolator_ = std::make_shared< interpolators::PiecewiseConstantInterpolator< double, Eigen::Vector6d > >(
                    timeLimits_, fullAerodynamicCoefficients_ );

        setCoefficientInterfaceClosure( );
    }

    //! Destructor.
    ~ArcWiseConstantDragCoefficient( ) { }

    //! Function to get the current value of the arc-wise drag coefficient that is to be estimated.
    /*!
     * Function to get the current value of the arc-wise drag coefficient that is to be estimated.
     * \return Current value of the arc-wise drag coefficients that is to be estimated.
     */
    Eigen::VectorXd getParameterValue( )
    {
        return utilities::convertStlVectorToEigenVector( dragCoefficients_ );
    }

    //! Function to reset the value of the arc-wise drag coefficient that is to be estimated.
    /*!
     * Function to reset the value of the arc-wise drag coefficient that is to be estimated.
     * \param parameterValue New value of the arc-wise drag coefficient that is to be estimated.
     */
    void setParameterValue( Eigen::VectorXd parameterValue )
    {
        if( static_cast< int >( dragCoefficients_.size( ) ) !=
                static_cast< int >( parameterValue.rows( ) ) )
        {
            throw std::runtime_error( "Error when resetting arc-wise drag coefficients, sizes are incompatible" );
        }

        dragCoefficients_ = utilities::convertEigenVectorToStlVector( parameterValue );

        for( unsigned int i = 0; i < dragCoefficients_.size( ); i++ )
        {
            fullAerodynamicCoefficients_[ i ]( 0 ) = dragCoefficients_[ i ];
        }
        fullAerodynamicCoefficients_[ dragCoefficients_.size( ) ] =
                fullAerodynamicCoefficients_[ dragCoefficients_.size( ) - 1 ];
        coefficientInterpolator_->resetDependentValues( fullAerodynamicCoefficients_ );
    }

    //! Function to retrieve the size of the parameter
    /*!
     *  Function to retrieve the size of the parameter
     *  \return Size of parameter value
     */
    int getParameterSize( ){ return dragCoefficients_.size( ); }

    void updateCurrentCoefficients( const double time )
    {
        currentCoefficients_ = coefficientInterpolator_->interpolate( time );
    }

    Eigen::Vector6d getCurrentCoefficients( )
    {
        return currentCoefficients_;
    }

    std::shared_ptr< interpolators::LookUpScheme< double > > getArcTimeLookupScheme( )
    {
        return coefficientInterpolator_->getLookUpScheme( );
    }

    int getNumberOfArcs( )
    {
        return dragCoefficients_.size( );
    }

protected:

private:

    void setCoefficientInterfaceClosure( )
    {
        coefficientInterface_->setTimeDependentCoefficientClosure(
                    std::bind( &ArcWiseConstantDragCoefficient::getCurrentCoefficients, this ),
                    std::bind( &ArcWiseConstantDragCoefficient::updateCurrentCoefficients, this, std::placeholders::_1 ) );
    }

    //! Object containing the drag coefficient to be estimated.
    const std::shared_ptr< aerodynamics::CustomAerodynamicCoefficientInterface > coefficientInterface_;

    //! Times at which the arcs are to start (including end time at maximum double value).
    std::vector< double > timeLimits_;

    //! Values of drag coefficients in each arc.
    std::vector< double > dragCoefficients_;

    //! Values of drag coefficients in each arc, with additional value copied at end.
    std::vector< Eigen::Vector6d > fullAerodynamicCoefficients_;

    //! Interpolator that returns the drag coefficient as a function of time.
    std::shared_ptr< interpolators::PiecewiseConstantInterpolator< double, Eigen::Vector6d > > coefficientInterpolator_;

    Eigen::Vector6d currentCoefficients_;
};

} // namespace estimatable_parameters

} // namespace tudat


#endif // TUDAT_CONSTANTDRAGCOEFFICIENT_H
