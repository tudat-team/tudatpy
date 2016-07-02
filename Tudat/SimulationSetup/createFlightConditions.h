/*    Copyright (c) 2010-2016, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEFLIGHTCONDITIONS_H
#define TUDAT_CREATEFLIGHTCONDITIONS_H

#include <boost/multi_array.hpp>

#include <vector>

#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"
#include "Tudat/Astrodynamics/Aerodynamics/flightConditions.h"
#include "Tudat/Astrodynamics/Aerodynamics/customAerodynamicCoefficientInterface.h"
#include "Tudat/SimulationSetup/body.h"

namespace tudat
{

namespace simulation_setup
{


//! List of aerodynamic coefficient models available in simulations
/*!
 *  List of aerodynamic coefficient models available in simulations. Aerodynamic coefficient models
 *  not defined by this given enum cannot be used for automatic model setup.
 */
enum AerodynamicCoefficientTypes
{
    constant_aerodynamic_coefficients,
    hypersonic_local_inclincation_coefficients,
    tabulated_coefficients
};

//! Class for providing settings for aerodynamic coefficient model.
/*!
 *  Class for providing settings for automatic aerodynamic coefficient model creation. This class is
 *  a functional (base) class for settings of aerodynamic coefficient models that require no
 *  information in addition to their type. Aerodynamic coefficient model classes defining requiring
 *  additional information must be created using an object derived from this class.
 */
class AerodynamicCoefficientSettings
{
public:

    //! Constructor, sets type of aerodynamic coefficient model.
    /*!
     *  Constructor, sets type of aerodynamic coefficient model. Settings for aerodynamic
     *  coefficient models requiring additional information should be defined in a derived class.
     *  \param aerodynamicCoefficientTypes Type of aerodynamic coefficient model that is to be
     * created.
     *  \param referenceLength Reference length with which aerodynamic moments
     *  (about x- and z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
     *  is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
     *  coefficients are typically defined in negative direction.
     */
    AerodynamicCoefficientSettings(
            const AerodynamicCoefficientTypes aerodynamicCoefficientTypes,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = 1,
            const bool areCoefficientsInNegativeAxisDirection = 1 ):
        aerodynamicCoefficientTypes_( aerodynamicCoefficientTypes ),
        referenceLength_( referenceLength ), referenceArea_( referenceArea ),
        lateralReferenceLength_( lateralReferenceLength ),
        momentReferencePoint_( momentReferencePoint ),
        independentVariableNames_( independentVariableNames ),
        areCoefficientsInAerodynamicFrame_( areCoefficientsInAerodynamicFrame ),
        areCoefficientsInNegativeAxisDirection_( areCoefficientsInNegativeAxisDirection ){ }

    //! Destructor
    virtual ~AerodynamicCoefficientSettings( ){ }

    //! Function to return type of aerodynamic coefficient model that is to be created.
    /*!
     *  Function to return type of aerodynamic coefficient model that is to be created.
     *  \return Type of aerodynamic coefficient model that is to be created.
     */
    AerodynamicCoefficientTypes getAerodynamicCoefficientType( ){
        return aerodynamicCoefficientTypes_; }

    //! Get reference area.
    /*!
     * Returns reference area used to non-dimensionalize aerodynamic forces and moments.
     * \return Aerodynamic reference area.
     */
    double getReferenceArea( ) { return referenceArea_; }

    //! Get reference length.
    /*!
     * Returns reference length used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic reference length.
     */
    double getReferenceLength( ) { return referenceLength_; }

    //! Get lateral reference length.
    /*!
     * Returns lateral reference length used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic lateral reference length.
     */
    double getLateralReferenceLength( ) { return lateralReferenceLength_; }

    //! Get moment reference point.
    /*!
     * Returns the point w.r.t. which the arm of the aerodynamic moment on a vehicle panel is
     * determined.
     * \return Aerodynamic reference point.
     */
    Eigen::VectorXd getMomentReferencePoint( ) { return momentReferencePoint_; }

    //! Function to return identifiers of physical meaning of independent variables.
    /*!
     *  Function to return identifiers of physical meaning of independent variables.
     *  \return Identifiers of physical meaning of independent variables.
     */
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
    getIndependentVariableNames( )
    {
        return independentVariableNames_;
    }

    //! Function to return whether coefficients are in aerodynamic frame.
    /*!
     *  Function to return whether coefficients are in aerodynamic frame.
     *  \return Boolean defining whether coefficients are in aerodynamic frame.
     */
    bool getAreCoefficientsInAerodynamicFrame( )
    {
        return areCoefficientsInAerodynamicFrame_;
    }

    //! Function to return whether coefficients are positive along positive axes.
    /*!
     *  Function to return whether coefficients are positive along positive axes.
     *  \return Boolean defining whether coefficients are positive along positive axes.
     */
    bool getAreCoefficientsInNegativeAxisDirection( )
    {
        return areCoefficientsInNegativeAxisDirection_;
    }

private:

    //!  Type of atmosphere model that is to be created.
    AerodynamicCoefficientTypes aerodynamicCoefficientTypes_;

    //! Aerodynamic reference length.
    /*!
     * Reference length with which aerodynamic moments are non-dimensionalized.
     */
    double referenceLength_;

    //! Aerodynamic reference area.
    /*!
     * Reference area with which aerodynamic forces and moments are non-dimensionalized.
     */
    double referenceArea_;

    //! Lateral aerodynamic reference length.
    /*!
     * Lateral reference length with which aerodynamic moments are non-dimensionalized.
     */
    double lateralReferenceLength_;

    //! Aerodynamic moment reference point.
    /*!
     * Point w.r.t. which the arm of the moment on a vehicle panel is determined.
     */
    Eigen::Vector3d momentReferencePoint_;

    //! Vector with identifiers of the physical meaning of each independent variable of the
    //! aerodynamic coefficients.
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
    independentVariableNames_;

    //! Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame.
    /*!
     *  Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
     *  (lift, drag, side force) or in the body frame (typically denoted as Cx, Cy, Cz).
     */
    bool areCoefficientsInAerodynamicFrame_;

    //! Boolean to define whether the coefficients are positive along the positive axes.
    /*!
     *  Boolean to define whether the aerodynamic coefficients are positive along the positive
     *  axes of the body or aerodynamic frame (see areCoefficientsInAerodynamicFrame).
     *  Note that for (lift, drag, side force), the coefficients are typically defined in
     *  negative direction.
     */
    bool areCoefficientsInNegativeAxisDirection_;
};

//! AerodynamicCoefficientSettings for defining a constant aerodynamic coefficients
class ConstantAerodynamicCoefficientSettings: public AerodynamicCoefficientSettings
{
public:
    //! Constructor.
    /*!
     *  Constructor.
     *  \param constantForceCoefficient Constant force coefficients.
     *  \param constantMomentCoefficient Constant moment coefficients.
     *  \param referenceLength Reference length with which aerodynamic moments
     * (about x- and z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     * non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments
     * (about y-axis) is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
     *  coefficients are typically defined in negative direction.
     */
    ConstantAerodynamicCoefficientSettings(
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const Eigen::Vector3d& constantForceCoefficient,
            const Eigen::Vector3d& constantMomentCoefficient = Eigen::Vector3d::Zero( ),
            const bool areCoefficientsInAerodynamicFrame = 0,
            const bool areCoefficientsInNegativeAxisDirection = 1  ):
        AerodynamicCoefficientSettings(
            constant_aerodynamic_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint,
            std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
            areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ),
        constantForceCoefficient_( constantForceCoefficient ),
        constantMomentCoefficient_( constantMomentCoefficient )
    { }

   //! Constructor.
    /*!
    *  Constructor, omitting all moment coefficient data.
    *  \param constantForceCoefficient Constant force coefficients.
    *  \param referenceArea Reference area with which aerodynamic forces and moments are
    *  non-dimensionalized.
    *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
    *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
    *  frame (typically denoted as Cx, Cy, Cz).
    *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
    *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
    *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
    *  coefficients are typically defined in negative direction.
    */
    ConstantAerodynamicCoefficientSettings(
            const double referenceArea,
            const Eigen::Vector3d& constantForceCoefficient,
            const bool areCoefficientsInAerodynamicFrame = 0,
            const bool areCoefficientsInNegativeAxisDirection = 1 ):
        AerodynamicCoefficientSettings(
            constant_aerodynamic_coefficients, TUDAT_NAN, referenceArea,
            TUDAT_NAN, Eigen::Vector3d::Zero( ),
            std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
            areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ),
        constantForceCoefficient_( constantForceCoefficient ),
        constantMomentCoefficient_( Eigen::Vector3d::Zero( ) ){ }


    //! Function to return constant force coefficients.
    /*!
     *  Function to return constant force coefficients.
     *  \return Cnstant force coefficients.
     */
    Eigen::Vector3d getConstantForceCoefficient( )
    {
        return  constantForceCoefficient_;
    }

    //! Function to return constant moment coefficients.
    /*!
     *  Function to return constant moment coefficients.
     *  \return Cnstant force coefficients.
     */
    Eigen::Vector3d getConstantMomentCoefficient( )
    {
        return constantMomentCoefficient_;
    }

private:

    //! Constant moment coefficients.
    Eigen::Vector3d constantForceCoefficient_;

    //! Constant force coefficients.
    Eigen::Vector3d constantMomentCoefficient_;
};


//! Object for setting aerodynamic coefficients from a user-defined N-dimensional table.
/*!
 *  Object for setting aerodynamic coefficients from a user-defined N-dimensional table.
 *  The user must provide the force and moment coefficients in boost multi_arrays, and
 *  define the physical meaning of each of the independent variables.
 */
template< int NumberOfDimensions >
class TabulatedAerodynamicCoefficientSettings: public AerodynamicCoefficientSettings
{
public:

    //! Constructor, sets properties of aerodynamic coefficients.
    /*!
     *  Constructor, sets properties of aerodynamic coefficients.
     *  \param independentVariables Values of indepependent variables at which the coefficients
     *  in the input multi arrays are defined.
     *  \param forceCoefficients Values of force coefficients at independent variables defined
     *  by independentVariables.
     *  \param momentCoefficients Values of moment coefficients at independent variables defined
     *  by independentVariables.
     *  \param referenceLength Reference length with which aerodynamic moments
     *  (about x- and z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
     *  is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
     *  coefficients are typically defined in negative direction.
     */
    TabulatedAerodynamicCoefficientSettings(
            const std::vector< std::vector< double > > independentVariables,
            const boost::multi_array< Eigen::Vector3d, NumberOfDimensions > forceCoefficients,
            const boost::multi_array< Eigen::Vector3d, NumberOfDimensions > momentCoefficients,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = 1,
            const bool areCoefficientsInNegativeAxisDirection = 1 ):
        AerodynamicCoefficientSettings(
            tabulated_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint,
            independentVariableNames, areCoefficientsInAerodynamicFrame,
            areCoefficientsInNegativeAxisDirection ),
        independentVariables_( independentVariables ),
        forceCoefficients_( forceCoefficients ),
        momentCoefficients_( momentCoefficients ){ }

    //! Function to return the values of the indepependent variables of tables of coefficients.
    /*!
     *  Function to return the values of the indepependent variables of tables of coefficients.
     *  \return Values of the indepependent variables of tables of coefficients.
     */
    std::vector< std::vector< double > > getIndependentVariables( )
    {
        return independentVariables_;
    }

    //! Function to return values of force coefficients in table.
    /*!
     * Function to return values of force coefficients in table.
     * \return Values of force coefficients in table.
     */
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > getForceCoefficients( )
    {
        return forceCoefficients_;
    }

    //! Function to return values of moment coefficients in table.
    /*!
     * Function to return values of moment coefficients in table.
     * \return Values of moment coefficients in table.
     */
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > getMomentCoefficients( )
    {
        return momentCoefficients_;
    }

private:

    //! Values of indepependent variables at which the coefficients in the tables are defined.
    /*!
     *  Values of indepependent variables at which the coefficients in the forceCoefficients_ and
     *  momentCoefficients_ tables are defined.
     */
    std::vector< std::vector< double > > independentVariables_;

    //! Values of force coefficients at independent variables defined  by independentVariables_.
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > forceCoefficients_;

    //! Values of moment coefficients at independent variables defined  by independentVariables_.
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > momentCoefficients_;


};

//! Function to create an aerodynamic coefficient interface containing constant coefficients.
/*!
 *  Function to create an aerodynamic coefficient interface containing constant coefficients,
 *  As a result, the generated coefficient interface depends on zero parameters.
 *  \param constantForceCoefficient Constant force coefficients.
 *  \param constantMomentCoefficient Constant moment coefficients.
 *  \param referenceLength Reference length with which aerodynamic moments
 *  (about x- and z- axes) are non-dimensionalized.
 *  \param referenceArea Reference area with which aerodynamic forces and moments are
 *  non-dimensionalized.
 *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
 *  is non-dimensionalized.
 *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
 *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
 *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
 *  coefficients are typically defined in negative direction.
 *  \return Aerodynamic coefficient interface with constant coefficients.
 */
boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createConstantCoefficientAerodynamicCoefficientInterface(
        const Eigen::Vector3d constantForceCoefficient,
        const Eigen::Vector3d constantMomentCoefficient,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame = 0,
        const bool areCoefficientsInNegativeAxisDirection = 1 );

//! Factory function for tabulated aerodynamic coefficient interface.
/*!
 *  Factory function for tabulated aerodynamic coefficient interface.
 *  \param independentVariables Values of indepependent variables at which the coefficients
 *  in the input multi arrays are defined.
 *  \param forceCoefficients Values of force coefficients at independent variables defined
 *  by independentVariables.
 *  \param momentCoefficients Values of moment coefficients at independent variables defined
 *  by independentVariables.
 *  \param referenceLength Reference length with which aerodynamic moments
 *  (about x- and z- axes) are non-dimensionalized.
 *  \param referenceArea Reference area with which aerodynamic forces and moments are
 *  non-dimensionalized.
 *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
 *  is non-dimensionalized.
 *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
 *  \param independentVariableNames Vector with identifiers the physical meaning of each
 *  independent variable of the aerodynamic coefficients.
 *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
 *  coefficients are defined in the aerodynamic frame (lift, drag, side force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positive along the positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (lift, drag, side force), the
 *  coefficients are typically defined in negative direction.
 *  \return Tabulated aerodynamic coefficient interface pointer.
 */
template< int NumberOfDimensions >
boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createTabulatedCoefficientAerodynamicCoefficientInterface(
        const std::vector< std::vector< double > > independentVariables,
        const boost::multi_array< Eigen::Vector3d, NumberOfDimensions > forceCoefficients,
        const boost::multi_array< Eigen::Vector3d, NumberOfDimensions > momentCoefficients,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
        independentVariableNames,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame = 0,
        const bool areCoefficientsInNegativeAxisDirection = 1 )
{
    // Check input consistency.
    if( independentVariables.size( ) != NumberOfDimensions )
    {
        throw std::runtime_error( "Error when creating tabulated aerodynamic coefficient interface, inconsistent variable vector dimensioning" );
    }

    if( independentVariableNames.size( ) != NumberOfDimensions )
    {
       throw std::runtime_error( "Error when creating tabulated aerodynamic coefficient interface, inconsistent variable name vector dimensioning" );

    }

    // Create interpolators for coefficients.
    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > > forceInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > >(
                independentVariables, forceCoefficients );
    boost::shared_ptr< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > > momentInterpolator =
            boost::make_shared< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > >(
                independentVariables, momentCoefficients );

    // Create aerodynamic coefficient interface.
    return  boost::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                boost::bind( &interpolators::MultiLinearInterpolator
                             < double, Eigen::Vector3d, NumberOfDimensions >::interpolate,
                             forceInterpolator, _1 ),
                boost::bind( &interpolators::MultiLinearInterpolator
                             < double, Eigen::Vector3d, NumberOfDimensions >::interpolate,
                             momentInterpolator, _1 ),
                referenceLength, referenceArea, lateralReferenceLength, momentReferencePoint,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
}

//! Factory function for tabulated aerodynamic coefficient interface from coefficient settings.
/*!
 *  Factory function for tabulated aerodynamic coefficient interface from coefficient settings.
 *  This function is included to allow easier interface between the non-templated general
 *  createAerodynamicCoefficientInterface and the templated
 *  createTabulatedCoefficientAerodynamicCoefficientInterface.
 *  \param coefficientSettings Settings for aerodynamic coefficient interface, must be of derived
 *  type TabulatedAerodynamicCoefficientSettings< NumberOfDimensions >/
 *  \param body Name of body for which coefficient interface is to be made.
 *  \return Tabulated aerodynamic coefficient interface pointer.
 */
template< int NumberOfDimensions >
boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createTabulatedCoefficientAerodynamicCoefficientInterface(
        const boost::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    // Check consistency of type.
    boost::shared_ptr< TabulatedAerodynamicCoefficientSettings< NumberOfDimensions > > tabulatedCoefficientSettings =
            boost::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< NumberOfDimensions > >(
                coefficientSettings );
    if( tabulatedCoefficientSettings == NULL )
    {
        throw std::runtime_error(
                    "Error, expected tabulated aerodynamic coefficients of size " +
                    boost::lexical_cast<  std::string >( NumberOfDimensions ) + "for body " + body );
    }
    else
    {
        return createTabulatedCoefficientAerodynamicCoefficientInterface< NumberOfDimensions >(
                    tabulatedCoefficientSettings->getIndependentVariables( ),
                    tabulatedCoefficientSettings->getForceCoefficients( ),
                    tabulatedCoefficientSettings->getMomentCoefficients( ),
                    tabulatedCoefficientSettings->getIndependentVariableNames( ),
                    tabulatedCoefficientSettings->getReferenceLength( ),
                    tabulatedCoefficientSettings->getReferenceArea( ),
                    tabulatedCoefficientSettings->getReferenceLength( ),
                    tabulatedCoefficientSettings->getMomentReferencePoint( ),
                    tabulatedCoefficientSettings->getAreCoefficientsInAerodynamicFrame( ),
                    tabulatedCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( ) );
    }
}


//! Function to create an aerodynamic coefficient interface.
/*!
 * Function to create an aerodynamic coefficient interface from interface settings.
 * \param coefficientSettings Settings for the aerodynamic coefficient interface.
 * \param body Name of body for which aerodynamic coefficients are to be made.
 * \return Aerodynamic coefficient interface pointer of reqyested type and settings.
 */
boost::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createAerodynamicCoefficientInterface(
        const boost::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body );

//! Function to create a flight conditions object
/*!
 * Function to create a flight conditions object, which is responsible for calculating the various
 * dependent variables required for calculation of the aerodynamic acceleration
 * \param bodyWithFlightConditions Body for which flight conditions are to be created.
 * \param centralBody Body in  the atmosphere of which bodyWithFlightConditions is flying
 * \param nameOfBodyUndergoingAcceleration Name of body undergoing acceleration.
 * \param nameOfBodyExertingAcceleration Name of body with the atmosphere causing acceleration.
 * \param angleOfAttackFunction Function returning the current angle of attack (default 0).
 * \param angleOfSideslipFunction Function returning the current angle of sideslip (default 0).
 * \param bankAngleFunction Function returning the current bank angle (default 0).
 * \return Flight conditions object for given bodies and settings.
 */
boost::shared_ptr< aerodynamics::FlightConditions > createFlightConditions(
        const boost::shared_ptr< Body > bodyWithFlightConditions,
        const boost::shared_ptr< Body > centralBody,
        const std::string& nameOfBodyUndergoingAcceleration,
        const std::string& nameOfBodyExertingAcceleration,
        const boost::function< double( ) > angleOfAttackFunction =
        boost::lambda::constant ( 0.0 ),
        const boost::function< double( ) > angleOfSideslipFunction =
        boost::lambda::constant ( 0.0 ),
        const boost::function< double( ) > bankAngleFunction =
        boost::lambda::constant ( 0.0 ) );


} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATEACCELERATIONMODELS_H
