/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEAERODYNAMICCOEFFICIENTINTERFACE_H
#define TUDAT_CREATEAERODYNAMICCOEFFICIENTINTERFACE_H

#include <memory>
#include <boost/make_shared.hpp>

#include "tudat/astro/aerodynamics/aerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/customAerodynamicCoefficientInterface.h"
#include "tudat/simulation/environment_setup/createAerodynamicControlSurfaces.h"
#include "tudat/math/interpolators/multiLinearInterpolator.h"
#include "tudat/math/interpolators/createInterpolator.h"
#include "tudat/paths.hpp"

namespace tudat
{

namespace simulation_setup
{

//  Class for providing settings for aerodynamic coefficient model.
/*  
 *  Class for providing settings for automatic aerodynamic coefficient model creation. This class is
 *  a functional (base) class for settings of aerodynamic coefficient models that require no
 *  information in addition to their type. Aerodynamic coefficient model classes defining requiring
 *  additional information must be created using an object derived from this class.
 */

//! @get_docstring(AerodynamicCoefficientSettings.__docstring__)
class AerodynamicCoefficientSettings
{
public:

    //  Constructor, sets type of aerodynamic coefficient model.
    /*  
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
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction.
     *  \param interpolatorSettings Pointer to an interpolator settings object, where the
     *  conditions for interpolation are saved.
     */
    AerodynamicCoefficientSettings(
            const AerodynamicCoefficientTypes aerodynamicCoefficientTypes,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr ) :
        aerodynamicCoefficientTypes_( aerodynamicCoefficientTypes ),
        referenceLength_( referenceLength ), referenceArea_( referenceArea ),
        lateralReferenceLength_( lateralReferenceLength ),
        momentReferencePoint_( momentReferencePoint ),
        independentVariableNames_( independentVariableNames ),
        areCoefficientsInAerodynamicFrame_( areCoefficientsInAerodynamicFrame ),
        areCoefficientsInNegativeAxisDirection_( areCoefficientsInNegativeAxisDirection ),
        interpolatorSettings_( interpolatorSettings )
    { }

    //  Destructor
    virtual ~AerodynamicCoefficientSettings( ){ }

    //  Function to return type of aerodynamic coefficient model that is to be created.
    /*  
     *  Function to return type of aerodynamic coefficient model that is to be created.
     *  \return Type of aerodynamic coefficient model that is to be created.
     */
    AerodynamicCoefficientTypes getAerodynamicCoefficientType( ){
        return aerodynamicCoefficientTypes_; }

    //  Get reference area.
    /*  
     * Returns reference area used to non-dimensionalize aerodynamic forces and moments.
     * \return Aerodynamic reference area.
     */
    double getReferenceArea( ) { return referenceArea_; }

    //  Get reference length.
    /*  
     * Returns reference length used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic reference length.
     */
    double getReferenceLength( ) { return referenceLength_; }

    //  Get lateral reference length.
    /*  
     * Returns lateral reference length used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic lateral reference length.
     */
    double getLateralReferenceLength( ) { return lateralReferenceLength_; }

    //  Get moment reference point.
    /*  
     * Returns the point w.r.t. which the arm of the aerodynamic moment on a vehicle panel is
     * determined.
     * \return Aerodynamic reference point.
     */
    Eigen::VectorXd getMomentReferencePoint( ) { return momentReferencePoint_; }

    //  Function to return identifiers of physical meaning of independent variables.
    /*  
     *  Function to return identifiers of physical meaning of independent variables.
     *  \return Identifiers of physical meaning of independent variables.
     */
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
    getIndependentVariableNames( ) { return independentVariableNames_; }

    //  Function to return whether coefficients are in aerodynamic frame.
    /*  
     *  Function to return whether coefficients are in aerodynamic frame.
     *  \return Boolean defining whether coefficients are in aerodynamic frame.
     */
    bool getAreCoefficientsInAerodynamicFrame( ) { return areCoefficientsInAerodynamicFrame_; }

    //  Function to return whether coefficients are positive along positive axes.
    /*  
     *  Function to return whether coefficients are positive along positive axes.
     *  \return Boolean defining whether coefficients are positive along positive axes.
     */
    bool getAreCoefficientsInNegativeAxisDirection( ) { return areCoefficientsInNegativeAxisDirection_; }

    //  Function to return settings to be used for creating the interpoaltor of data.
    /*  
     *  Function to return settings to be used for creating the interpoaltor of data.
     *  \return Settings to be used for creating the one-dimensional interpoaltor of data.
     */
    std::shared_ptr< interpolators::InterpolatorSettings > getInterpolatorSettings( )
    {
        return interpolatorSettings_;
    }

    std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > >
    getControlSurfaceSettings( )
    {
        return controlSurfaceSettings_;
    }

    //  Function to define settings for the aerodynamic coefficients of a single control surface
    /*  
     * Function to define settings for the aerodynamic coefficients of a single control surface
     * \param controlSurfaceSetting Settings for the arodynamic coefficients of control surface.
     * \param controlSurfaceName Id of control surface.
     */
    void setControlSurfaceSettings(
            const std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > controlSurfaceSetting,
            const std::string controlSurfaceName )
    {
        controlSurfaceSettings_[ controlSurfaceName ] = controlSurfaceSetting;
    }

private:

    //  Type of atmosphere model that is to be created.
    AerodynamicCoefficientTypes aerodynamicCoefficientTypes_;

    //  Aerodynamic reference length.
    /*  
     * Reference length with which aerodynamic moments are non-dimensionalized.
     */
    double referenceLength_;

    //  Aerodynamic reference area.
    /*  
     * Reference area with which aerodynamic forces and moments are non-dimensionalized.
     */
    double referenceArea_;

    //  Lateral aerodynamic reference length.
    /*  
     * Lateral reference length with which aerodynamic moments are non-dimensionalized.
     */
    double lateralReferenceLength_;

    //  Aerodynamic moment reference point.
    /*  
     * Point w.r.t. which the arm of the moment on a vehicle panel is determined.
     */
    Eigen::Vector3d momentReferencePoint_;

    //  Vector with identifiers of the physical meaning of each independent variable of the
    //  aerodynamic coefficients.
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
    independentVariableNames_;

    //  Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame.
    /*  
     *  Boolean to define whether the aerodynamic coefficients are defined in the aerodynamic frame
     *  (drag, side, lift force) or in the body frame (typically denoted as Cx, Cy, Cz).
     */
    bool areCoefficientsInAerodynamicFrame_;

    //  Boolean to define whether the coefficients are positive along the positive axes.
    /*  
     *  Boolean to define whether the aerodynamic coefficients are positive along the positive
     *  axes of the body or aerodynamic frame (see areCoefficientsInAerodynamicFrame).
     *  Note that for (drag, side, lift force), the coefficients are typically defined in
     *  negative direction.
     */
    bool areCoefficientsInNegativeAxisDirection_;

    //  Settings for interpolation.
    /*  
     *  Settings for interpolation of aerodynamic coefficients, used to define an interpolator
     *  object, such that the coefficients are avaiable for a continuous set of independent variables.
     */
    std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings_;

    //  Settings for the aerodynamic coefficients of control surfaces, with map key denoting surface ID.
    std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > >
    controlSurfaceSettings_;
};


class ScaledAerodynamicCoefficientInterfaceSettings: public AerodynamicCoefficientSettings
{
public:
    ScaledAerodynamicCoefficientInterfaceSettings(
            const std::shared_ptr< AerodynamicCoefficientSettings > baseSettings,
            const Eigen::Vector3d forceScaling,
            const Eigen::Vector3d momentScaling,
            const bool isScalingAbsolute ):
        AerodynamicCoefficientSettings(
            scaled_coefficients, baseSettings->getReferenceLength( ), baseSettings->getReferenceArea( ),
            baseSettings->getLateralReferenceLength( ), baseSettings->getMomentReferencePoint( ), baseSettings->getIndependentVariableNames( ),
            baseSettings->getAreCoefficientsInAerodynamicFrame( ), baseSettings->getAreCoefficientsInNegativeAxisDirection(),
            baseSettings->getInterpolatorSettings( ) ),
        baseSettings_( baseSettings ),
        forceScaling_( [=]( const double ){ return forceScaling; } ),
        momentScaling_( [=]( const double ){ return momentScaling; } ),
        isScalingAbsolute_( isScalingAbsolute ){ }

    ScaledAerodynamicCoefficientInterfaceSettings(
            const std::shared_ptr< AerodynamicCoefficientSettings > baseSettings,
            const std::function< Eigen::Vector3d( const double ) > forceScaling,
            const std::function< Eigen::Vector3d( const double ) > momentScaling,
            const bool isScalingAbsolute ):
        AerodynamicCoefficientSettings(
            scaled_coefficients, baseSettings->getReferenceLength( ), baseSettings->getReferenceArea( ),
            baseSettings->getLateralReferenceLength( ), baseSettings->getMomentReferencePoint( ), baseSettings->getIndependentVariableNames( ),
            baseSettings->getAreCoefficientsInAerodynamicFrame( ), baseSettings->getAreCoefficientsInNegativeAxisDirection(),
            baseSettings->getInterpolatorSettings( ) ),
        baseSettings_( baseSettings ), forceScaling_( forceScaling ), momentScaling_( momentScaling ),
        isScalingAbsolute_( isScalingAbsolute ){ }

    std::shared_ptr< AerodynamicCoefficientSettings > getBaseSettings( )
    {
        return baseSettings_;
    }

    std::function< Eigen::Vector3d( const double ) > getForceScaling( )
    {
        return forceScaling_;
    }

    std::function< Eigen::Vector3d( const double ) > getMomentScaling( )
    {
        return momentScaling_;
    }

    bool getIsScalingAbsolute( )
    {
        return isScalingAbsolute_;
    }

protected:

    std::shared_ptr< AerodynamicCoefficientSettings > baseSettings_;

    const std::function< Eigen::Vector3d( const double ) > forceScaling_;

    const std::function< Eigen::Vector3d( const double ) > momentScaling_;

    bool isScalingAbsolute_;
};

//  AerodynamicCoefficientSettings for defining a constant aerodynamic coefficients
//! @get_docstring(ConstantAerodynamicCoefficientSettings.__docstring__)
class ConstantAerodynamicCoefficientSettings: public AerodynamicCoefficientSettings
{
public:

    //  Constructor.
    /*  
     *  Constructor.
     *  \param constantForceCoefficient Constant force coefficients.
     *  \param constantMomentCoefficient Constant moment coefficients.
     *  \param referenceLength Reference length with which aerodynamic moments
     *  (about x- and z- axes) are non-dimensionalized.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param lateralReferenceLength Reference length with which aerodynamic moments
     *  (about y-axis) is non-dimensionalized.
     *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction.
     *  \param interpolatorSettings Pointer to an interpolator settings object, where the
     *  conditions for interpolation are saved.
     */
    ConstantAerodynamicCoefficientSettings(
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const Eigen::Vector3d& constantForceCoefficient,
            const Eigen::Vector3d& constantMomentCoefficient = Eigen::Vector3d::Zero( ),
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr ) :
        AerodynamicCoefficientSettings(
            constant_aerodynamic_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint,
            std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
            areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection, interpolatorSettings ),
        constantForceCoefficient_( constantForceCoefficient ),
        constantMomentCoefficient_( constantMomentCoefficient )
    { }

    //  Constructor.
    /*  
    *  Constructor, omitting all moment coefficient data.
    *  \param constantForceCoefficient Constant force coefficients.
    *  \param referenceArea Reference area with which aerodynamic forces and moments are
    *  non-dimensionalized.
    *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
    *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
    *  frame (typically denoted as Cx, Cy, Cz).
    *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
    *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
    *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
    *  coefficients are typically defined in negative direction.
    */
    ConstantAerodynamicCoefficientSettings(
            const double referenceArea,
            const Eigen::Vector3d& constantForceCoefficient,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true ):
        AerodynamicCoefficientSettings(
            constant_aerodynamic_coefficients, TUDAT_NAN, referenceArea,
            TUDAT_NAN, Eigen::Vector3d::Zero( ),
            std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >( ),
            areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection, nullptr ),
        constantForceCoefficient_( constantForceCoefficient ),
        constantMomentCoefficient_( Eigen::Vector3d::Zero( ) ){ }

    //  Function to return constant force coefficients.
    /*  
     *  Function to return constant force coefficients.
     *  \return Cnstant force coefficients.
     */
    Eigen::Vector3d getConstantForceCoefficient( )
    {
        return  constantForceCoefficient_;
    }

    //  Function to return constant moment coefficients.
    /*  
     *  Function to return constant moment coefficients.
     *  \return Cnstant force coefficients.
     */
    Eigen::Vector3d getConstantMomentCoefficient( )
    {
        return constantMomentCoefficient_;
    }

private:

    //  Constant moment coefficients.
    Eigen::Vector3d constantForceCoefficient_;

    //  Constant force coefficients.
    Eigen::Vector3d constantMomentCoefficient_;

};

//  AerodynamicCoefficientSettings for defining a constant aerodynamic coefficients
class CustomAerodynamicCoefficientSettings: public AerodynamicCoefficientSettings
{
public:

    CustomAerodynamicCoefficientSettings(
            const std::function< Eigen::Vector3d( const std::vector< double >& ) > forceCoefficientFunction,
            const std::function< Eigen::Vector3d( const std::vector< double >& ) > momentCoefficientFunction,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true ) :
        AerodynamicCoefficientSettings(
            custom_aerodynamic_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint,
            independentVariableNames,
            areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection ),
        forceCoefficientFunction_( forceCoefficientFunction ),
        momentCoefficientFunction_( momentCoefficientFunction )
    { }

    std::function< Eigen::Vector3d( const std::vector< double >& ) > getForceCoefficientFunction( )
    {
        return forceCoefficientFunction_;
    }

    std::function< Eigen::Vector3d( const std::vector< double >& ) > getMomentCoefficientFunction( )
    {
        return momentCoefficientFunction_;
    }


private:
    std::function< Eigen::Vector3d( const std::vector< double >& ) > forceCoefficientFunction_;

    std::function< Eigen::Vector3d( const std::vector< double >& ) > momentCoefficientFunction_;


};

//! @get_docstring(constantAerodynamicCoefficientSettings)
inline std::shared_ptr< AerodynamicCoefficientSettings > constantAerodynamicCoefficientSettings(
        const double referenceArea,
        const Eigen::Vector3d& constantForceCoefficient,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true )
{
    return std::make_shared< ConstantAerodynamicCoefficientSettings >(
                referenceArea, constantForceCoefficient,
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
}

//! @get_docstring(scaledAerodynamicCoefficientSettings)
inline std::shared_ptr< AerodynamicCoefficientSettings > scaledAerodynamicCoefficientSettings(
        const std::shared_ptr< AerodynamicCoefficientSettings > baseSettings,
        const double forceScaling,
        const double momentScaling,
        const bool isScalingAbsolute )
{
    std::function< Eigen::Vector3d( const double ) > forceScalingFunction =
            [=]( const double ){ return Eigen::Vector3d::Constant( forceScaling ); };
    std::function< Eigen::Vector3d( const double ) > momentScalingFunction =
            [=]( const double ){ return Eigen::Vector3d::Constant( momentScaling ); };
    return std::make_shared< ScaledAerodynamicCoefficientInterfaceSettings >(
                baseSettings, forceScalingFunction, momentScalingFunction, isScalingAbsolute );
}

//! @get_docstring(scaledAerodynamicCoefficientSettings, 1)
inline std::shared_ptr< AerodynamicCoefficientSettings > scaledAerodynamicCoefficientSettings(
        const std::shared_ptr< AerodynamicCoefficientSettings > baseSettings,
        const Eigen::Vector3d forceScaling,
        const Eigen::Vector3d momentScaling,
        const bool isScalingAbsolute )
{
    std::function< Eigen::Vector3d( const double ) > forceScalingFunction =
            [=]( const double ){ return forceScaling; };
    std::function< Eigen::Vector3d( const double ) > momentScalingFunction =
            [=]( const double ){ return momentScaling; };
    return std::make_shared< ScaledAerodynamicCoefficientInterfaceSettings >(
                baseSettings, forceScalingFunction, momentScalingFunction, isScalingAbsolute );
}

//! @get_docstring(scaledAerodynamicCoefficientSettings, 2)
inline std::shared_ptr< AerodynamicCoefficientSettings > scaledAerodynamicCoefficientSettings(
        const std::shared_ptr< AerodynamicCoefficientSettings > baseSettings,
        const std::function< Eigen::Vector3d( const double ) > forceScaling,
        const std::function< Eigen::Vector3d( const double ) > momentScaling,
        const bool isScalingAbsolute )
{
    return std::make_shared< ScaledAerodynamicCoefficientInterfaceSettings >(
                baseSettings, forceScaling, momentScaling, isScalingAbsolute );
}

inline std::shared_ptr< AerodynamicCoefficientSettings > customAerodynamicCoefficientSettings(
        const std::function< Eigen::Vector3d( const std::vector< double >& ) > forceCoefficientFunction,
        const std::function< Eigen::Vector3d( const std::vector< double >& ) > momentCoefficientFunction,
        const double referenceLength,
        const double referenceArea,
        const Eigen::Vector3d& momentReferencePoint,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
        independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true )
{
    return std::make_shared< CustomAerodynamicCoefficientSettings >(
                forceCoefficientFunction, momentCoefficientFunction, referenceLength, referenceArea, referenceLength,
                momentReferencePoint, independentVariableNames, areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
}

//! @get_docstring(customAerodynamicCoefficientSettings)
inline std::shared_ptr< AerodynamicCoefficientSettings > customAerodynamicCoefficientSettings(
        const std::function< Eigen::Vector3d( const std::vector< double >& ) > forceCoefficientFunction,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
        independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true )
{
    return std::make_shared< CustomAerodynamicCoefficientSettings >(
                forceCoefficientFunction, [=](const std::vector< double >& ){ return Eigen::Vector3d::Constant( TUDAT_NAN ); },
    TUDAT_NAN, referenceArea, TUDAT_NAN, Eigen::Vector3d::Constant( TUDAT_NAN ),
    independentVariableNames, areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
}


//  Base class (non-functional) for the different classes of TabulatedAerodynamicCoefficientSettings.
/*  
 * Base class (non-functional) for the different classes of TabulatedAerodynamicCoefficientSettings.
 */
class TabulatedAerodynamicCoefficientSettingsBase: public AerodynamicCoefficientSettings
{
public:

    // Inherit constructors
    using AerodynamicCoefficientSettings::AerodynamicCoefficientSettings;

    //  Function to return files for force coefficients.
    /*  
     * Function to return files for force coefficients.
     * \return Files for force coefficients.
     */
    std::map< int, std::string > getForceCoefficientsFiles( )
    {
        return forceCoefficientsFiles_;
    }

    //  Function to return files for moment coefficients.
    /*  
     * Function to return files for moment coefficients.
     * \return Files for moment coefficients.
     */
    std::map< int, std::string > getMomentCoefficientsFiles( )
    {
        return momentCoefficientsFiles_;
    }

    //  Function to set the force coefficients files.
    /*  
     * Function to set the force coefficients files.
     * \param forceCoefficientsFiles The force coefficients files.
     */
    void setForceCoefficientsFiles( const std::map< int, std::string >& forceCoefficientsFiles )
    {
        forceCoefficientsFiles_ = forceCoefficientsFiles;
    }

    //  Function to set the moment coefficients files.
    /*  
     * Function to set the moment coefficients files.
     * \param momentCoefficientsFiles The moment coefficients files.
     */
    void setMomentCoefficientsFiles( const std::map< int, std::string >& momentCoefficientsFiles )
    {
        momentCoefficientsFiles_ = momentCoefficientsFiles;
    }

private:

    //  Files from which the force coefficients should be loaded.
    std::map< int, std::string > forceCoefficientsFiles_;

    //  Files from which the moment coefficients should be loaded.
    std::map< int, std::string > momentCoefficientsFiles_;

};

//  Object for setting aerodynamic coefficients from a user-defined N-dimensional table (with N>1).
/*  
 *  Object for setting aerodynamic coefficients from a user-defined N-dimensional table (with N>1). The N=1 case has its
 *  own template specialization.
 *  The user must provide the force (and moment) coefficients in boost multi_arrays, and
 *  define the physical meaning of each of the independent variables.
 */
template< unsigned int NumberOfDimensions >
class TabulatedAerodynamicCoefficientSettings: public TabulatedAerodynamicCoefficientSettingsBase
{
public:

    //  Constructor, sets properties of aerodynamic coefficients.
    /*  
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
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction.
     *  \param interpolatorSettings Pointer to an interpolator settings object, where the
     *  conditions for interpolation are saved.
     */
    TabulatedAerodynamicCoefficientSettings(
            const std::vector< std::vector< double > > independentVariables,
            const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > forceCoefficients,
            const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > momentCoefficients,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr ) :
        TabulatedAerodynamicCoefficientSettingsBase(
            tabulated_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint,
            independentVariableNames, areCoefficientsInAerodynamicFrame,
            areCoefficientsInNegativeAxisDirection, interpolatorSettings ),
        independentVariables_( independentVariables ),
        forceCoefficients_( forceCoefficients ),
        momentCoefficients_( momentCoefficients )
    { }

    //  Constructor, sets properties of aerodynamic force coefficients, zero moment coefficients.
    /*  
     *  Constructor, sets properties of aerodynamic force coefficients, zero moment coefficients
     *  \param independentVariables Values of indepependent variables at which the coefficients
     *  in the input multi arrays are defined.
     *  \param forceCoefficients Values of force coefficients at independent variables defined
     *  by independentVariables.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction.
     *  \param interpolatorSettings Pointer to an interpolator settings object, where the
     *  conditions for interpolation are saved.
     */
    TabulatedAerodynamicCoefficientSettings(
            const std::vector< std::vector< double > > independentVariables,
            const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > forceCoefficients,
            const double referenceArea,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr ) :
        TabulatedAerodynamicCoefficientSettingsBase(
            tabulated_coefficients, TUDAT_NAN, referenceArea,
            TUDAT_NAN, Eigen::Vector3d::Constant( TUDAT_NAN ),
            independentVariableNames, areCoefficientsInAerodynamicFrame,
            areCoefficientsInNegativeAxisDirection, interpolatorSettings ),
        independentVariables_( independentVariables ),
        forceCoefficients_( forceCoefficients )
    {
        std::vector< size_t > sizeVector;
        const size_t* arrayShape = forceCoefficients_.shape( );
        sizeVector.assign( arrayShape, arrayShape+ forceCoefficients_.num_dimensions( ) );

        momentCoefficients_.resize( sizeVector );

        std::fill( momentCoefficients_.data( ),
                   momentCoefficients_.data( ) + momentCoefficients_.num_elements( ), Eigen::Vector3d::Zero( ) );
    }

    //  Destructor
    ~TabulatedAerodynamicCoefficientSettings( ){ }

    //  Function to return the values of the indepependent variables of tables of coefficients.
    /*  
     *  Function to return the values of the indepependent variables of tables of coefficients.
     *  \return Values of the indepependent variables of tables of coefficients.
     */
    std::vector< std::vector< double > > getIndependentVariables( )
    {
        return independentVariables_;
    }

    //  Function to return values of force coefficients in table.
    /*  
     * Function to return values of force coefficients in table.
     * \return Values of force coefficients in table.
     */
    boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > getForceCoefficients( )
    {
        return forceCoefficients_;
    }

    //  Function to return values of moment coefficients in table.
    /*  
     * Function to return values of moment coefficients in table.
     * \return Values of moment coefficients in table.
     */
    boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > getMomentCoefficients( )
    {
        return momentCoefficients_;
    }

private:

    //  Values of indepependent variables at which the coefficients in the tables are defined.
    /*  
     *  Values of indepependent variables at which the coefficients in the forceCoefficients_ and
     *  momentCoefficients_ tables are defined.
     */
    std::vector< std::vector< double > > independentVariables_;

    //  Values of force coefficients at independent variables defined  by independentVariables_.
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > forceCoefficients_;

    //  Values of moment coefficients at independent variables defined  by independentVariables_.
    boost::multi_array< Eigen::Vector3d, NumberOfDimensions > momentCoefficients_;
};

//  Object for setting aerodynamic coefficients from a user-defined 1-dimensional table.
/*  
 *  Object for setting aerodynamic coefficients from a user-defined 1-dimensional table.
 *  The user must provide the force (and moment) coefficients in std::vectors, and
 *  define the physical meaning of the independent variables.
 */
template< >
class TabulatedAerodynamicCoefficientSettings< 1 >: public TabulatedAerodynamicCoefficientSettingsBase
{
public:

    //  Constructor, sets properties of aerodynamic coefficients.
    /*  
     *  Constructor, sets properties of aerodynamic coefficients.
     *  \param independentVariables Values of indepependent variables at which the coefficients
     *  in the input multi vector are defined.
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
     *  \param independentVariableName Identifiers the of physical meaning of the
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction.
     *  \param interpolatorSettings Pointer to an interpolator settings object, where the
     *  conditions for interpolation are saved.
     */
    TabulatedAerodynamicCoefficientSettings< 1 >(
            const std::vector< double > independentVariables,
            const std::vector< Eigen::Vector3d > forceCoefficients,
            const std::vector< Eigen::Vector3d > momentCoefficients,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const aerodynamics::AerodynamicCoefficientsIndependentVariables independentVariableName,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr ) :
        TabulatedAerodynamicCoefficientSettingsBase(
            tabulated_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint, { independentVariableName }, areCoefficientsInAerodynamicFrame,
            areCoefficientsInNegativeAxisDirection, interpolatorSettings )
    {
        if( forceCoefficients.size( ) != independentVariables.size( ) )
        {
            throw std::runtime_error( "Error, force coefficient size is inconsistent in TabulatedAerodynamicCoefficientSettings< 1 >" );
        }

        if( momentCoefficients.size( ) != independentVariables.size( ) )
        {
            throw std::runtime_error( "Error, moment coefficient size is inconsistent in TabulatedAerodynamicCoefficientSettings< 1 >" );
        }

        for( unsigned int i = 0; i < independentVariables.size( ); i++ )
        {
            forceCoefficients_[ independentVariables.at( i ) ] = forceCoefficients.at( i );
            momentCoefficients_[ independentVariables.at( i ) ] = momentCoefficients.at( i );
        }
    }

    //  Constructor, sets properties of aerodynamic coefficients.
    /*  
     *  Constructor, sets properties of aerodynamic coefficients.
     *  \param independentVariables Values of indepependent variables at which the coefficients
     *  in the input multi vector are defined.
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
     *  \param independentVariableName Identifiers the of physical meaning of the
     *  independent variable of the aerodynamic coefficients  (size 1).
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction.
     *  \param interpolatorSettings Pointer to an interpolator settings object, where the
     *  conditions for interpolation are saved.
     */
    TabulatedAerodynamicCoefficientSettings< 1 >(
            const std::vector< std::vector< double > > independentVariables,
            const boost::multi_array< Eigen::Vector3d, 1 > forceCoefficients,
            const boost::multi_array< Eigen::Vector3d, 1 > momentCoefficients,
            const double referenceLength,
            const double referenceArea,
            const double lateralReferenceLength,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableName,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr ) :
        TabulatedAerodynamicCoefficientSettingsBase(
            tabulated_coefficients, referenceLength, referenceArea,
            lateralReferenceLength, momentReferencePoint,
            independentVariableName, areCoefficientsInAerodynamicFrame,
            areCoefficientsInNegativeAxisDirection, interpolatorSettings )
    {
        if( forceCoefficients.size( ) != independentVariables.size( ) )
        {
            throw std::runtime_error( "Error, force coefficient size is inconsistent in TabulatedAerodynamicCoefficientSettings< 1 >" );
        }

        if( momentCoefficients.size( ) != independentVariables.size( ) )
        {
            throw std::runtime_error( "Error, moment coefficient size is inconsistent in TabulatedAerodynamicCoefficientSettings< 1 >" );
        }

        for( unsigned int i = 0; i < independentVariables.size( ); i++ )
        {
            forceCoefficients_[ independentVariables.at( 0 ).at( i ) ] = forceCoefficients[ i ];
            momentCoefficients_[ independentVariables.at( 0 ).at( i ) ] = momentCoefficients[ i ];
        }
    }

    //  Constructor, sets properties of aerodynamic force coefficients, zero moment coefficients.
    /*  
     *  Constructor, sets properties of aerodynamic force coefficients, zero moment coefficients.
     *  \param independentVariables Values of indepependent variables at which the coefficients
     *  in the input multi vector are defined.
     *  \param forceCoefficients Values of force coefficients at independent variables defined
     *  by independentVariables.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param independentVariableName Identifiers the of physical meaning of the
     *  independent variable of the aerodynamic coefficients.
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction.
     *  \param interpolatorSettings Pointer to an interpolator settings object, where the
     *  conditions for interpolation are saved.
     */
    TabulatedAerodynamicCoefficientSettings< 1 >(
            const std::vector< double > independentVariables,
            const std::vector< Eigen::Vector3d > forceCoefficients,
            const double referenceArea,
            const aerodynamics::AerodynamicCoefficientsIndependentVariables independentVariableName,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr ) :
        TabulatedAerodynamicCoefficientSettingsBase(
            tabulated_coefficients, TUDAT_NAN, referenceArea,
            TUDAT_NAN, Eigen::Vector3d::Constant( TUDAT_NAN ), { independentVariableName }, areCoefficientsInAerodynamicFrame,
            areCoefficientsInNegativeAxisDirection, interpolatorSettings )
    {
        if( forceCoefficients.size( ) != independentVariables.size( ) )
        {
            throw std::runtime_error( "Error, force coefficient size is inconsistent in TabulatedAerodynamicCoefficientSettings< 1 >" );
        }

        for( unsigned int i = 0; i < independentVariables.size( ); i++ )
        {
            forceCoefficients_[ independentVariables.at( i ) ] = forceCoefficients.at( i );
            momentCoefficients_[ independentVariables.at( i ) ] = Eigen::Vector3d::Zero( );
        }
    }

    //  Constructor, sets properties of aerodynamic force coefficients, zero moment coefficients.
    /*  
     *  Constructor, sets properties of aerodynamic force coefficients, zero moment coefficients.
     *  \param independentVariables Values of indepependent variables at which the coefficients
     *  in the input multi vector are defined (size 1).
     *  \param forceCoefficients Values of force coefficients at independent variables defined
     *  by independentVariables.
     *  \param referenceArea Reference area with which aerodynamic forces and moments are
     *  non-dimensionalized.
     *  \param independentVariableNames Identifiers the of physical meaning of the
     *  independent variable of the aerodynamic coefficients (size 1).
     *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
     *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
     *  frame (typically denoted as Cx, Cy, Cz).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along the positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction.
     *  \param interpolatorSettings Pointer to an interpolator settings object, where the
     *  conditions for interpolation are saved.
     */
    TabulatedAerodynamicCoefficientSettings< 1 >(
            const std::vector< std::vector< double > > independentVariables,
            const boost::multi_array< Eigen::Vector3d, 1 > forceCoefficients,
            const double referenceArea,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
            const bool areCoefficientsInAerodynamicFrame = true,
            const bool areCoefficientsInNegativeAxisDirection = true,
            const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr ) :
        TabulatedAerodynamicCoefficientSettingsBase(
            tabulated_coefficients, TUDAT_NAN, referenceArea,
            TUDAT_NAN, Eigen::Vector3d::Constant( TUDAT_NAN ),
            independentVariableNames, areCoefficientsInAerodynamicFrame,
            areCoefficientsInNegativeAxisDirection, interpolatorSettings )
    {
        if( forceCoefficients.shape( )[ 0 ] != independentVariables.size( ) )
        {
            throw std::runtime_error( "Error, force coefficient size is inconsistent in TabulatedAerodynamicCoefficientSettings< 1 >" );
        }

        for( unsigned int i = 0; i < independentVariables.size( ); i++ )
        {
            forceCoefficients_[ independentVariables.at( 0 ).at( i ) ] = forceCoefficients[ i ];
            momentCoefficients_[ independentVariables.at( 0 ).at( i ) ] = Eigen::Vector3d::Zero( );
        }
    }

    //  Destructor
    ~TabulatedAerodynamicCoefficientSettings< 1 >( ){ }

    //  Function to return values of force coefficients in table.
    /*  
     * Function to return values of force coefficients in table.
     * \return Values of force coefficients in table.
     */
    std::map< double, Eigen::Vector3d >  getForceCoefficients( )
    {
        return forceCoefficients_;
    }

    //  Function to return values of moment coefficients in table.
    /*  
     * Function to return values of moment coefficients in table.
     * \return Values of moment coefficients in table.
     */
    std::map< double, Eigen::Vector3d >  getMomentCoefficients( )
    {
        return momentCoefficients_;
    }

    //  Function to return settings to be used for creating the one-dimensional interpoaltor of data.
    /*  
     * Function to return settings to be used for creating the one-dimensional interpoaltor of data.
     * \return Settings to be used for creating the one-dimensional interpoaltor of data.
     */
    std::shared_ptr< interpolators::InterpolatorSettings > getInterpolationSettings( )
    {
        return interpolationSettings_;
    }

private:

    //  Values of force coefficients at independent variables defined  by independentVariables_.
    std::map< double, Eigen::Vector3d > forceCoefficients_;

    //  Values of moment coefficients at independent variables defined  by independentVariables_.
    std::map< double, Eigen::Vector3d > momentCoefficients_;

    //  Settings to be used for creating the one-dimensional interpolator of data.
    std::shared_ptr< interpolators::InterpolatorSettings > interpolationSettings_;

};

// 1-dimensional case
//! @get_docstring(oneDimensionalTabulatedAerodynamicCoefficientSettings)
inline std::shared_ptr< AerodynamicCoefficientSettings > oneDimensionalTabulatedAerodynamicCoefficientSettings(
        const std::vector< double > independentVariables,
        const std::vector< Eigen::Vector3d > forceCoefficients,
        const std::vector< Eigen::Vector3d > momentCoefficients,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const aerodynamics::AerodynamicCoefficientsIndependentVariables independentVariableName,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr
        )
{
    return std::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                independentVariables, forceCoefficients, momentCoefficients, referenceLength,
                referenceArea, lateralReferenceLength, momentReferencePoint, independentVariableName,
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection,
                interpolatorSettings );
}

//! @get_docstring(oneDimensionalTabulatedAerodynamicCoefficientSettings, 1)
inline std::shared_ptr< AerodynamicCoefficientSettings > oneDimensionalTabulatedAerodynamicCoefficientSettings(
        const std::vector< double > independentVariables,
        const std::vector< Eigen::Vector3d > forceCoefficients,
        const double referenceArea,
        const aerodynamics::AerodynamicCoefficientsIndependentVariables independentVariableName,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr
        )
{
    return std::make_shared< TabulatedAerodynamicCoefficientSettings< 1 > >(
                independentVariables, forceCoefficients, referenceArea, independentVariableName,
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection,
                interpolatorSettings );
}


//  Function to create aerodynamic coefficient settings from coefficients stored in data files
/*  
 *  Function to create aerodynamic coefficient settings from coefficients stored in data files. Separate files are defined for
 *  the three components of the force coefficients.  The file format is discussed in the Tudat wiki
 *  Note that this function requires the number of independent variables in the coefficient files to be known. If this is not
 *  the case, the readTabulatedAerodynamicCoefficientsFromFiles function should be used.
 *  \param forceCoefficientFiles List (size 3) of files containing the aerodynamic force coefficients
 *  \param momentCoefficientFiles List (size 3) of files containing the aerodynamic moment coefficients
 *  \param referenceLength Reference length with which aerodynamic moments
 *  (about x- and z- axes) are non-dimensionalized.
 *  \param referenceArea Reference area with which aerodynamic forces and moments are
 *  non-dimensionalized.
 *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
 *  is non-dimensionalized.
 *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
 *  \param independentVariableNames Physical meaning of the independent variables of the aerodynamic coefficients
 *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
 *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positive along the positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
 *  coefficients are typically defined in negative direction.
 *  \param interpolatorSettings Pointer to an interpolator settings object, where the
 *  conditions for interpolation are saved.
 *  \return Settings for creation of aerodynamic coefficient interface, based on contents read from files defined in
 *  forceCoefficientFiles and reference data given as input.
 */
template< unsigned int NumberOfDimensions >
std::shared_ptr< AerodynamicCoefficientSettings >
readGivenSizeTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr )
{
    std::pair< boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) >,
            std::vector< std::vector< double > > >
            aerodynamicForceCoefficients =
            input_output::readAerodynamicCoefficients< NumberOfDimensions >( forceCoefficientFiles );
    std::pair< boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) >,
            std::vector< std::vector< double > > >
            aerodynamicMomentCoefficients =
            input_output::readAerodynamicCoefficients< NumberOfDimensions >( momentCoefficientFiles );

    if( !input_output::compareIndependentVariables(
                aerodynamicForceCoefficients.second, aerodynamicMomentCoefficients.second ) )
    {
        throw std::runtime_error( "Error when creating aerodynamic coefficient settings from file, "
                                  "force and moment independent variables are inconsistent" );
    }

    if( independentVariableNames.size( ) != NumberOfDimensions )
    {
        throw std::runtime_error( "Error when creating aerodynamic coefficient settings from file, input sizes are inconsistent" );
    }

    // Create coefficient settings.
    std::shared_ptr< TabulatedAerodynamicCoefficientSettings< NumberOfDimensions > > tabulatedCoefficients =
            std::make_shared< TabulatedAerodynamicCoefficientSettings< NumberOfDimensions > >(
                aerodynamicForceCoefficients.second, aerodynamicForceCoefficients.first, aerodynamicMomentCoefficients.first,
                referenceLength, referenceArea, lateralReferenceLength, momentReferencePoint, independentVariableNames,
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection, interpolatorSettings );
    tabulatedCoefficients->setForceCoefficientsFiles( forceCoefficientFiles );
    tabulatedCoefficients->setMomentCoefficientsFiles( momentCoefficientFiles );
    return tabulatedCoefficients;
}

//  Function to create aerodynamic coefficient settings from coefficients stored in data files
/*  
 *  Function to create aerodynamic coefficient settings from coefficients stored in data files. Separate files are defined for
 *  the three components of the force coefficients. From this function, no moment coefficients are read (set to zero for all
 *  cases). The file format is discussed in the Tudat wiki
 *  Note that this function requires the number of independent variables in the coefficient files to be known. If this is not
 *  the case, the readTabulatedAerodynamicCoefficientsFromFiles function should be used.
 *  \param forceCoefficientFiles List (size 3) of files containing the aerodynamic coefficients
 *  \param referenceArea Reference area of aerodynamic coefficients
 *  \param independentVariableNames Physical meaning of the independent variables of the aerodynamic coefficients
 *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
 *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positive along the positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
 *  coefficients are typically defined in negative direction.
 *  \param interpolatorSettings Pointer to an interpolator settings object, where the
 *  conditions for interpolation are saved.
 *  \return Settings for creation of aerodynamic coefficient interface, based on contents read from files defined in
 *  forceCoefficientFiles and reference data given as input.
 */
template< unsigned int NumberOfDimensions >
std::shared_ptr< AerodynamicCoefficientSettings >
readGivenSizeTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr )
{
    std::pair< boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) >,
            std::vector< std::vector< double > > >
            aerodynamicCoefficients = input_output::readAerodynamicCoefficients< NumberOfDimensions >( forceCoefficientFiles );

    // Check input consistency
    if( independentVariableNames.size( ) != NumberOfDimensions )
    {
        throw std::runtime_error( "Error when creating aerodynamic coefficient settings from file, input sizes are inconsistent" );
    }

    // Create coefficient settings.
    std::shared_ptr< TabulatedAerodynamicCoefficientSettings< NumberOfDimensions > > tabulatedCoefficients =
            std::make_shared< TabulatedAerodynamicCoefficientSettings< NumberOfDimensions > >(
                aerodynamicCoefficients.second, aerodynamicCoefficients.first, referenceArea, independentVariableNames,
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection, interpolatorSettings );
    tabulatedCoefficients->setForceCoefficientsFiles( forceCoefficientFiles );
    return tabulatedCoefficients;
}

//  Function to create aerodynamic coefficient settings from coefficients stored in data files
/*  
 *  Function to create aerodynamic coefficient settings from coefficients stored in data files. Separate files are defined for
 *  the three components of the force coefficients.  The file format is discussed in the Tudat wiki
 *  \param forceCoefficientFiles List (size 3) of files containing the aerodynamic force coefficients
 *  \param momentCoefficientFiles List (size 3) of files containing the aerodynamic moment coefficients
 *  \param referenceLength Reference length with which aerodynamic moments
 *  (about x- and z- axes) are non-dimensionalized.
 *  \param referenceArea Reference area with which aerodynamic forces and moments are
 *  non-dimensionalized.
 *  \param lateralReferenceLength Reference length with which aerodynamic moments (about y-axis)
 *  is non-dimensionalized.
 *  \param momentReferencePoint Point w.r.t. aerodynamic moment is calculated
 *  \param independentVariableNames Physical meaning of the independent variables of the aerodynamic coefficients
 *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
 *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positive along the positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
 *  coefficients are typically defined in negative direction.
 *  \param interpolatorSettings Pointer to an interpolator settings object, where the
 *  conditions for interpolation are saved.
 *  \return Settings for creation of aerodynamic coefficient interface, based on contents read from files defined in
 *  forceCoefficientFiles and reference data given as input.
 */
std::shared_ptr< AerodynamicCoefficientSettings > readTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr );

//  Function to create aerodynamic coefficient settings from coefficients stored in data files
/*  
 * Function to create aerodynamic coefficient settings from coefficients stored in data files. Separate files are defined for
 * the three components of the force coefficients. From this function, no moment coefficients are read (set to zero for all
 * cases). The file format is discussed in the Tudat wiki
 * \param forceCoefficientFiles List (size 3) of files containing the aerodynamic coefficients
 * \param referenceArea Reference area of aerodynamic coefficients
 * \param independentVariableNames Physical meaning of the independent variables of the aerodynamic coefficients
 *  \param areCoefficientsInAerodynamicFrame Boolean to define whether the aerodynamic
 *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positive along the positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
 *  coefficients are typically defined in negative direction.
 *  \param interpolatorSettings Pointer to an interpolator settings object, where the
 *  conditions for interpolation are saved.
 *  \return Settings for creation of aerodynamic coefficient interface, based on contents read from files defined in
 *  forceCoefficientFiles and reference data given as input.
 */
std::shared_ptr< AerodynamicCoefficientSettings >
readTabulatedAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const double referenceArea,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const bool areCoefficientsInAerodynamicFrame = true,
        const bool areCoefficientsInNegativeAxisDirection = true,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr );

//  Function to create an aerodynamic coefficient interface containing constant coefficients.
/*  
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
 *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positiver along tyhe positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
 *  coefficients are typically defined in negative direction.
 *  \return Aerodynamic coefficient interface with constant coefficients.
 */
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createConstantCoefficientAerodynamicCoefficientInterface(
        const Eigen::Vector3d constantForceCoefficient,
        const Eigen::Vector3d constantMomentCoefficient,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame = false,
        const bool areCoefficientsInNegativeAxisDirection = true );

std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createZeroParameterAerodynamicCoefficientInterface(
        const std::function< Eigen::Vector3d( ) > constantForceCoefficientFunction,
        const std::function< Eigen::Vector3d( ) > constantMomentCoefficientFunction,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame = false,
        const bool areCoefficientsInNegativeAxisDirection = true );

//  Factory function for tabulated (N-D independent variables) aerodynamic coefficient interface.
/*  
 *  Factory function for tabulated (N-D independent variables) aerodynamic coefficient interface.
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
 *  coefficients are defined in the aerodynamic frame (drag, side, lift force) or in the body
 *  frame (typically denoted as Cx, Cy, Cz).
 *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
 *  coefficients are positive along the positive axes of the body or aerodynamic frame
 *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
 *  coefficients are typically defined in negative direction.
 *  \param interpolatorSettings Pointer to an interpolator settings object, where the
 *  conditions for interpolation are saved.
 *  \return Tabulated aerodynamic coefficient interface pointer.
 */
template< unsigned int NumberOfDimensions >
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createTabulatedCoefficientAerodynamicCoefficientInterface(
        const std::vector< std::vector< double > > independentVariables,
        const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > forceCoefficients,
        const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > momentCoefficients,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames,
        const double referenceLength,
        const double referenceArea,
        const double lateralReferenceLength,
        const Eigen::Vector3d& momentReferencePoint,
        const bool areCoefficientsInAerodynamicFrame = false,
        const bool areCoefficientsInNegativeAxisDirection = true,
        const std::shared_ptr< interpolators::InterpolatorSettings > interpolatorSettings = nullptr )
{
    using namespace tudat::interpolators;

    // Check input consistency.
    if( independentVariables.size( ) != NumberOfDimensions )
    {
        throw std::runtime_error( "Error when creating tabulated aerodynamic coefficient interface, "
                                  "inconsistent variable vector dimensioning" );
    }

    if( independentVariableNames.size( ) != NumberOfDimensions )
    {
        throw std::runtime_error( "Error when creating tabulated aerodynamic coefficient interface, "
                                  "inconsistent variable name vector dimensioning" );
    }

    // Create interpolators for coefficients.
    std::shared_ptr< MultiDimensionalInterpolator < double, Eigen::Vector3d, NumberOfDimensions > > forceInterpolator;
    std::shared_ptr< MultiDimensionalInterpolator< double, Eigen::Vector3d, NumberOfDimensions > > momentInterpolator;
    if ( interpolatorSettings == nullptr )
    {
        forceInterpolator = createMultiDimensionalInterpolator< double, Eigen::Vector3d, NumberOfDimensions >(
                    independentVariables, forceCoefficients,
                    std::make_shared< InterpolatorSettings >( multi_linear_interpolator, huntingAlgorithm, false,
                                                              std::vector< BoundaryInterpolationType >( NumberOfDimensions,
                                                                                                        use_boundary_value ) ) );
        momentInterpolator = createMultiDimensionalInterpolator< double, Eigen::Vector3d, NumberOfDimensions >(
                    independentVariables, momentCoefficients,
                    std::make_shared< InterpolatorSettings >( multi_linear_interpolator, huntingAlgorithm, false,
                                                              std::vector< BoundaryInterpolationType >( NumberOfDimensions,
                                                                                                        use_boundary_value ) ) );
    }
    else
    {
        forceInterpolator = createMultiDimensionalInterpolator< double, Eigen::Vector3d, NumberOfDimensions >(
                    independentVariables, forceCoefficients, interpolatorSettings );
        momentInterpolator = createMultiDimensionalInterpolator< double, Eigen::Vector3d, NumberOfDimensions >(
                    independentVariables, momentCoefficients, interpolatorSettings );
    }

    // Create aerodynamic coefficient interface.
    return std::make_shared< aerodynamics::CustomAerodynamicCoefficientInterface >(
                std::bind( &MultiDimensionalInterpolator< double, Eigen::Vector3d, NumberOfDimensions >::interpolate,
                           forceInterpolator, std::placeholders::_1 ),
                std::bind( &MultiDimensionalInterpolator< double, Eigen::Vector3d, NumberOfDimensions >::interpolate,
                           momentInterpolator, std::placeholders::_1 ),
                referenceLength, referenceArea, lateralReferenceLength, momentReferencePoint,
                independentVariableNames,
                areCoefficientsInAerodynamicFrame, areCoefficientsInNegativeAxisDirection );
}

//  Factory function for tabulated (1-D independent variables) aerodynamic coefficient interface from coefficient settings.
/*  
 *  Factory function for tabulated (1-D independent variables) aerodynamic coefficient interface from coefficient settings.
 *  \param coefficientSettings Settings for aerodynamic coefficient interface, must be of derived
 *  type TabulatedAerodynamicCoefficientSettings< 1 >
 *  \param body Name of body for which coefficient interface is to be made.
 *  \return Tabulated aerodynamic coefficient interface pointer.
 */
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createUnivariateTabulatedCoefficientAerodynamicCoefficientInterface(
        const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body );

//  Factory function for tabulated aerodynamic coefficient interface from coefficient settings.
/*  
 *  Factory function for tabulated aerodynamic coefficient interface from coefficient settings.
 *  This function is included to allow easier interface between the non-templated general
 *  createAerodynamicCoefficientInterface and the templated
 *  createTabulatedCoefficientAerodynamicCoefficientInterface.
 *  \param coefficientSettings Settings for aerodynamic coefficient interface, must be of derived
 *  type TabulatedAerodynamicCoefficientSettings< NumberOfDimensions >/
 *  \param body Name of body for which coefficient interface is to be made.
 *  \return Tabulated aerodynamic coefficient interface pointer.
 */
template< unsigned int NumberOfDimensions >
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createTabulatedCoefficientAerodynamicCoefficientInterface(
        const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    // Check consistency of type.
    std::shared_ptr< TabulatedAerodynamicCoefficientSettings< NumberOfDimensions > > tabulatedCoefficientSettings =
            std::dynamic_pointer_cast< TabulatedAerodynamicCoefficientSettings< NumberOfDimensions > >(
                coefficientSettings );
    if( tabulatedCoefficientSettings == nullptr )
    {
        throw std::runtime_error(
                    "Error, expected tabulated aerodynamic coefficients of size " +
                    std::to_string( NumberOfDimensions ) + "for body " + body );
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
                    tabulatedCoefficientSettings->getAreCoefficientsInNegativeAxisDirection( ),
                    tabulatedCoefficientSettings->getInterpolatorSettings( ) );
    }
}

//  Function to create an aerodynamic coefficient interface.
/*  
 * Function to create an aerodynamic coefficient interface from interface settings.
 * \param coefficientSettings Settings for the aerodynamic coefficient interface.
 * \param body Name of body for which aerodynamic coefficients are to be made.
 * \return Aerodynamic coefficient interface pointer of reqyested type and settings.
 */
std::shared_ptr< aerodynamics::AerodynamicCoefficientInterface >
createAerodynamicCoefficientInterface(
        const std::shared_ptr< AerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body );


} // simulation_setup

} // tudat

#endif // TUDAT_CREATEAERODYNAMICCOEFFICIENTINTERFACE_H
