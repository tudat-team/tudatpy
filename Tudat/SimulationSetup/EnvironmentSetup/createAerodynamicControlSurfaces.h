/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEAERODYNAMICCONTROLSURFACES_H
#define TUDAT_CREATEAERODYNAMICCONTROLSURFACES_H

#include <vector>

#include <memory>
#include <boost/make_shared.hpp>

#include "Tudat/InputOutput/aerodynamicCoefficientReader.h"
#include "Tudat/Astrodynamics/Aerodynamics/controlSurfaceAerodynamicCoefficientInterface.h"
#include "Tudat/Mathematics/Interpolators/multiLinearInterpolator.h"

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

//! Class for providing settings for aerodynamic coefficient model of control surfaces.
/*!
 *  Class for providing settings for automatic  control surfaces aerodynamic coefficient model creation. This class is
 *  a functional (base) class for settings of aerodynamic coefficient models that require no
 *  information in addition to their type. Aerodynamic coefficient model classes defining requiring
 *  additional information must be created using an object derived from this class. The
 *  aerodynamic reference quantities of this interface are not stored here, but are instead implicitly assumed to be equal
 *  to the AerodynamicCoefficientSettings hosting it.
 */
class ControlSurfaceIncrementAerodynamicCoefficientSettings
{
public:

    //! Constructor, sets type of aerodynamic coefficient model.
    /*!
     *  Constructor, sets type of aerodynamic coefficient model. Settings for aerodynamic
     *  coefficient models requiring additional information should be defined in a derived class.
     *  \param aerodynamicCoefficientType Type of aerodynamic coefficient model that is to be
     *  created.
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients (must contain exactly one entry of
     *  control_surface_deflection_dependent).
     */
    ControlSurfaceIncrementAerodynamicCoefficientSettings(
            const AerodynamicCoefficientTypes aerodynamicCoefficientType,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames ):
        aerodynamicCoefficientType_( aerodynamicCoefficientType ),
        independentVariableNames_( independentVariableNames )
    {
        // Check input type consistency
        if( aerodynamicCoefficientType_ == hypersonic_local_inclincation_coefficients )
        {
            throw std::runtime_error( "Error, hypersonic local inclination control surface increments not yet available" );
        }

        if( aerodynamicCoefficientType_ == constant_aerodynamic_coefficients )
        {
            throw std::runtime_error( "Error, constant control surface increments not available, should depend at least on deflection angle" );
        }

        if( std::count( independentVariableNames.begin( ), independentVariableNames.end( ),
                        aerodynamics::control_surface_deflection_dependent ) != 1 )
        {
            std::cerr << "Warning when creating ControlSurfaceIncrementAerodynamicCoefficientSettings, expected single dependency on control surface deflections, found " <<
                std::count( independentVariableNames.begin( ), independentVariableNames.end( ),
                            aerodynamics::control_surface_deflection_dependent ) << ", dependencies" << std::endl;
        }
    }

    //! Destructor
    virtual ~ControlSurfaceIncrementAerodynamicCoefficientSettings( ){ }

    //! Function to return type of aerodynamic coefficient model that is to be created.
    /*!
     *  Function to return type of aerodynamic coefficient model that is to be created.
     *  \return Type of aerodynamic coefficient model that is to be created.
     */
    AerodynamicCoefficientTypes getAerodynamicCoefficientType( )
    {
        return aerodynamicCoefficientType_;
    }

       //! Function to return identifiers of physical meaning of independent variables.
    /*!
     *  Function to return identifiers of physical meaning of independent variables.
     *  \return Identifiers of physical meaning of independent variables.
     */
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > getIndependentVariableNames( )
    {
        return independentVariableNames_;
    }


protected:

    //! Type of aerodynamic coefficient model that is to be created.
    AerodynamicCoefficientTypes aerodynamicCoefficientType_;

    //! Identifiers of physical meaning of independent variables.
    std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames_;
};

//! Object for setting aerodynamic coefficients from a user-defined N-dimensional table
/*!
 *  Object for setting aerodynamic coefficients from a user-defined N-dimensional table (with N>1).
 *  The user must provide the force (and moment) coefficients in boost multi_arrays, and
 *  define the physical meaning of each of the independent variables.
 */
template< unsigned int NumberOfDimensions >
class TabulatedControlSurfaceIncrementAerodynamicCoefficientSettings:
        public ControlSurfaceIncrementAerodynamicCoefficientSettings
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
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients (must contain exactly one entry of
     *  control_surface_deflection_dependent).
     */
    TabulatedControlSurfaceIncrementAerodynamicCoefficientSettings(
            const std::vector< std::vector< double > > independentVariables,
            const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > forceCoefficients,
            const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > momentCoefficients,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >  independentVariableNames ):
        ControlSurfaceIncrementAerodynamicCoefficientSettings( tabulated_coefficients, independentVariableNames ),
        independentVariables_( independentVariables ),
        forceCoefficients_( forceCoefficients ),
        momentCoefficients_( momentCoefficients )
    {

    }

    //! Constructor, sets properties of force coefficients (and zero moment coefficients).
    /*!
     *  Constructor, sets properties of force coefficients (and zero moment coefficients).
     *  \param independentVariables Values of indepependent variables at which the coefficients
     *  in the input multi arrays are defined.
     *  \param forceCoefficients Values of force coefficients at independent variables defined
     *  by independentVariables.
     *  \param independentVariableNames Vector with identifiers the physical meaning of each
     *  independent variable of the aerodynamic coefficients (must contain exactly one entry of
     *  control_surface_deflection_dependent).
     */
    TabulatedControlSurfaceIncrementAerodynamicCoefficientSettings(
            const std::vector< std::vector< double > > independentVariables,
            const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > forceCoefficients,
            const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >  independentVariableNames ):
        ControlSurfaceIncrementAerodynamicCoefficientSettings( tabulated_coefficients, independentVariableNames ),
        independentVariables_( independentVariables ),
        forceCoefficients_( forceCoefficients )
    {
        // Create zero moment coefficients at same data points as force coefficients
        std::vector< size_t > sizeVector;
        const size_t* arrayShape = forceCoefficients_.shape( );
        sizeVector.assign( arrayShape, arrayShape + forceCoefficients_.num_dimensions( ) );
        momentCoefficients_.resize( sizeVector );
        std::fill( momentCoefficients_.data( ), momentCoefficients_.data( ) + momentCoefficients_.num_elements( ),
                   Eigen::Vector3d::Zero( ) );
    }

    //! Destructor
    ~TabulatedControlSurfaceIncrementAerodynamicCoefficientSettings( ){ }

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
    boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > getForceCoefficients( )
    {
        return forceCoefficients_;
    }

    //! Function to return values of moment coefficients in table.
    /*!
     * Function to return values of moment coefficients in table.
     * \return Values of moment coefficients in table.
     */
    boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > getMomentCoefficients( )
    {
        return momentCoefficients_;
    }

protected:

    //! Values of the indepependent variables of tables of coefficients.
    const std::vector< std::vector< double > > independentVariables_;

    //! Values of force coefficients in table.
    boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > forceCoefficients_;

    //! Values of moment coefficients in table.
    boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > momentCoefficients_;
};

//! Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files
/*!
 *  Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files.
 *  Separate files are defined for
 *  the three components of the force coefficients.  The file format is discussed in the Tudat wiki
 *  Note that this function requires the number of independent variables in the coefficient files to be known. If this is not
 *  the case, the readTabulatedAerodynamicCoefficientsFromFiles function should be used.
 *  \param forceCoefficientFiles List (size 3) of files containing the aerodynamic force coefficients
 *  \param momentCoefficientFiles List (size 3) of files containing the aerodynamic moment coefficients
 *  \param independentVariableNames Physical meaning of the independent variables of the aerodynamic coefficient
 *  \return Settings for creation of control surface aerodynamic coefficient interface, based on contents read from files
 *  defined in forceCoefficientFiles and momentCoefficientFiles
 */
template< unsigned int NumberOfIndependentVariables >
std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings >
readGivenSizeTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames )
{
    std::pair< boost::multi_array< Eigen::Vector3d, NumberOfIndependentVariables >, std::vector< std::vector< double > > >
            aerodynamicForceCoefficients = input_output::readAerodynamicCoefficients< NumberOfIndependentVariables >(
                forceCoefficientFiles );
    std::pair< boost::multi_array< Eigen::Vector3d, NumberOfIndependentVariables >, std::vector< std::vector< double > > >
            aerodynamicMomentCoefficients = input_output::readAerodynamicCoefficients< NumberOfIndependentVariables >(
                momentCoefficientFiles );

    if( !input_output::compareIndependentVariables(
                aerodynamicForceCoefficients.second, aerodynamicMomentCoefficients.second ) )
    {
        throw std::runtime_error( "Error when creating aerodynamic coefficient settings from file, force and moment independent variables are inconsistent" );
    }

    if( independentVariableNames.size( ) != NumberOfIndependentVariables )
    {
        throw std::runtime_error( "Error when creating aerodynamic coefficient settings from file, input sizes are inconsistent" );
    }

    return std::make_shared< TabulatedControlSurfaceIncrementAerodynamicCoefficientSettings<
            NumberOfIndependentVariables > >(
                aerodynamicForceCoefficients.second, aerodynamicForceCoefficients.first, aerodynamicMomentCoefficients.first,
                independentVariableNames );
}

//! Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files
/*!
 *  Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files. Separate files
 *  are defined for the three components of the force coefficients. From this function, no moment coefficients are read
 *  (set to zero for all cases). The file format is discussed in the Tudat wiki
 *  Note that this function requires the number of independent variables in the coefficient files to be known. If this is not
 *  the case, the readTabulatedAerodynamicCoefficientsFromFiles function should be used.
 *  \param forceCoefficientFiles List (size 3) of files containing the aerodynamic coefficients
 *  \param independentVariableNames Physical meaning of the independent variables of the aerodynamic coefficients
 *  \return Settings for creation of control surface aerodynamic coefficient interface, based on contents read from files
 *   defined in forceCoefficientFiles and reference data given as input
 */
template< unsigned int NumberOfIndependentVariables >
std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings >
readGivenSizeTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames )
{
    std::pair< boost::multi_array< Eigen::Vector3d, NumberOfIndependentVariables >, std::vector< std::vector< double > > >
            aerodynamicCoefficients = input_output::readAerodynamicCoefficients< NumberOfIndependentVariables >(
                forceCoefficientFiles );

    // Check input consistency
    if( independentVariableNames.size( ) != NumberOfIndependentVariables )
    {
        throw std::runtime_error( "Error when creating aerodynamic coefficient settings from file, input sizes are inconsistent" );
    }

    // Create coefficient settings.
    return std::make_shared< TabulatedControlSurfaceIncrementAerodynamicCoefficientSettings<
            NumberOfIndependentVariables > >(
                aerodynamicCoefficients.second, aerodynamicCoefficients.first, independentVariableNames );
}

//! Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files
/*!
 *  Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files.
 *  Separate files are defined for the three components of the force coefficients.  The file format is discussed in
 *  the Tudat wiki
 *  \param forceCoefficientFiles List  (containing entries at key = 0,1 and/or 2, denoting coefficients in x, y and z
 *  direction) of files containing the aerodynamic force coefficients
 *  \param momentCoefficientFiles List (containing entries at key = 0,1 and/or 2, denoting coefficients in x, y and z
 *  direction) of files containing the aerodynamic moment coefficients
 *  \param independentVariableNames Physical meaning of the independent variables of the aerodynamic coefficients
 *  \return Settings for creation of control surface aerodynamic coefficient interface, based on contents read from files
 *  defined in forceCoefficientFiles and momentCoefficientFiles.
 */
std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings >
readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::map< int, std::string > momentCoefficientFiles,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames );

//! Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files
/*!
 * Function to create control surface aerodynamic coefficient settings fom coefficients stored in data files. Separate files
 * are defined for the three components of the force coefficients. From this function, no moment coefficients are read
 * (set to zero for all cases). The file format is discussed in the Tudat wiki
 * \param forceCoefficientFiles List  (containing entries at key = 0,1 and/or 2, denoting coefficients in x, y and z
 * direction) of files containing the aerodynamic force coefficients.
 * \param independentVariableNames Physical meaning of the independent variables of the aerodynamic coefficients
 * \return Settings for creation of control surface aerodynamic coefficient interface, based on contents read from
 * files defined in forceCoefficientFiles
 */
std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings >
readTabulatedControlIncrementAerodynamicCoefficientsFromFiles(
        const std::map< int, std::string > forceCoefficientFiles,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables > independentVariableNames );

//! Function to create control surface aerodynamic coefficient settings from user-defined coefficients.
/*!
 * Function to create control surface aerodynamic coefficient settings from crom user-defined coefficients.
 * \param independentVariables Values of indepependent variables at which the coefficients
 * in the input multi arrays are defined.
 * \param forceCoefficients Values of force coefficients at independent variables defined
 * by independentVariables.
 * \param momentCoefficients Values of moment coefficients at independent variables defined
 * by independentVariables.
 * \param independentVariableNames Vector with identifiers the physical meaning of each
 * independent variable of the aerodynamic coefficients (must contain exactly one entry of
 * control_surface_deflection_dependent)
 * \return Settings for creation of control surface aerodynamic coefficient interface, based on contents read from
 * files defined in forceCoefficientFiles
 */
template< unsigned int NumberOfDimensions >
std::shared_ptr< aerodynamics::ControlSurfaceIncrementAerodynamicInterface >
createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface(
        const std::vector< std::vector< double > > independentVariables,
        const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > forceCoefficients,
        const boost::multi_array< Eigen::Vector3d, static_cast< size_t >( NumberOfDimensions ) > momentCoefficients,
        const std::vector< aerodynamics::AerodynamicCoefficientsIndependentVariables >
        independentVariableNames )
{
    // Check input consistency.
    if( independentVariables.size( ) != NumberOfDimensions )
    {
        throw std::runtime_error( "Error when creating tabulated control surface increment aerodynamic coefficient interface, inconsistent variable vector dimensioning" );
    }

    if( independentVariableNames.size( ) != NumberOfDimensions )
    {
       throw std::runtime_error( "Error when creating tabulated control surface increment aerodynamic coefficient interface, inconsistent variable name vector dimensioning" );

    }

    // Create interpolators for coefficients.
    std::shared_ptr< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > > forceInterpolator =
            std::make_shared< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > >(
                independentVariables, forceCoefficients );
    std::shared_ptr< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > > momentInterpolator =
            std::make_shared< interpolators::MultiLinearInterpolator
            < double, Eigen::Vector3d, NumberOfDimensions > >(
                independentVariables, momentCoefficients );

    // Create aerodynamic coefficient interface.
    return  std::make_shared< aerodynamics::CustomControlSurfaceIncrementAerodynamicInterface >(
                std::bind( &interpolators::MultiLinearInterpolator
                             < double, Eigen::Vector3d, NumberOfDimensions >::interpolate,
                             forceInterpolator, std::placeholders::_1 ),
                std::bind( &interpolators::MultiLinearInterpolator
                             < double, Eigen::Vector3d, NumberOfDimensions >::interpolate,
                             momentInterpolator, std::placeholders::_1 ),
                independentVariableNames );
}

//! Function to create tabulated control surface aerodynamic coefficients from associated settings object
/*!
 *  Function to create tabulated control surface aerodynamic coefficients from associated settings object
 *  \param coefficientSettings Object containing settings for tabulated control surface aerodynamic coefficients
 *  \param body Name of body for which coefficients are to be created.
 *  \return Object used to compute/update control surface aerodynamics during propagation.
 */
template< unsigned int NumberOfDimensions >
std::shared_ptr< aerodynamics::ControlSurfaceIncrementAerodynamicInterface >
createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface(
        const std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body )
{
    // Check consistency of type.
    std::shared_ptr< TabulatedControlSurfaceIncrementAerodynamicCoefficientSettings< NumberOfDimensions > >
            tabulatedCoefficientSettings = std::dynamic_pointer_cast<
            TabulatedControlSurfaceIncrementAerodynamicCoefficientSettings< NumberOfDimensions > >( coefficientSettings );
    if( tabulatedCoefficientSettings == nullptr )
    {
        throw std::runtime_error(
                    "Error, expected tabulated control surface increment aerodynamic coefficients of size " +
                    std::to_string( NumberOfDimensions ) + "for body " + body );
    }
    else
    {
        return createTabulatedControlSurfaceIncrementAerodynamicCoefficientInterface< NumberOfDimensions >(
                    tabulatedCoefficientSettings->getIndependentVariables( ),
                    tabulatedCoefficientSettings->getForceCoefficients( ),
                    tabulatedCoefficientSettings->getMomentCoefficients( ),
                    tabulatedCoefficientSettings->getIndependentVariableNames( ) );
    }
}

//! Function to tabulated control surface aerodynamic coefficients from associated settings object
/*!
 *  Function to tabulated control surface aerodynamic coefficients from associated settings object
 *  \param coefficientSettings Object containing settings for control surface aerodynamic coefficients
 *  \param body Name of body for which coefficients are to be created.
 *  \return Object used to compute/update control surface aerodynamics during propagation.
 */
std::shared_ptr< aerodynamics::ControlSurfaceIncrementAerodynamicInterface >
createControlSurfaceIncrementAerodynamicCoefficientInterface(
        const std::shared_ptr< ControlSurfaceIncrementAerodynamicCoefficientSettings > coefficientSettings,
        const std::string& body );

}

}

#endif // TUDAT_CREATEAERODYNAMICCONTROLSURFACES_H
