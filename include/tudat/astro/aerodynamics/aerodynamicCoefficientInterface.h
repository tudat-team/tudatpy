/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    Notes
 *      The computeCoefficients() function is not yet implemented in any derived classes and is
 *      therefore not pure virtual in this base class.
 *
 */

#ifndef TUDAT_AERODYNAMIC_COEFFICIENT_INTERFACE_H
#define TUDAT_AERODYNAMIC_COEFFICIENT_INTERFACE_H

#include <vector>
#include <map>

#include <memory>
#include <Eigen/Core>

#include "tudat/astro/aerodynamics/controlSurfaceAerodynamicCoefficientInterface.h"
#include "tudat/astro/aerodynamics/aerodynamics.h"
#include "tudat/astro/reference_frames/referenceFrameTransformations.h"
#include "tudat/basics/utilities.h"

namespace tudat
{
namespace aerodynamics
{

enum AerodynamicCoefficientFrames
{
    body_fixed_frame_coefficients,
    negative_body_fixed_frame_coefficients,
    negative_aerodynamic_frame_coefficients,
    positive_aerodynamic_frame_coefficients,
    undefined_frame_coefficients
};

inline std::pair< reference_frames::AerodynamicsReferenceFrames, double > convertCoefficientFrameToGeneralAerodynamicFrame(
        const AerodynamicCoefficientFrames coefficientFrame )
{
    std::pair< reference_frames::AerodynamicsReferenceFrames, double > frameConversion;
    switch( coefficientFrame )
    {
        case body_fixed_frame_coefficients:
            frameConversion = std::make_pair( reference_frames::body_frame, 1.0 );
            break;
        case negative_body_fixed_frame_coefficients:
            frameConversion = std::make_pair( reference_frames::body_frame, -1.0 );
            break;
        case negative_aerodynamic_frame_coefficients:
            frameConversion = std::make_pair( reference_frames::aerodynamic_frame, -1.0 );
            break;
        case positive_aerodynamic_frame_coefficients:
            frameConversion = std::make_pair( reference_frames::aerodynamic_frame, 1.0 );
            break;
        default:
            throw std::runtime_error( "Error when retrieving aerodynamic coefficient frame defintion; undefined frame used as input" );
    }
    return frameConversion;
}

inline reference_frames::AerodynamicsReferenceFrames getCompleteFrameForCoefficients(
        const AerodynamicCoefficientFrames coefficientsFrame )
{
    reference_frames::AerodynamicsReferenceFrames completeFrame;
    switch( coefficientsFrame )
    {
        case body_fixed_frame_coefficients:
            completeFrame = reference_frames::body_frame;
            break;
        case negative_body_fixed_frame_coefficients:
            completeFrame = reference_frames::body_frame;
            break;
        case negative_aerodynamic_frame_coefficients:
            completeFrame = reference_frames::aerodynamic_frame;
            break;
        case positive_aerodynamic_frame_coefficients:
            completeFrame = reference_frames::aerodynamic_frame;
            break;
        default:
            throw std::runtime_error( "Error when getting aerodynamic coefficient frame; frame not rezognized" );

    }
    return completeFrame;
}

inline bool areCoefficientsInNegativeDirection(
        const AerodynamicCoefficientFrames coefficientsFrame )
{
    bool areCoefficientsInNegativeDirection;
    switch( coefficientsFrame )
    {
        case body_fixed_frame_coefficients:
            areCoefficientsInNegativeDirection = false;
            break;
        case negative_body_fixed_frame_coefficients:
            areCoefficientsInNegativeDirection = true;
            break;
        case negative_aerodynamic_frame_coefficients:
            areCoefficientsInNegativeDirection = true;
            break;
        case positive_aerodynamic_frame_coefficients:
            areCoefficientsInNegativeDirection = false;
            break;
        default:
            throw std::runtime_error( "Error when getting aerodynamic coefficient frame; frame not rezognized" );

    }
    return areCoefficientsInNegativeDirection;
}


inline AerodynamicCoefficientFrames getAerodynamicCoefficientFrame(
    const bool areCoefficientsInAerodynamicFrame = true,
    const bool areCoefficientsInNegativeAxisDirection = true )
{
    AerodynamicCoefficientFrames coefficientsFrame = undefined_frame_coefficients;
    if( areCoefficientsInAerodynamicFrame && areCoefficientsInNegativeAxisDirection )
    {
        coefficientsFrame = negative_aerodynamic_frame_coefficients;
    }
    else if( areCoefficientsInAerodynamicFrame && !areCoefficientsInNegativeAxisDirection )
    {
        coefficientsFrame = positive_aerodynamic_frame_coefficients;
    }
    else if( !areCoefficientsInAerodynamicFrame && !areCoefficientsInNegativeAxisDirection )
    {
        coefficientsFrame = body_fixed_frame_coefficients;
    }
    else if( !areCoefficientsInAerodynamicFrame && areCoefficientsInNegativeAxisDirection )
    {
        coefficientsFrame = negative_body_fixed_frame_coefficients;
    }
    else
    {
        throw std::runtime_error( "Error, could not define aerodynamic coefficient frame from booleans." );
    }
    return coefficientsFrame;
}

struct AerodynamicMomentContributionInterface
{
    AerodynamicMomentContributionInterface( const std::function< Eigen::Matrix3d( ) > forceToMomentFrameRotation,
                                            const std::function< Eigen::Matrix3d( ) > bodyFixedToMomentFrameRotation,
                                            const std::function< Eigen::Vector3d( ) > centerOfMassPosition,
                                            const double signMultiplier ):
            forceToMomentFrameRotation_( forceToMomentFrameRotation ),
            bodyFixedToMomentFrameRotation_( bodyFixedToMomentFrameRotation ),
            centerOfMassPosition_( centerOfMassPosition ),
            signMultiplier_( signMultiplier )
            { }

    Eigen::Vector3d getMomentCoefficientsCorrection(
        const Eigen::Vector3d& momentReferencePoint,
        const Eigen::Vector3d& forceCoefficients,
        const double referenceLength );

    std::function< Eigen::Matrix3d( ) > forceToMomentFrameRotation_;
    std::function< Eigen::Matrix3d( ) > bodyFixedToMomentFrameRotation_;
    std::function< Eigen::Vector3d( ) > centerOfMassPosition_;
    double signMultiplier_;
};

//! Base class to hold an aerodynamic coefficient interface.
/*!
 * This interface can, for instance, be a database of coefficients or an aerodynamic analysis code
 * which generates coefficients. The aerodynamic coefficients are defined as a function of any number of independent
 * variables, the physical meaning of which is stored in the coefficient interface.
 */
class AerodynamicCoefficientInterface
{
public:

    friend class ScaledAerodynamicCoefficientInterface;

    //! Constructor.
    /*!
     *  Constructor, sets quantities common to all derived class aerodynamic coefficient interfaces.
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
     *  frame (typically denoted as Cx, Cy, Cz) (default true).
     *  \param areCoefficientsInNegativeAxisDirection Boolean to define whether the aerodynamic
     *  coefficients are positive along tyhe positive axes of the body or aerodynamic frame
     *  (see areCoefficientsInAerodynamicFrame). Note that for (drag, side, lift force), the
     *  coefficients are typically defined in negative direction (default true).
     */
    AerodynamicCoefficientInterface(
            const double referenceLength,
            const double referenceArea,
            const Eigen::Vector3d& momentReferencePoint,
            const std::vector< AerodynamicCoefficientsIndependentVariables >
            independentVariableNames,
            const AerodynamicCoefficientFrames forceCoefficientsFrame = negative_aerodynamic_frame_coefficients,
            const AerodynamicCoefficientFrames momentCoefficientsFrame = body_fixed_frame_coefficients):
        referenceLength_( referenceLength ),
        referenceArea_( referenceArea ),
        momentReferencePoint_( momentReferencePoint ),
        independentVariableNames_( independentVariableNames ),
        forceCoefficientsFrame_( forceCoefficientsFrame ),
        momentCoefficientsFrame_( momentCoefficientsFrame ),
        momentContributionInterface_( nullptr )
    {
        numberOfIndependentVariables_ = static_cast< unsigned int >( independentVariableNames.size( ) );
        referenceLengths_ = Eigen::Vector3d::Constant( referenceLength_ );
    }

    //! Default destructor.
    virtual ~AerodynamicCoefficientInterface( ) { }

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


    //! Get reference lengths.
    /*!
     * Returns the all reference lengths used to non-dimensionalize aerodynamic moments.
     * \return Aerodynamic reference lengths.
     */
    Eigen::Vector3d getReferenceLengths( ) { return referenceLengths_; }

    //! Get moment reference point.
    /*!
     * Returns the point w.r.t. which the arm of the aerodynamic moment on a vehicle panel is
     * determined.
     * \return Aerodynamic reference point.
     */
    Eigen::VectorXd getMomentReferencePoint( ) { return momentReferencePoint_; }


    //! Compute the aerodynamic coefficients of the body itself (without control surfaces) at current flight condition.
    /*!
     *  Computes the current force and moment coefficients of the body itself (without control surfaces) and is to be
     *  implemented in derived classes. Input is a set of independent variables
     *  (doubles) which represent the variables from which the coefficients are calculated
     *  \param independentVariables Independent variables of force and moment coefficient
     *  determination implemented by derived class
     *  \param currentTime Time to which coefficients are to be updated
     */
    virtual void updateCurrentCoefficients(
            const std::vector< double >& independentVariables,
            const double currentTime = TUDAT_NAN ) = 0;

    //! Function to update the aerodynamic coefficients of the full body with control surfaces
    /*!
     *  Function to update the aerodynamic coefficients of the full body with control surfaces. The full body coefficients
     *  are cimputed first, after which the control surfaces are updated and the results added to the full coefficients.
     *  \param independentVariables Independent variables of force and moment coefficient of body without control surfaces
     *  \param controlSurfaceIndependentVariables Map of independent variables of force and moment coefficient of
     *  control surfaces, with map key denoting the control surface identifier.
     *  \param currentTime Time to which coefficients are to be updated.
     */
    void updateFullCurrentCoefficients(
            const std::vector< double >& independentVariables,
            const std::map< std::string, std::vector< double > >& controlSurfaceIndependentVariables =
            std::map< std::string, std::vector< double > > ( ),
            const double currentTime = TUDAT_NAN,
            const bool addMomentContributionIfPresent = true );

    /*!
    //! Function for calculating and returning aerodynamic force coefficients
     *  Function for calculating and returning aerodynamic force coefficients.
     *  \return Force coefficients at current independent variables
     */
    Eigen::Vector3d getCurrentForceCoefficients( )
    {
        return currentForceCoefficients_;
    }

    Eigen::Vector3d& getCurrentForceCoefficientsReference( )
    {
        return currentForceCoefficients_;
    }

    Eigen::Vector3d getCurrentControlSurfaceFreeForceCoefficients( )
    {
        if( controlSurfaceIncrementInterfaces_.size( ) > 0 )
        {
            return currentControlSurfaceFreeForceCoefficients_;
        }
        else
        {
            return currentForceCoefficients_;
        }
    }

    Eigen::Vector3d getCurrentForceCoefficientIncrement( const std::string& controlSurfaceName )
    {
        if( currentControlSurfaceForceCoefficient_.count( controlSurfaceName ) == 0 )
        {
            throw std::runtime_error( "Error when getting current control surface force increment of surface "
                + controlSurfaceName + ", no such control surface found" );
        }
        else
        {
            return currentControlSurfaceForceCoefficient_.at( controlSurfaceName );
        }
    }

    Eigen::Vector3d getCurrentMomentCoefficientIncrement( const std::string& controlSurfaceName )
    {
        if( currentControlSurfaceMomentCoefficient_.count( controlSurfaceName ) == 0 )
        {
            throw std::runtime_error( "Error when getting current control surface moment increment of surface "
                                      + controlSurfaceName + ", no such control surface found" );
        }
        else
        {
            return currentControlSurfaceMomentCoefficient_.at( controlSurfaceName );
        }
    }



    //! Function for calculating and returning aerodynamic moment coefficients
    /*!
     *  Function for calculating and returning aerodynamic moment coefficients.
     *  \return Moment coefficients at current independent variables
     */
    Eigen::Vector3d getCurrentMomentCoefficients( )
    {
        return currentMomentCoefficients_;
    }

    Eigen::Vector3d getCurrentControlSurfaceFreeMomentCoefficients( )
    {
        if( controlSurfaceIncrementInterfaces_.size( ) > 0 )
        {
            return currentControlSurfaceFreeMomentCoefficients_;
        }
        else
        {
            return currentMomentCoefficients_;
        }
    }
    //! Function for calculating and returning aerodynamic force and moment coefficients
    /*!
     *  Function for calculating and returning aerodynamic force and moment coefficients
     *  \return Force and moment coefficients at given independent variables
     */
    Eigen::Vector6d getCurrentAerodynamicCoefficients(  )
    {
        Eigen::Vector6d coefficients;
        coefficients.segment( 0, 3 ) = getCurrentForceCoefficients( );
        coefficients.segment( 3, 3 ) = getCurrentMomentCoefficients( );
        return coefficients;
    }

    //! Function to return the identifiers of the physical meaning of each independent variable.
    /*!
     *  Function to return the identifiers of the physical meaning of each independent variable
     *  of the aerodynamic coefficient interface.
     *  \return A vector with the identifiers of the physical meaning of each independent variable.
     */
    std::vector< AerodynamicCoefficientsIndependentVariables > getIndependentVariableNames( )
    {
        return independentVariableNames_;
    }

    //! Function to return a single identifier of the physical meaning of one independent variable.
    /*!
     *  Function to return a single identifier of the physical meaning of one of the independent
     *  independent variable of the coefficient interface. The index of the variable is defined
     *  by the input variable.
     *  \param index Index of list of identfiers to return
     *  \return The identifiers of the physical meaning of the independent variable at the position
     *  of the input variable.
     */
    AerodynamicCoefficientsIndependentVariables getIndependentVariableName(
            const unsigned int index )
    {
        if( index >= numberOfIndependentVariables_ )
        {
            throw std::runtime_error(
                        std::string( "Error when retrieving aerodynamic coefficient interface " ) +
                        ( " variable name, requested variable index " ) +
                        std::to_string( index ) +
                        ", but only " + std::to_string(
                            numberOfIndependentVariables_ ) + " variables available." );
        }

        return independentVariableNames_.at( index );
    }

    //! Function to return the number of independent variables upon which the coeficients depend.
    /*!
     *  Function to return the number of independent variables upon which the coeficients depend.
     *  The size of the vector used as input for updateCurrentCoefficients should always have the
     *  size returned by this variable.
     *  \return Number of independent variables upon which the coeficients depend
     */
    unsigned int getNumberOfIndependentVariables( )
    {
        return numberOfIndependentVariables_;
    }


    AerodynamicCoefficientFrames getForceCoefficientsFrame( )
    {
        return forceCoefficientsFrame_;
    }

    AerodynamicCoefficientFrames getMomentCoefficientsFrame( )
    {
        return momentCoefficientsFrame_;
    }


    void resetForceCoefficientsFrame( const AerodynamicCoefficientFrames forceCoefficientsFrame )
    {
        forceCoefficientsFrame_ = forceCoefficientsFrame;
    }

    void resetMomentCoefficientsFrame( const AerodynamicCoefficientFrames momentCoefficientsFrame )
    {
        momentCoefficientsFrame_ = momentCoefficientsFrame;
    }


    //! Function to set the list of control surface aerodynamic coefficient interfaces
    /*!
     * Function to set the list of control surface aerodynamic coefficient interfaces
     * \param controlSurfaceIncrementInterfaces Map of ControlSurfaceIncrementAerodynamicInterface pointers, wach pointer
     * denoting the coefficient interface of a single control sureface, where the map key denotes the surface's name.
     */
    void setControlSurfaceIncrements(
            const std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > >
            controlSurfaceIncrementInterfaces )
    {
        if( controlSurfaceIncrementInterfaces_.size( ) != 0 )
        {
            std::cerr<<"Warning, overwriting existing control surface increments"<<std::endl;
        }
        controlSurfaceIncrementInterfaces_ = controlSurfaceIncrementInterfaces;
        controlSurfaceNames_ = utilities::createVectorFromMapKeys( controlSurfaceIncrementInterfaces_ );
    }

    void addControlSurface(
            const std::string& controlSurfaceName,
            const std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > controlSurfaceIncrementInterface )
    {
        controlSurfaceNames_.push_back( controlSurfaceName );
        controlSurfaceIncrementInterfaces_[ controlSurfaceName ] = controlSurfaceIncrementInterface;
    }

    //! Function to get control surface name at given index in list of control surfaces
    /*!
     * Function to get control surface name at given index in list of control surfaces
     * \param index Index in list of control surfaces (controlSurfaceNames_)
     * \return Name of requested control surfaces.
     */
    std::string getControlSurfaceName( const int index )
    {
        return controlSurfaceNames_.at( index );
    }

    std::vector< std::string > getControlSurfaceNames( )
    {
        return controlSurfaceNames_;
    }

    //! Function to return the number of control surfaces in current coefficient interface.
    /*!
     * Function to return the number of control surfaces in current coefficient interface.
     * \return Number of control surfaces in current coefficient interface.
     */
    unsigned int getNumberOfControlSurfaces( )
    {
        return static_cast< unsigned int >( controlSurfaceIncrementInterfaces_.size( ) );
    }

    //! Function to return the number of independent variables for a given control surface.
    /*!
     * Function to return the number of independent variables for a given control surface.
     * \param controlSurface Name of control surface for which the number of independent variables is to be retrieved.
     * \return Number of independent variables for selected control surface.
     */
    unsigned int getNumberOfControlSurfaceIndependentVariables( const std::string controlSurface )
    {
        return controlSurfaceIncrementInterfaces_.at( controlSurface )->getNumberOfIndependentVariables( );
    }

    //! Function to get the list of independent variables for all control surfaces
    /*!
     * Function to get the list of independent variables for all control surfaces
     * \return List of independent variables for all control surfaces, with map key denoting control surface identifier
     */
    std::map< std::string, std::vector< AerodynamicCoefficientsIndependentVariables > >
    getControlSurfaceIndependentVariables( )
    {
        std::map< std::string, std::vector< AerodynamicCoefficientsIndependentVariables > >
                controlSurfaceIndependentVariables;
        for( std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > >::iterator
             contolSurfaceIterator = controlSurfaceIncrementInterfaces_.begin( );
             contolSurfaceIterator != controlSurfaceIncrementInterfaces_.end( ); contolSurfaceIterator++ )
        {
            controlSurfaceIndependentVariables[ contolSurfaceIterator->first ] =
                    contolSurfaceIterator->second->getIndependentVariableNames( );
        }

        return controlSurfaceIndependentVariables;
    }

    //! Function to get the identifier for a single independent variable for a single control surface
    /*!
     * Function to get the identifier for a single independent variable for a single control surface
     * \param controlSurface Control surface identifier.
     * \param index Index of independent variable of control surface
     * \return  Identifier for a requested independent variable for requested control surface
     */
    AerodynamicCoefficientsIndependentVariables getControlSurfaceIndependentVariableName(
            const std::string& controlSurface,
            const unsigned int index )
    {
        if( controlSurfaceIncrementInterfaces_.count( controlSurface ) == 0 )
        {
            throw std::runtime_error(
                        std::string( "Error when retrieving control surface aerodynamic coefficient interface variable name, requested surface " ) +
                        controlSurface + " , requested surface not found." );
        }
        else
        {
            return controlSurfaceIncrementInterfaces_.at( controlSurface )->getIndependentVariableName(
                        index );
        }
    }

    std::shared_ptr< AerodynamicMomentContributionInterface > getMomentContributionInterface( )
    {
        return momentContributionInterface_;
    }

    void setMomentContributionInterface( const std::shared_ptr< AerodynamicMomentContributionInterface > momentContributionInterface )
    {
        momentContributionInterface_ = momentContributionInterface;
    }





protected:


     //! Compute the aerodynamic coefficients for a single control surface, and add to full configuration coefficients.
     /*!
      *  Compures the aerodynamic coefficients for a single control surface at the current flight conditions, and adds these
      *  to the coefficients of the full configuration, stored in class instance by currentForceCoefficients_ and
      *  currentMomentCoefficients_ variables.
      *  \param currentControlSurface Name of control surface that is to be updated.
      *  \param controlSurfaceIndependentVariables Current values of independent variables of force and moment coefficient
      *  of given control surface.
      */
     void updateCurrentControlSurfaceCoefficientsCoefficients(
             const std::string& currentControlSurface,
             std::vector< double > controlSurfaceIndependentVariables )
     {
         if( controlSurfaceIncrementInterfaces_.count( currentControlSurface ) == 0 )
         {
             throw std::runtime_error( "Error when updating coefficients, could not fid control surface " + currentControlSurface );
         }
         controlSurfaceIncrementInterfaces_.at( currentControlSurface )->updateCurrentCoefficients(
                 controlSurfaceIndependentVariables );

         currentControlSurfaceForceCoefficient_[ currentControlSurface ] =
                 controlSurfaceIncrementInterfaces_.at( currentControlSurface )->getCurrentForceCoefficients( );
         currentForceCoefficients_ += currentControlSurfaceForceCoefficient_.at( currentControlSurface );

         currentControlSurfaceMomentCoefficient_[ currentControlSurface ] =
                 controlSurfaceIncrementInterfaces_.at( currentControlSurface )->getCurrentMomentCoefficients( );
         currentMomentCoefficients_ += currentControlSurfaceMomentCoefficient_.at( currentControlSurface );
     }

    //! The current force coefficients.
    /*!
     * The force coefficients at the current flight condition.
     */
    Eigen::Vector3d currentForceCoefficients_;

    Eigen::Vector3d currentControlSurfaceFreeForceCoefficients_;

    std::map< std::string, Eigen::Vector3d > currentControlSurfaceForceCoefficient_;

    //! The current moment coefficients.
    /*!
     * The moment coefficients at the current flight condition.
     */
    Eigen::Vector3d currentMomentCoefficients_;

    Eigen::Vector3d currentControlSurfaceFreeMomentCoefficients_;

    std::map< std::string, Eigen::Vector3d > currentControlSurfaceMomentCoefficient_;

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

    //! Aerodynamic reference lengths.
    /*!
     * All reference lengths with which aerodynamic moments are non-dimensionalized.
     */
    Eigen::Vector3d referenceLengths_;

    //! Aerodynamic moment reference point.
    /*!
     * Point w.r.t. which the arm of the moment on a vehicle panel is determined.
     */
    Eigen::Vector3d momentReferencePoint_;

    //! Vector with identifiers for the physical meaning of each independent variable of the
    //! aerodynamic coefficients.
    std::vector< AerodynamicCoefficientsIndependentVariables > independentVariableNames_;

    //! Number of independent variables upon which the force and moment coefficients depend.
    /*!
     *  Number of independent variables upon which the force and moment coefficients depend, i.e.
     *  the length of the vectors that should be used as input to updateCurrentCoefficients
     */
    unsigned int numberOfIndependentVariables_;

    AerodynamicCoefficientFrames forceCoefficientsFrame_;

    AerodynamicCoefficientFrames momentCoefficientsFrame_;

    std::shared_ptr< AerodynamicMomentContributionInterface > momentContributionInterface_;

    //! List of control surface aerodynamic coefficient interfaces
    std::map< std::string, std::shared_ptr< ControlSurfaceIncrementAerodynamicInterface > >
    controlSurfaceIncrementInterfaces_;

    //! Explicit list of control surface names, in same order as iterator over controlSurfaceIncrementInterfaces_
    std::vector< std::string > controlSurfaceNames_;

private:
};

class ScaledAerodynamicCoefficientInterface: public AerodynamicCoefficientInterface
{
public:
    ScaledAerodynamicCoefficientInterface(
            const std::shared_ptr< AerodynamicCoefficientInterface > baseCoefficientInterface,
            const std::function< Eigen::Vector3d( const double ) > forceCoefficientScalingFunction,
            const std::function< Eigen::Vector3d( const double ) > momentCoefficientScalingFunction,
            const bool isScalingRelative = true ):
        AerodynamicCoefficientInterface(
            baseCoefficientInterface->getReferenceLength( ),
            baseCoefficientInterface->getReferenceArea( ),
            baseCoefficientInterface->getMomentReferencePoint( ),
            baseCoefficientInterface->getIndependentVariableNames( ),
            baseCoefficientInterface->getForceCoefficientsFrame( ),
            baseCoefficientInterface->getMomentCoefficientsFrame( ) ),
    baseCoefficientInterface_( baseCoefficientInterface ),
    forceCoefficientScalingFunction_( forceCoefficientScalingFunction ),
    momentCoefficientScalingFunction_( momentCoefficientScalingFunction ),
    isScalingRelative_( isScalingRelative){ }

    void updateCurrentCoefficients(
            const std::vector< double >& independentVariables,
            const double currentTime = TUDAT_NAN )
    {
        baseCoefficientInterface_->updateCurrentCoefficients( independentVariables, currentTime );

        currentForceCoefficients_ = baseCoefficientInterface_->getCurrentForceCoefficients( );
        currentMomentCoefficients_ = baseCoefficientInterface_->getCurrentMomentCoefficients( );

        if( isScalingRelative_ )
        {
            currentForceCoefficients_.array() *= forceCoefficientScalingFunction_( currentTime ).array();
            currentMomentCoefficients_.array() *= momentCoefficientScalingFunction_( currentTime ).array();
        }
        else
        {
            currentForceCoefficients_ += forceCoefficientScalingFunction_( currentTime );
            currentMomentCoefficients_ += momentCoefficientScalingFunction_( currentTime );
        }

    }

private:

    std::shared_ptr< AerodynamicCoefficientInterface > baseCoefficientInterface_;

    std::function< Eigen::Vector3d( const double ) > forceCoefficientScalingFunction_;

    std::function< Eigen::Vector3d( const double ) > momentCoefficientScalingFunction_;

    bool isScalingRelative_;

};

//! Typedef for shared-pointer to AerodynamicCoefficientInterface object.
typedef std::shared_ptr< AerodynamicCoefficientInterface >
AerodynamicCoefficientInterfacePointer;

} // namespace aerodynamics

} // namespace tudat

#endif // TUDAT_AERODYNAMIC_COEFFICIENT_INTERFACE_H
