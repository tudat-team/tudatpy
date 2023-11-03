/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
//RP-OLD

#ifndef TUDAT_CREATERADIATIONPRESSUREINTERFACE_H
#define TUDAT_CREATERADIATIONPRESSUREINTERFACE_H

#include <memory>
#include <functional>


#include "tudat/simulation/environment_setup/body.h"
#include "tudat/astro/electromagnetism/radiationPressureInterface.h"


namespace tudat
{

namespace simulation_setup
{

//  Default values for radiation pressure.
static const std::map< std::string, double > defaultRadiatedPowerValues =
{ { "Sun",  3.839E26 } };

//  List of radiation pressure model types.
//! @get_docstring(RadiationPressureType.__docstring__)
enum RadiationPressureType
{
    cannon_ball_radiation_pressure_interface
//    panelled_radiation_pressure_interface,
//    solar_sailing_radiation_pressure_interface
};

//  Base class for radiation pressure interface settings.
/* 
 *  Base class for providing settings for automatic radiation pressure properties creation.  This is
 *  a non-functional base class, specific implementations must be defined in derived classes.
 */
//! @get_docstring(RadiationPressureInterfaceSettings.__docstring__)
class RadiationPressureInterfaceSettings
{
public:

    //  Constructor.
    /* 
     * Constructor.
     * \param radiationPressureType Type of radiation pressure interface that is to be made.
     * \param sourceBody Name of body emitting the radiation.
     * \param occultingBodies List of bodies causing (partial) occultation
     */
    RadiationPressureInterfaceSettings(
            const RadiationPressureType radiationPressureType,
            const std::string& sourceBody,
            const std::vector< std::string > occultingBodies = std::vector< std::string >( ) ):
        radiationPressureType_( radiationPressureType ), sourceBody_( sourceBody ),
        occultingBodies_( occultingBodies ){  }

    //  Destructor
    virtual ~RadiationPressureInterfaceSettings( ){ }

    //  Function returning type of radiation pressure interface that is to be made.
    /* 
     *  Function returning type of radiation pressure interface that is to be made.
     *  \return Type of radiation pressure interface that is to be made.
     */
    RadiationPressureType getRadiationPressureType( ){ return radiationPressureType_; }

    //  Function returning name of body emitting the radiation.
    /* 
     *  Function returning name of body emitting the radiation.
     *  \return Name of body emitting the radiation.
     */
    std::string getSourceBody( ){ return sourceBody_; }

    //  Function returning list of bodies causing (partial) occultation.
    /* 
     *  Function returning list of bodies causing (partial) occultation.
     *  \return List of bodies causing (partial) occultation.
     */
    std::vector< std::string > getOccultingBodies( ){ return occultingBodies_; }

protected:

    //  Type of radiation pressure interface that is to be made.
    RadiationPressureType radiationPressureType_;

    //  Name of body emitting the radiation.
    std::string sourceBody_;

    //  List of bodies causing (partial) occultation
    std::vector< std::string > occultingBodies_;
};

//  Class providing settings for the creation of a cannonball radiation pressure interface
//! @get_docstring(CannonBallRadiationPressureInterfaceSettings.__docstring__)
class CannonBallRadiationPressureInterfaceSettings: public RadiationPressureInterfaceSettings
{
public:

    /*  Constructor
     * Constructor
     * \param sourceBody Name of body emitting the radiation.
     * \param area Surface area that undergoes radiation pressure.
     * \param radiationPressureCoefficient Radiation pressure coefficient.
     * \param occultingBodies List of bodies causing (partial) occultation.
     */
    CannonBallRadiationPressureInterfaceSettings(
            const std::string& sourceBody, const double area, const double radiationPressureCoefficient,
            const std::vector< std::string >& occultingBodies = std::vector< std::string >( ) ):
        RadiationPressureInterfaceSettings( cannon_ball_radiation_pressure_interface, sourceBody, occultingBodies ),
        area_( area ), radiationPressureCoefficient_( radiationPressureCoefficient ),
        radiationPressureCoefficientFunction_( [=]( const double ){ return radiationPressureCoefficient;} ){ }

    CannonBallRadiationPressureInterfaceSettings(
            const std::string& sourceBody, const double area, const std::function< double(  const double ) > radiationPressureCoefficientFunction,
            const std::vector< std::string >& occultingBodies = std::vector< std::string >( ) ):
        RadiationPressureInterfaceSettings( cannon_ball_radiation_pressure_interface, sourceBody, occultingBodies ),
        area_( area ), radiationPressureCoefficient_( TUDAT_NAN ),
        radiationPressureCoefficientFunction_( radiationPressureCoefficientFunction ){ }

    //  Function to return surface area that undergoes radiation pressure.
    /* 
     *  Function to return surface area that undergoes radiation pressure.
     *  \return Surface area that undergoes radiation pressure.
     */
    double getArea( ){ return area_; }

    //  Function to set surface area that undergoes radiation pressure.
    /* 
     *  Function to set surface area that undergoes radiation pressure.
     *  \param area Surface area that undergoes radiation pressure.
     */
    void setArea( double area ){ area_ = area; }

    //  Function to return radiation pressure coefficient.
    /* 
     *  Function to return radiation pressure coefficient.
     *  \return Radiation pressure coefficient.
     */
    double getRadiationPressureCoefficient( ){ return radiationPressureCoefficient_; }

    std::function< double( const double ) >  getRadiationPressureCoefficientFunction( ){ return radiationPressureCoefficientFunction_; }


private:

    //  Surface area that undergoes radiation pressure.
    double area_;

    //  Radiation pressure coefficient.
    double radiationPressureCoefficient_;

    std::function< double( const double ) > radiationPressureCoefficientFunction_;
};
//! @get_docstring(cannonBallRadiationPressureSettings)
inline std::shared_ptr< RadiationPressureInterfaceSettings > cannonBallRadiationPressureSettings(
        const std::string& sourceBody, const double area, const double radiationPressureCoefficient,
        const std::vector< std::string >& occultingBodies )
{
    return std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                sourceBody, area, radiationPressureCoefficient, occultingBodies );
}

inline std::shared_ptr< RadiationPressureInterfaceSettings > cannonBallRadiationPressureSettings(
        const std::string& sourceBody, const double area, const std::function< double(  const double ) >  radiationPressureCoefficientFunction,
        const std::vector< std::string >& occultingBodies )
{
    return std::make_shared< CannonBallRadiationPressureInterfaceSettings >(
                sourceBody, area, radiationPressureCoefficientFunction, occultingBodies );
}


//  Function to obtain (by reference) the position functions and radii of occulting bodies
/* 
 * Function to obtain (by reference) the position functions and radii of occulting bodies.
 * \param bodies List of body objects.
 * \param occultingBodies List of bodies causing occultation.
 * \param occultingBodyPositions List of position functions of occulting bodies (return by reference
 * output variable).
 * \param occultingBodyRadii List of radii of occulting bodies (return by reference
 * output variable).
 */
void getOccultingBodiesInformation(
        const SystemOfBodies& bodies, const std::vector< std::string >& occultingBodies,
        std::vector< std::function< Eigen::Vector3d( ) > >& occultingBodyPositions,
        std::vector< double >& occultingBodyRadii );


//  Function to obtain (by reference) the position functions and velocity of the central body.
/* 
 * Function to obtain (by reference) the position functions and velocity of the central body.
 * \param bodies List of body objects.
 * \param centralBody Name of the central body.
 * \param centralBodyPosition Central body's position function (return by reference output variable).
 * \param centralBodyVelocity Central body's velocity function (return by reference output variable).
 */
void getCentralBodyInformation(
    const SystemOfBodies& bodies, const std::string& centralBody,
    std::function< Eigen::Vector3d( ) >& centralBodyPosition,
    std::function< Eigen::Vector3d( ) >& centralBodyVelocity);

//  Function to create a radiation pressure interface.
/* 
 *  Function to create a radiation pressure interface.
 *  \param radiationPressureInterfaceSettings Settings for the radiation pressure interface.
 *  \param bodyName Name of body for which radiation pressure interface.
 *  \param bodies List of body objects to use for creation of radiation pressure interface.
 *  \return Radiation pressure interface pointer of requested type and settings.
 */
std::shared_ptr< electromagnetism::RadiationPressureInterface > createRadiationPressureInterface(
        const std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureInterfaceSettings,
        const std::string& bodyName, const SystemOfBodies& bodies );

std::function< double( const double ) > getOccultationFunction(
        const SystemOfBodies& bodyMap,
        const std::string& sourceBody,
        const std::string& occultingBody,
        const std::string& shadowedBody );

} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATERADIATIONPRESSUREINTERFACE_H
