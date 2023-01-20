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
    cannon_ball_radiation_pressure_interface,
    panelled_radiation_pressure_interface,
    solar_sailing_radiation_pressure_interface
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

class PanelledRadiationPressureInterfaceSettings: public RadiationPressureInterfaceSettings
{
public:

    /*  Constructor with constant panels orientation in body-fixed frame.
     * Constructor with constant panels orientation in body-fixed frame.
     * \param sourceBody Name of body emitting the radiation.
     * \param emissivities Vector containing the panels emissivities.
     * \param areas Vector containing the panels areas.
     * \param diffusionCoefficients Vector containing the diffuse reflection coefficients of the panels.
     * \param surfaceNormalsInBodyFixedFrame Vector containing the (constant) surface normals of the panels, expressed in the body-fixed frame.
     * \param occultingBodies List of bodies causing (partial) occultation.
     */
    PanelledRadiationPressureInterfaceSettings(
            const std::string& sourceBody,
            const std::vector< double >& emissivities,
            const std::vector< double >& areas,
            const std::vector< double >& diffusionCoefficients,
            const std::vector< Eigen::Vector3d >& surfaceNormalsInBodyFixedFrame,
            const std::vector< std::string >& occultingBodies = std::vector< std::string >( ) ):
        RadiationPressureInterfaceSettings( panelled_radiation_pressure_interface, sourceBody, occultingBodies ),
        emissivities_( emissivities ), areas_( areas ),  diffusionCoefficients_( diffusionCoefficients )
    {

        for ( unsigned int i = 0 ; i < surfaceNormalsInBodyFixedFrame.size() ; i++ )
        {
            surfaceNormalsInBodyFixedFrameFunctions_.push_back( [ = ]( const double ){ return surfaceNormalsInBodyFixedFrame[ i ]; } );
        }

    }


    /*  Constructor with time-varying panels orientation in body-fixed frame.
     * Constructor with time-varying panels orientation in body-fixed frame.
     * \param sourceBody Name of body emitting the radiation.
     * \param emissivities Vector containing the panels emissivities.
     * \param areas Vector containing the panels areas.
     * \param diffusionCoefficients Vector containing the diffuse reflection coefficients of the panels.
     * \param surfaceNormalsInBodyFixedFrameFunctions Vector containing the functions describing the panels surface normals,
     * which depend on time. The surface normals are expressed in the body-fixed frame.
     * \param occultingBodies List of bodies causing (partial) occultation.
     */
    PanelledRadiationPressureInterfaceSettings(
            const std::string& sourceBody,
            const std::vector< double >& emissivities,
            const std::vector< double >& areas,
            const std::vector< double >& diffusionCoefficients,
            const std::vector< std::function< Eigen::Vector3d( const double ) > >& surfaceNormalsInBodyFixedFrameFunctions,
            const std::vector< std::string >& occultingBodies = std::vector< std::string >( ) ):
        RadiationPressureInterfaceSettings( panelled_radiation_pressure_interface, sourceBody, occultingBodies ),
        emissivities_( emissivities ), areas_( areas ),  diffusionCoefficients_( diffusionCoefficients ),
        surfaceNormalsInBodyFixedFrameFunctions_( surfaceNormalsInBodyFixedFrameFunctions ){ }


    //  Function to return the panels emissivity vector.
    /* 
     *  Function to return the panels emissivity vector.
     *  \return Vector containing the panels emissivities.
     */
    std::vector< double > getEmissivities( )
    {
        return emissivities_;
    }


    //  Function to return the panels area vector.
    /* 
     *  Function to return the panels area vector.
     *  \return Vector containing the panels areas.
     */
    std::vector< double > getAreas( )
    {
        return areas_;
    }

    //  Function to return the panels diffuse reflection coefficients vector.
    /* 
     *  Function to return the panels diffuse reflection coefficients vector.
     *  \return Vector containing the panels diffuse reflection coefficients.
     */
    std::vector< double > getDiffusionCoefficients( )
    {
        return diffusionCoefficients_;
    }

    //  Function to return the panels surface normal vector.
    /* 
     *  Function to return the panels surface normal vector.
     *  \return Vector containing the panels surface normals in body-fixed frame.
     */
    std::vector< std::function< Eigen::Vector3d( const double ) > > getSurfaceNormalsInBodyFixedFrameFunctions( )
    {
        return surfaceNormalsInBodyFixedFrameFunctions_;
    }

private:

    //  Vector containing the emissivitie for all panels.
    std::vector< double > emissivities_;

    //  Vector containing the area for all panels.
    std::vector< double > areas_;

    //  Vector containing the diffuse reflection coefficient for all panels.
    std::vector< double > diffusionCoefficients_;

    //  Vector containing the time-dependent functions that return the surface normal for each panel,
    //  in body-fixed frame.
    std::vector< std::function< Eigen::Vector3d( const double ) > > surfaceNormalsInBodyFixedFrameFunctions_;

};

//! @get_docstring(panelledRadiationPressureInterfaceSettings)
inline std::shared_ptr< RadiationPressureInterfaceSettings > panelledRadiationPressureInterfaceSettings(
		const std::string& sourceBody,
		const std::vector< double >& emissivities,
		const std::vector< double >& areas,
		const std::vector< double >& diffusionCoefficients,
		const std::vector< Eigen::Vector3d >& surfaceNormalsInBodyFixedFrame,
		const std::vector< std::string >& occultingBodies = std::vector< std::string >( )
		)
{
	return std::make_shared< PanelledRadiationPressureInterfaceSettings >(
			sourceBody, emissivities, areas, diffusionCoefficients, surfaceNormalsInBodyFixedFrame,
			occultingBodies );
}

//  Class providing settings for the creation of a solar sail radiation pressure interface.
class SolarSailRadiationInterfaceSettings: public RadiationPressureInterfaceSettings
{
public:

    /*  Constructor
     * Constructor
     * \param sourceBody Name of body emitting the radiation.
     * \param area Surface area that undergoes radiation pressure [m^2].
     * \param coneAngle Function returning the current sail cone angle (in rad).
     * \param clockAngle Function returning the current sail clock angle (in rad).
     * \param frontEmissivityCoefficient Emissivity coefficient of the front of the sail [-].
     * \param backEmissivityCoefficient Emissivity coefficient of the back of the sail [-].
     * \param frontLambertianCoefficient Lambertian coefficient of the front of the sail [-].
     * \param backLambertianCoefficient Lambertian coefficient of the back of the sail [-].
     * \param reflectivityCoefficient Reflectivity coefficient of the sail [-].
     * \param specularReflectionCoefficient Specular reflection coefficient of the sail [-].
     * \param occultingBodies List of bodies causing (partial) occultation.
     * \param centralBody Name of the central body.
     */

    SolarSailRadiationInterfaceSettings(
        const std::string& sourceBody,
        const double area,
        const std::function< double( const double ) > coneAngle,
        const std::function< double( const double ) > clockAngle,
        const double frontEmissivityCoefficient,
        const double backEmissivityCoefficient,
        const double frontLambertianCoefficient,
        const double backLambertianCoefficient,
        const double reflectivityCoefficient,
        const double specularReflectionCoefficient,
        const std::vector< std::string >& occultingBodies = std::vector< std::string >( ),
        const std::string& centralBody = std::string( ) ):
        RadiationPressureInterfaceSettings( solar_sailing_radiation_pressure_interface, sourceBody, occultingBodies ),
        area_( area ), coneAngleFunction_( coneAngle ),
        clockAngleFunction_( clockAngle ), frontEmissivityCoefficient_( frontEmissivityCoefficient ),
        backEmissivityCoefficient_( backEmissivityCoefficient ),
        frontLambertianCoefficient_( frontLambertianCoefficient ),
        backLambertianCoefficient_( backLambertianCoefficient ),
        reflectivityCoefficient_( reflectivityCoefficient ),
        specularReflectionCoefficient_( specularReflectionCoefficient ),
        centralBody_( centralBody ){ }

    //  Function to return surface area that undergoes radiation pressure.
    /* 
     *  Function to return surface area that undergoes radiation pressure.
     *  \return Surface area that undergoes radiation pressure [m^2].
     */
    double getArea( ){ return area_; }

    //  Function to set surface area that undergoes radiation pressure.
    /* 
     *  Function to set surface area that undergoes radiation pressure.
     *  \param area Surface area that undergoes radiation pressure [m^2].
     */
    void setArea( double area ){ area_ = area; }

    //  Returns the function returning the current sail cone angle.
    /* 
     *  Returns the function returning the current sail cone angle.
     *  \return Function returning the current sail cone angle (in rad).
     */
    std::function< double( const double ) > getConeAngle(  ) {return coneAngleFunction_;}

    //  Resets the function returning the current sail cone angle.
    /* 
     *  Resets the function returning the current sail cone angle.
     *  \param coneAngle Function returning the current sail cone angle (in rad).
     */
    void setConeAngle( std::function< double( const double ) > coneAngle ){ coneAngleFunction_ = coneAngle; }

    //  Returns the function returning the current sail clock angle.
    /* 
     *  Returns the function returning the current sail clock angle.
     *  \return Function returning the current clock angle (in rad).
     */
    std::function< double( const double ) > getClockAngle(  ) { return clockAngleFunction_; }

    //  Resets the function returning the current sail clock angle.
    /* 
     *  Resets the function returning the current sail clock angle.
     *  \param clockAngle Function returning the current clock angle (in rad).
     */
    void setClockAngle( std::function< double( const double ) > clockAngle ){ clockAngleFunction_ = clockAngle; }

    //  Returns front emissivity coefficient.
    /* 
     *  Returns front emissivity coefficient.
     *  \return Front emissivity coefficient [-].
     */
    double getFrontEmissivityCoefficient() { return frontEmissivityCoefficient_; }

    //  Resets front emissivity coefficient.
    /* 
     *  Resets front emissivity coefficient.
     *  \param frontEmissivityCoefficient Front emissivity coefficient [-].
     */
    void setFrontEmissivityCoefficient( double frontEmissivityCoefficient ){ frontEmissivityCoefficient_ = frontEmissivityCoefficient; }

    //  Returns back emissivity coefficient.
    /* 
     *  Returns back emissivity coefficient.
     *  \return Back emissivity coefficient [-].
     */
    double getBackEmissivityCoefficient() { return backEmissivityCoefficient_; }

    //  Resets back emissivity coefficient.
    /* 
     *  Resets back emissivity coefficient.
     *  \param backEmissivityCoefficient Back emissivity coefficient [-].
     */
    void setBackEmissivityCoefficient( double backEmissivityCoefficient ){ backEmissivityCoefficient_ = backEmissivityCoefficient; }

    //  Returns front Lambertian coefficient.
    /* 
     *  Returns front Lambertian coefficient.
     *  \return Front Lambertian coefficient [-].
     */
    double getFrontLambertianCoefficient() { return frontLambertianCoefficient_; }

    //  Resets front Lambertian coefficient.
    /* 
     *  Resets front Lambertian coefficient.
     *  \param frontLambertianCoefficient Front Lambertian coefficient [-].
     */
    void setFrontLambertianCoefficient( double frontLambertianCoefficient ){ frontLambertianCoefficient_ = frontLambertianCoefficient; }

    //  Returns back Lambertian coefficient.
    /* 
     *  Returns back Lambertian coefficient.
     *  \return Back Lambertian coefficient [-].
     */
    double getBackLambertianCoefficient() { return backLambertianCoefficient_; }

    //  Resets back Lambertian coefficient.
    /* 
     *  Resets back Lambertian coefficient.
     *  \param backLambertianCoefficient Back Lambertian coefficient [-].
     */
    void setBackLambertianCoefficient( double backLambertianCoefficient ){ backLambertianCoefficient_ = backLambertianCoefficient; }

    //  Returns reflectivity coefficient.
    /* 
     *  Returns reflectivity coefficient.
     *  \return Reflectivity coefficient [-].
     */
    double getReflectivityCoefficient() { return reflectivityCoefficient_; }

    //  Resets reflectivity coefficient.
    /* 
     *  Resets reflectivity coefficient.
     *  \param reflectivityCoefficient Reflectivity coefficient [-].
     */
    void setReflectivityCoefficient( double reflectivityCoefficient ){ reflectivityCoefficient_ = reflectivityCoefficient; }

    //  Returns specular reflection coefficient.
    /* 
     *  Returns specular reflection coefficient.
     *  \return Specular reflection coefficient [-].
     */
    double getSpecularReflectionCoefficient() { return specularReflectionCoefficient_; }

    //  Resets specular reflection coefficient.
    /* 
     *  Resets specular reflection coefficient.
     *  \param specularReflectionCoefficient Specular reflection coefficient [-].
     */
    void setSpecularReflectionCoefficient( double specularReflectionCoefficient ){ specularReflectionCoefficient_ = specularReflectionCoefficient; }

    //  Returns the name of the central body.
    /* 
     *  Returns the name of the central body.
     *  \return Name of the central body.
     */
    std::string getCentralBody( ){ return centralBody_; }


private:

    //  Surface area that undergoes radiation pressure [m^2].
    double area_;

    //  Cone angle of the sail [rad].
    std::function< double( const double ) > coneAngleFunction_;

    //  Clock angle of the sail [rad].
    std::function< double( const double ) > clockAngleFunction_;

    //  Front emissivity coefficient of the sail [-].
    double frontEmissivityCoefficient_;

    //  Back emissivity coefficient of the sail [-].
    double backEmissivityCoefficient_;

    //  Front Lambertian coefficient of the sail [-].
    double frontLambertianCoefficient_;

    //  Back lambertian coefficient of the sail [-].
    double backLambertianCoefficient_;

    //  Reflectivity coefficient of the sail [-].
    double reflectivityCoefficient_;

    //  Specular reflection coefficient of the sail [-].
    double specularReflectionCoefficient_;

    //  Name of the central body.
    std::string centralBody_;

};


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
