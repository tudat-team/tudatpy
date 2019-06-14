/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATERADIATIONPRESSUREINTERFACE_H
#define TUDAT_CREATERADIATIONPRESSUREINTERFACE_H

#include <memory>
#include <boost/make_shared.hpp>

#include "Tudat/SimulationSetup/EnvironmentSetup/body.h"
#include "Tudat/Astrodynamics/ElectroMagnetism/radiationPressureInterface.h"


namespace tudat
{

namespace simulation_setup
{

//! Default values for radiation pressure.
static const std::map< std::string, double > defaultRadiatedPowerValues =
{ { "Sun",  3.839E26 } };

//! List of radiation pressure model types.
enum RadiationPressureType
{
    cannon_ball_radiation_pressure_interface,
    panelled_radiation_pressure_interface
};

//! Base class for radiation pressure interface settings.
/*!
 *  Base class for providing settings for automatic radiation pressure properties creation.  This is
 *  a non-functional base class, specific implementations must be defined in derived classes.
 */
class RadiationPressureInterfaceSettings
{
public:

    //! Constructor.
    /*!
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

    //! Destructor
    virtual ~RadiationPressureInterfaceSettings( ){ }

    //! Function returning type of radiation pressure interface that is to be made.
    /*!
     *  Function returning type of radiation pressure interface that is to be made.
     *  \return Type of radiation pressure interface that is to be made.
     */
    RadiationPressureType getRadiationPressureType( ){ return radiationPressureType_; }

    //! Function returning name of body emitting the radiation.
    /*!
     *  Function returning name of body emitting the radiation.
     *  \return Name of body emitting the radiation.
     */
    std::string getSourceBody( ){ return sourceBody_; }

    //! Function returning list of bodies causing (partial) occultation.
    /*!
     *  Function returning list of bodies causing (partial) occultation.
     *  \return List of bodies causing (partial) occultation.
     */
    std::vector< std::string > getOccultingBodies( ){ return occultingBodies_; }

protected:

    //! Type of radiation pressure interface that is to be made.
    RadiationPressureType radiationPressureType_;

    //! Name of body emitting the radiation.
    std::string sourceBody_;

    //! List of bodies causing (partial) occultation
    std::vector< std::string > occultingBodies_;
};

//! Class providing settings for the creation of a cannonball radiation pressure interface
class CannonBallRadiationPressureInterfaceSettings: public RadiationPressureInterfaceSettings
{
public:

    /*! Constructor
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
        area_( area ), radiationPressureCoefficient_( radiationPressureCoefficient ){ }

    //! Function to return surface area that undergoes radiation pressure.
    /*!
     *  Function to return surface area that undergoes radiation pressure.
     *  \return Surface area that undergoes radiation pressure.
     */
    double getArea( ){ return area_; }

    //! Function to set surface area that undergoes radiation pressure.
    /*!
     *  Function to set surface area that undergoes radiation pressure.
     *  \param area Surface area that undergoes radiation pressure.
     */
    void setArea( double area ){ area_ = area; }

    //! Function to return radiation pressure coefficient.
    /*!
     *  Function to return radiation pressure coefficient.
     *  \return Radiation pressure coefficient.
     */
    double getRadiationPressureCoefficient( ){ return radiationPressureCoefficient_; }


private:

    //! Surface area that undergoes radiation pressure.
    double area_;

    //! Radiation pressure coefficient.
    double radiationPressureCoefficient_;
};

class PanelledRadiationPressureInterfaceSettings: public RadiationPressureInterfaceSettings
{
public:

    /*! Constructor with constant panels orientation in body-fixed frame.
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


    /*! Constructor with time-varying panels orientation in body-fixed frame.
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


    //! Function to return the panels emissivity vector.
    /*!
     *  Function to return the panels emissivity vector.
     *  \return Vector containing the panels emissivities.
     */
    std::vector< double > getEmissivities( )
    {
        return emissivities_;
    }


    //! Function to return the panels area vector.
    /*!
     *  Function to return the panels area vector.
     *  \return Vector containing the panels areas.
     */
    std::vector< double > getAreas( )
    {
        return areas_;
    }

    //! Function to return the panels diffuse reflection coefficients vector.
    /*!
     *  Function to return the panels diffuse reflection coefficients vector.
     *  \return Vector containing the panels diffuse reflection coefficients.
     */
    std::vector< double > getDiffusionCoefficients( )
    {
        return diffusionCoefficients_;
    }

    //! Function to return the panels surface normal vector.
    /*!
     *  Function to return the panels surface normal vector.
     *  \return Vector containing the panels surface normals in body-fixed frame.
     */
    std::vector< std::function< Eigen::Vector3d( const double ) > > getSurfaceNormalsInBodyFixedFrameFunctions( )
    {
        return surfaceNormalsInBodyFixedFrameFunctions_;
    }

private:

    //! Vector containing the emissivitie for all panels.
    std::vector< double > emissivities_;

    //! Vector containing the area for all panels.
    std::vector< double > areas_;

    //! Vector containing the diffuse reflection coefficient for all panels.
    std::vector< double > diffusionCoefficients_;

    //! Vector containing the time-dependent functions that return the surface normal for each panel,
    //! in body-fixed frame.
    std::vector< std::function< Eigen::Vector3d( const double ) > > surfaceNormalsInBodyFixedFrameFunctions_;

};



//! Function to obtain (by reference) the position functions and radii of occulting bodies
/*!
 * Function to obtain (by reference) the position functions and radii of occulting bodies.
 * \param bodyMap List of body objects.
 * \param occultingBodies List of bodies causing occultation.
 * \param occultingBodyPositions List of position functions of occulting bodies (return by reference
 * output variable).
 * \param occultingBodyRadii List of radii of occulting bodies (return by reference
 * output variable).
 */
void getOccultingBodiesInformation(
        const NamedBodyMap& bodyMap, const std::vector< std::string >& occultingBodies,
        std::vector< std::function< Eigen::Vector3d( ) > >& occultingBodyPositions,
        std::vector< double >& occultingBodyRadii );


//! Function to create a radiation pressure interface.
/*!
 *  Function to create a radiation pressure interface.
 *  \param radiationPressureInterfaceSettings Settings for the radiation pressure interface.
 *  \param bodyName Name of body for which radiation pressure interface.
 *  \param bodyMap List of body objects to use for creation of radiation pressure interface.
 *  \return Radiation pressure interface pointer of requested type and settings.
 */
std::shared_ptr< electro_magnetism::RadiationPressureInterface > createRadiationPressureInterface(
        const std::shared_ptr< RadiationPressureInterfaceSettings > radiationPressureInterfaceSettings,
        const std::string& bodyName, const NamedBodyMap& bodyMap );



} // namespace simulation_setup

} // namespace tudat

#endif // TUDAT_CREATERADIATIONPRESSUREINTERFACE_H
