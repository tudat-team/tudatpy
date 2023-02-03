/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 *    References
 *      E. B. Saff, et al. "Distributing many points on a sphere".
 *          The Mathematical Intelligencer 19. 1(1997): 5–11.
 *      Frank G. Lemoine, et al. "High‒degree gravity models from GRAIL primary mission data".
 *          Journal of Geophysical Research: Planets 118. 8(2013): 1676–1698.
 *      Knocke, Philip et al. "Earth radiation pressure effects on satellites."
 *          Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988.
 *      Charles J. Wetterer, et al. "Refining Space Object Radiation Pressure Modeling with Bidirectional Reflectance
 *          Distribution Functions". Journal of Guidance, Control, and Dynamics 37. 1(2014): 185–196.
 */

#ifndef TUDAT_RADIATIONSOURCEMODEL_H
#define TUDAT_RADIATIONSOURCEMODEL_H

#include <utility>
#include <vector>
#include <map>
#include <utility>
#include <memory>

#include <Eigen/Core>
#include <Eigen/Geometry>

#include "tudat/astro/electromagnetism/luminosityModel.h"
#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/astro/electromagnetism/surfacePropertyDistribution.h"
#include "tudat/astro/basic_astro/bodyShapeModel.h"
#include "tudat/math/basic/mathematicalConstants.h"
#include "tudat/math/basic/coordinateConversions.h"


namespace tudat
{
namespace electromagnetism
{

typedef std::vector<std::pair<double, Eigen::Vector3d>> IrradianceWithSourceList;

/*!
 * Class modeling a body that emits electromagnetic radiation, used for radiation pressure calculations.
 *
 * A source may consist of multiple sub-sources. For example, in a paneled source, each panel acts as a sub-source, with
 * the radiation from each sub-source arriving at the target from a (slightly) different direction.
 *
 * Parameters to all functions are in local source-fixed coordinates. Proper translation and rotation is ensured by the
 * RadiationPressureAcceleration class.
 */
class RadiationSourceModel
{
public:
    explicit RadiationSourceModel() = default;

    virtual ~RadiationSourceModel() = default;

    /*!
     * Update class members.
     *
     * @param currentTime Current simulation time
     */
    void updateMembers(double currentTime);

    /*!
     * Evaluate the irradiance [W/m²] at a certain position due to this source.
     *
     * @param targetPosition Position where to evaluate the irradiance in local (i.e. source-fixed) coordinates
     * @param originalSourceIrradiance Irradiance from the original source (if applicable)
     * @param originalSourceToSourceDirection Direction of incoming radiation in local (i.e. source-fixed) coordinates
     * @return List of irradiances at the target position and their source-fixed origin due to all sub-sources.
     *         Contains a single element for point sources, multiple elements for paneled sources. Each element can be
     *         thought of as ray from the source to the target.
     */
     //TODO-DOMINIK move original source params to subclass
    virtual IrradianceWithSourceList evaluateIrradianceAtPosition(
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const = 0;

    /*!
     * Evaluate the total irradiance [W/m²] at a certain position due to this source. Individual rays are merged, as
     * opposed to evaluateIrradianceAtPosition().
     *
     * @param targetPosition Position where to evaluate the irradiance in local (i.e. source-fixed) coordinates
     * @param originalSourceIrradiance Irradiance from the original source (if applicable)
     * @param originalSourceToSourceDirection Direction of incoming radiation in local (i.e. source-fixed) coordinates
     * @return Total irradiance from this source at the target position [W/m²]
     */
    double evaluateTotalIrradianceAtPosition(
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const;

protected:
    virtual void updateMembers_(const double currentTime) {};

    double currentTime_{TUDAT_NAN};
};

//*********************************************************************************************
//   Isotropic point radiation source
//*********************************************************************************************

/*!
 * Class modeling a radiation source that emits isotropically in all directions from a single point. The source power
 * is given by the total luminosity.
 *
 * Far-away sources such as the Sun at 1 AU can be modeled as isotropic point source.
 */
class IsotropicPointRadiationSourceModel : public RadiationSourceModel
{
public:
    /*!
     * Constructor.
     *
     * @param luminosityModel Luminosity model of this source
     */
    explicit IsotropicPointRadiationSourceModel(
            const std::shared_ptr<LuminosityModel>& luminosityModel) :
        luminosityModel_(luminosityModel) {}

    IrradianceWithSourceList evaluateIrradianceAtPosition(
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    /*!
     * Evaluate the irradiance [W/m²] at a certain position due to this source.
     *
     * @param targetPosition Position where to evaluate the irradiance in local (i.e. source-fixed) coordinates
     * @return Irradiance from this source at the target position [W/m²]
     */
    double evaluateIrradianceAtPosition(const Eigen::Vector3d& targetPosition) const;

    const std::shared_ptr<LuminosityModel>& getLuminosityModel() const
    {
        return luminosityModel_;
    }

private:
    void updateMembers_(double currentTime) override;

    std::shared_ptr<LuminosityModel> luminosityModel_;
};

//*********************************************************************************************
//   Paneled radiation source
//*********************************************************************************************

/*!F
 * Class modeling a radiation source that is discretized into panels, each with its own radiosity model. A radiosity
 * model describes the radiation emitted (thermal) or reflected (albedo) by each panel. Since both depend on the
 * radiation from a third body (referred to as the original source), all irradiance evaluations involve the irradiance
 * from the original source (for now, always an isotropic point source) at the paneled source center.
 *
 * Sources close to the spacecraft with the radiosity depending on a third body, such as Earth and Moon reflecting
 * solar radiation due to surface albedo and thermal reradiation, can be modeled as a paneled source. For farther away
 * sources, all radiation can be treated as originating from the same point, which is modeles by a point radiation
 * source.
 */
class PaneledRadiationSourceModel : public RadiationSourceModel
{
public:
    class Panel;
    class PanelRadiosityModel; // cannot make this nested class of Panel since forward declaration is impossible

    /*!
     * Constructor for a shape-aware paneled source. The shape model is necessary for automatic panel generation.
     *
     * @param originalSourceName Name of the original source body
     * @param sourceBodyShapeModel Shape model of this source
     * @param originalSourceToSourceOccultingBodies Names of bodies to occult the original source as seen from this source
     */
    explicit PaneledRadiationSourceModel(
            const std::string& originalSourceName,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& sourceBodyShapeModel,
            const std::vector<std::string>& originalSourceToSourceOccultingBodies) :
        originalSourceName_(originalSourceName),
        sourceBodyShapeModel_(sourceBodyShapeModel),
        originalSourceToSourceOccultingBodies_(originalSourceToSourceOccultingBodies) {}

    /*!
     * Constructor for shape-oblivious paneled source. No knowledge of shape necessary if panels are given manually.
     *
     * @param originalSourceName Name of the original source body
     */
    explicit PaneledRadiationSourceModel(
        const std::string& originalSourceName,
        const std::vector<std::string>& originalSourceToSourceOccultingBodies = {}) :
            PaneledRadiationSourceModel(originalSourceName, nullptr, originalSourceToSourceOccultingBodies) {}

    IrradianceWithSourceList evaluateIrradianceAtPosition(
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    /*!
     * Get all panels comprising this paneled source.
     *
     * @return Panels comprising this paneled source
     */
    virtual const std::vector<Panel>& getPanels() const = 0;

    const std::string& getOriginalSourceName() const
    {
        return originalSourceName_;
    }

    std::vector<std::string> getOriginalSourceToSourceOccultingBodies() const
    {
        return originalSourceToSourceOccultingBodies_;
    }

    /*!
    * Get the number of panels comprising this paneled source. For dynamically generated paneling, this may be
    * the target number, not the currently existing one.
    */
    virtual unsigned int getNumberOfPanels() const = 0;

protected:
    std::string originalSourceName_; // needed for environment updater setup

    std::shared_ptr<basic_astrodynamics::BodyShapeModel> sourceBodyShapeModel_;

    // Only needed to transfer occultation settings from body setup to acceleration setup
    std::vector<std::string> originalSourceToSourceOccultingBodies_;
};

/*!
 * Class modeling a paneled source with paneling of the whole body that is constant and generated only once. Panel
 * properties such as albedo and emissivity are constant in time as well. The panels are generated in the first call
 * of updateMembers() using the staggered method.
 *
 * Static paneling is easier to implement than dynamic paneling, and possibly costly albedo/emissivity evaluations (e.g.,
 * if their distributions are given by spherical harmonics) only have to be done once. However, a large portion of the
 * panels will not be visible, especially for spacecraft in low orbits. Furthermore, for sufficient panels to be visible
 * for low spacecraft, a large total number of panels is necessary. Both issues increase computational effort
 * unnecessarily.
 */
class StaticallyPaneledRadiationSourceModel : public PaneledRadiationSourceModel
{
public:
    /*!
     * Constructor for given panels.
     *
     * @param originalSourceName Name of the original source body
     * @param panels Panels comprising this paneled source
     * @param originalSourceToSourceOccultingBodies Names of bodies to occult the original source as seen from this source
     */
    explicit StaticallyPaneledRadiationSourceModel(
            const std::string& originalSourceName,
            std::vector<Panel> panels,
            const std::vector<std::string>& originalSourceToSourceOccultingBodies = {}) :
            PaneledRadiationSourceModel(originalSourceName, originalSourceToSourceOccultingBodies),
            numberOfPanels(panels.size()),
            panels_(std::move(panels)) {}

    /*!
     * Constructor for automatic generation of panels.
     *
     * @param originalSourceName Name of the original source body
     * @param sourceBodyShapeModel Shape model of this source
     * @param baseRadiosityModels Radiosity models that will be copied for each panel
     * @param numberOfPanels Number of panels for automatic generation
     * @param originalSourceToSourceOccultingBodies Names of bodies to occult the original source as seen from this source
     */
    explicit StaticallyPaneledRadiationSourceModel(
            const std::string& originalSourceName,
            const std::shared_ptr<basic_astrodynamics::BodyShapeModel>& sourceBodyShapeModel,
            const std::vector<std::unique_ptr<PanelRadiosityModel>>& baseRadiosityModels,
            int numberOfPanels,
            const std::vector<std::string>& originalSourceToSourceOccultingBodies = {}) :
            PaneledRadiationSourceModel(
                originalSourceName, sourceBodyShapeModel, originalSourceToSourceOccultingBodies),
            numberOfPanels(numberOfPanels)
    {
        generatePanels(baseRadiosityModels);
    }

    const std::vector<Panel>& getPanels() const override
    {
        return panels_;
    }

    unsigned int getNumberOfPanels() const override
    {
        return numberOfPanels;
    }

private:
    void updateMembers_(double currentTime) override;

    void generatePanels(const std::vector<std::unique_ptr<PanelRadiosityModel>>& baseRadiosityModels);


    unsigned int numberOfPanels;

    std::vector<Panel> panels_;
};

///*!
// * Class modeling a paneled source with paneling created dynamically based on the spacecraft position. Panel
// * properties such as albedo and emissivity are evaluated anew with every panel generation.
// *
// * This is the classic paneling for albedo modeling introduced by Knocke (1998). Panels are only generated for a
// * spherical cap under the spacecraft, divided into 2 rings and a central panel.
// */
//// TODO-DOMINIK figure out how to implement dynamic paneling with support for multiple targets
//class DynamicallyPaneledRadiationSourceModel : public PaneledRadiationSourceModel
//{
//private:
//    // keep target-specific panel lists, otherwise may have to regenerate a lot
//    std::map<std::string, std::vector<Panel>> panels_;
//};

/*!
 * Class modeling a single panel on a paneled source.
 *
 * A panel can be given manually or automatically generated, and only contains geometric properties. All
 * radiation-related functionality is delegated to radiosity models. The irradiance at a position due to one panel is
 * the sum of the contributions of all radiosity models of that panel. A panel can have a single radiosity model (e.g.,
 * only albedo) or multiple (albedo and thermal).
 */
class PaneledRadiationSourceModel::Panel
{
public:
    EIGEN_MAKE_ALIGNED_OPERATOR_NEW

    /*!
     * Constructor.
     *
     * @param area Area of the panel [m²]
     * @param relativeCenter Center of the panel relative to the source center, in source-fixed Cartesian coordinates [m]
     * @param surfaceNormal Surface normal vector of the panel [-]
     * @param radiosityModels List of radiosity models of the panel
     */
    explicit Panel(
            double area,
            const Eigen::Vector3d& relativeCenter,
            const Eigen::Vector3d& surfaceNormal,
            std::vector<std::unique_ptr<PanelRadiosityModel>> radiosityModels) :
            area_(area),
            relativeCenter_(relativeCenter),
            surfaceNormal_(surfaceNormal),
            radiosityModels_(std::move(radiosityModels))
    {
        Eigen::Vector3d relativeCenterInSphericalCoords =
                coordinate_conversions::convertCartesianToSpherical(relativeCenter);
        latitude_ = M_PI_2 - relativeCenterInSphericalCoords[1];
        longitude_ = relativeCenterInSphericalCoords[2];
    }

    /*!
    * Update class members.
    *
    * @param currentTime Current simulation time
    */
    void updateMembers(double currentTime);

    double getArea() const
    {
        return area_;
    }

    double getLatitude() const
    {
        return latitude_;
    }

    double getLongitude() const
    {
        return longitude_;
    }

    const Eigen::Vector3d& getRelativeCenter() const
    {
        return relativeCenter_;
    }

    const Eigen::Vector3d& getSurfaceNormal() const
    {
        return surfaceNormal_;
    }

    const std::vector<std::unique_ptr<PanelRadiosityModel>>& getRadiosityModels() const
    {
        return radiosityModels_;
    }

private:
    double area_;
    double latitude_;
    double longitude_;
    Eigen::Vector3d relativeCenter_;
    Eigen::Vector3d surfaceNormal_;
    std::vector<std::unique_ptr<PanelRadiosityModel>> radiosityModels_;
};

/*!
 * Class modeling the radiosity emitted by a panel of a paneled radiation source.
 *
 * @see PaneledRadiationSourceModel::Panel
 */
class PaneledRadiationSourceModel::PanelRadiosityModel
{
public:
    virtual ~PanelRadiosityModel() = default;

    /*!
     * Evaluate the irradiance [W/m²] at a certain position due to this panel.
     *
     * @param targetPosition Position where to evaluate the irradiance in panel-local coordinates (source rotation,
     *        centered in panel)
     * @return Irradiance due to this radiosity model for single panel
     */
    virtual double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const = 0;

    /*!
     * Update class members.
     *
     * @param panelLatitude Latitude of the panel this radiosity model belongs to
     * @param panelLongitude Longitude of the panel this radiosity model belongs to
     * @param currentTime Current simulation time
     */
    void updateMembers(
            const double panelLatitude,
            const double panelLongitude,
            const double currentTime);

    /*!
     * Clone this object.
     *
     * @return A clone of this object
     */
    virtual std::unique_ptr<PanelRadiosityModel> clone() const = 0;

protected:
    virtual void updateMembers_(
            const double panelLatitude,
            const double panelLongitude,
            const double currentTime) {};

    /*!
     * Whether the radiosity model is invariant with time. If yes, its members will not be updated even if the time
     * changed, it will only be updated when the latitude/longitude change. The time invariance is usually determined
     * by the albedo/emissivity distribution.
     *
     * @return Whether the radiosity model is invariant with time
     */
    virtual bool isTimeInvariant() = 0;

    double currentTime_{TUDAT_NAN};
    double panelLatitude_{TUDAT_NAN};
    double panelLongitude_{TUDAT_NAN};
};

/*!
 * Panel radiosity model for albedo radiation with an arbitrary reflection law. This model was introduced in
 * Knocke (1988) for Earth thermal radiation, but only assuming Lambertian reflectance.
 *
 * For most cases, albedo radiation with a diffuse-only Lambertian reflection law is sufficient. Only a small fraction
 * of Earth's albedo is actually specular, since only calm bodies of water are truly specular and specular reflection
 * occurs mostly at low solar zenith angles where the reflectance is low anyway (Knocke, 1988).
 *
 * More sophisticated reflection laws (e.g., SpecularDiffuseMixReflectionLaw, or any BRDF via a custom ReflectionLaw)
 * can be implemented to take into account surface properties like different vegetation and ground types on Earth.
 */
class AlbedoPanelRadiosityModel : public PaneledRadiationSourceModel::PanelRadiosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param albedoDistribution Albedo distribution
     * @param withInstantaneousReradiation Whether to instantaneously reradiate absorbed radiation
     */
    explicit AlbedoPanelRadiosityModel(
            const std::shared_ptr<SurfacePropertyDistribution>& albedoDistribution,
            bool withInstantaneousReradiation = false) :
            albedoDistribution_(albedoDistribution),
            reflectionLaw_(std::make_shared<LambertianReflectionLaw>(TUDAT_NAN, withInstantaneousReradiation)) {}

    double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    std::unique_ptr<PanelRadiosityModel> clone() const override
    {
        return std::make_unique<AlbedoPanelRadiosityModel>(*this);
    }

    const std::shared_ptr<LambertianReflectionLaw>& getReflectionLaw() const
    {
        return reflectionLaw_;
    }

private:
    void updateMembers_(
            double panelLatitude,
            double panelLongitude,
            double currentTime) override;

    bool isTimeInvariant() override
    {
        return albedoDistribution_->isTimeInvariant();
    }

    std::shared_ptr<SurfacePropertyDistribution> albedoDistribution_;

    // Reflection law governing reflection of original source radiation
    std::shared_ptr<LambertianReflectionLaw> reflectionLaw_;
};

/*!
 * Panel radiosity model for thermal emissions based on delayed, isotropic and constant flux. This model was introduced
 * in Knocke (1988) for Earth thermal radiation.
 *
 * As opposed to instantaneous reradiation, the body is assumed to act as heat buffer, absorbing incoming radiation and
 * reradiating it in a delayed fashion as longwave infrared radiation. For most bodies, and especially Earth, this is a
 * good assumption (Knocke, 1988).
 */
class DelayedThermalPanelRadiosityModel : public PaneledRadiationSourceModel::PanelRadiosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param emissivityDistribution Emissivity distribution
     */
    explicit DelayedThermalPanelRadiosityModel(
            const std::shared_ptr<SurfacePropertyDistribution>& emissivityDistribution) :
            emissivityDistribution_(emissivityDistribution) {}

    double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;

    std::unique_ptr<PanelRadiosityModel> clone() const override
    {
        return std::make_unique<DelayedThermalPanelRadiosityModel>(*this);
    }

    double getEmissivity() const
    {
        return emissivity;
    }

private:
    void updateMembers_(
            double panelLatitude,
            double panelLongitude,
            double currentTime) override;

    bool isTimeInvariant() override
    {
        return emissivityDistribution_->isTimeInvariant();
    }

    std::shared_ptr<SurfacePropertyDistribution> emissivityDistribution_;
    double emissivity{TUDAT_NAN};
};

/*!
 * Panel radiosity model for thermal emissions based on the angle to sub-solar point and emissivity-corrected
 * black-body radiation. This model was introduced in Lemoine (2013) for lunar thermal radiation with
 * minTemperature = 100 K and maxTemperature = 375 K.
 *
 * The surface temperature is approximated by interpolation between minTemperature and maxTemperature based on the angle
 * to the subsolar point. From the surface temperature, the black-body radiation is then calculated and corrected for
 * emissivity. The radiation is maximum if the target is above the sub-solar point. At positions away from the sub-solar
 * point, the temperature drops, until it reaches minTemperature on the backside.
 */
class AngleBasedThermalPanelRadiosityModel : public PaneledRadiationSourceModel::PanelRadiosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param minTemperature Minimum surface temperature (in shade) [K]
     * @param maxTemperature Maximum surface temperature (at subsolar point) [K]
     * @param emissivityDistribution Emissivity distribution
     */
    explicit AngleBasedThermalPanelRadiosityModel(
            double minTemperature,
            double maxTemperature,
            const std::shared_ptr<SurfacePropertyDistribution>& emissivityDistribution) :
            minTemperature_(minTemperature),
            maxTemperature_(maxTemperature),
            emissivityDistribution_(emissivityDistribution) {}

    double evaluateIrradianceAtPosition(
            double panelArea,
            const Eigen::Vector3d& panelSurfaceNormal,
            const Eigen::Vector3d& targetPosition,
            double originalSourceIrradiance,
            const Eigen::Vector3d& originalSourceToSourceDirection) const override;


    std::unique_ptr<PanelRadiosityModel> clone() const override
    {
        return std::make_unique<AngleBasedThermalPanelRadiosityModel>(*this);
    }

    double getMinTemperature() const
    {
        return minTemperature_;
    }

    double getMaxTemperature() const
    {
        return maxTemperature_;
    }

    double getEmissivity() const
    {
        return emissivity;
    }

private:
    void updateMembers_(
            double panelLatitude,
            double panelLongitude,
            double currentTime) override;

    bool isTimeInvariant() override
    {
        return emissivityDistribution_->isTimeInvariant();
    }

    double minTemperature_;
    double maxTemperature_;
    std::shared_ptr<SurfacePropertyDistribution> emissivityDistribution_;
    double emissivity{TUDAT_NAN};
};

/*!
 * Generate evenly spaced points on a sphere using algorithm from Saff (1997). This is used for uniform panel generation
 * for the statically paneled source. Points are given from southmost to northmost. The poles are included.
 *
 * The generated points are more staggered with Wetterer's algorithm, while Saff's points are more aligned along zonal
 * and meridional lines (spiraling). Staggered points should be preferred.
 *
 * @param n number of points to generate
 * @return a pair of vectors, first vector are polar angles, second vector are azimuth angles (both in radians)
 */
std::pair<std::vector<double>, std::vector<double>> generateEvenlySpacedPoints_Spiraling(unsigned int n);

/*!
 * Generate evenly spaced points on a sphere using algorithm from Wetterer (2014). This is used for uniform panel
 * generation for the statically paneled source. Points are given from southmost to northmost. The poles are not included.
 *
 * The generated points are more staggered with Wetterer's algorithm, while Saff's points are more aligned along zonal
 * and meridional lines (spiraling). Staggered points should be preferred.
 *
 * @param n number of points to generate
 * @return a pair of vectors, first vector are polar angles, second vector are azimuth angles (both in radians)
 */
std::pair<std::vector<double>, std::vector<double>> generateEvenlySpacedPoints_Staggered(unsigned int n);


} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONSOURCEMODEL_H
