/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createRadiationSourceModel.h"

#include <memory>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/basic_astro/physicalConstants.h"


namespace tudat
{
namespace simulation_setup
{

SphericalHarmonicsSurfacePropertyDistributionSettings::SphericalHarmonicsSurfacePropertyDistributionSettings(
        SphericalHarmonicsSurfacePropertyDistributionModel model) :
        SphericalHarmonicsSurfacePropertyDistributionSettings(Eigen::MatrixXd(), Eigen::MatrixXd())
{
    model_ = model;

    switch(model)
    {
        case SphericalHarmonicsSurfacePropertyDistributionModel::albedo_dlam1:
        {
            // DLAM-1 lunar albedo model: Floberghagen, R. et al. "Lunar Albedo Force Modeling and its Effect on Low Lunar Orbit and Gravity Field Determination". ASR 23. 4(1999): 733-738.
            auto maximumDegree = 15;
            auto maximumOrder = 15;

            cosineCoefficients_ = Eigen::MatrixXd (maximumDegree + 1, maximumOrder + 1);
            cosineCoefficients_ << 0.19467246, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0037619635, -0.036562523, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.026013759, -0.023342472, -0.0015511737, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.010791359, 0.012984099, -0.00051947501, 0.00083037838, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0015490099, 0.0024150199, -0.00045980599, -5.8965344e-05, 0.00010396812, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    -0.0047391417, 0.00097388352, -0.00066144592, 3.629954e-05, 3.9759108e-06, 2.270628e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0021741056, 0.0019854158, 3.6293963e-05, -2.763806e-05, -4.6429902e-06, -3.8219121e-06, 3.8641269e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    -0.0026214758, 0.00096586546, -7.8726437e-05, 3.5085491e-05, 3.0497316e-06, -4.0231517e-08, -5.2535643e-08, -1.3545976e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    -0.0072392123, -3.6002321e-05, -0.00015504282, -4.2175334e-06, 2.0927038e-06, -2.1213522e-07, 5.7936412e-08, -2.8380071e-09, 1.0183988e-09, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    -0.004653015, 0.00077363813, 5.8165535e-05, -2.2643899e-05, -3.9981338e-07, -2.2805609e-08, -2.1918927e-08, 4.0524183e-09, 5.7278542e-10, 1.016917e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    -0.015925583, -1.4417769e-05, 5.8118004e-05, -1.0071993e-06, 3.8769952e-07, -1.2667752e-08, 8.7088087e-09, 1.0255584e-09, 3.500616e-10, 3.1475729e-11, 2.3275218e-12, 0.0, 0.0, 0.0, 0.0, 0.0,
                    -0.0037225793, -0.00041507339, -6.3009829e-05, -2.5513158e-06, 1.5175304e-07, -5.7092276e-08, 1.9616945e-09, 3.1055117e-10, 3.6901893e-11, 5.4058712e-12, 6.7325209e-14, 8.2735779e-14, 0.0, 0.0, 0.0, 0.0,
                    -0.0077462959, 0.00068495466, 0.00014207584, -1.0655919e-06, -3.481692e-07, 5.7062063e-09, 5.0926848e-11, 9.6192077e-11, 6.7701005e-12, 2.1793964e-13, -1.3931937e-13, -6.8635425e-14, -1.5758214e-14, 0.0, 0.0, 0.0,
                    -0.003354688, -0.00013212219, 4.5390329e-05, -2.6017423e-07, 7.6728373e-08, 3.466561e-08, -7.9035553e-10, 1.7610549e-11, -3.1363309e-13, -4.741638e-13, 5.5748058e-14, -1.390074e-14, -1.8879465e-15, -1.0161655e-16, 0.0, 0.0,
                    -0.0003644111, -9.7545225e-05, 1.5509838e-05, 4.0179933e-06, 1.6956168e-07, 1.7457361e-09, 7.2904379e-10, 2.9305495e-11, 1.3563638e-13, -3.3491741e-13, 4.1828106e-14, 4.4231192e-15, -2.2486444e-17, -2.6053296e-17, 1.6362147e-18, 0.0,
                    -0.00052703131, -0.00016753741, -2.3173518e-06, 1.9601325e-06, 1.3155985e-07, 8.7049262e-11, -1.5396847e-10, 3.0639418e-11, 3.3787168e-12, 1.8243357e-15, 5.1578908e-15, 6.2635779e-16, 7.1664911e-17, 4.9580752e-17, -6.7885954e-19, 1.5503219e-19;

            sineCoefficients_ = Eigen::MatrixXd(maximumDegree + 1, maximumOrder + 1);
            sineCoefficients_ << 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0056001365, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0081157719, 0.0038666951, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0044228069, 0.00084562088, 0.00056399842, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.00040835265, 0.00013355072, 0.00036099797, 2.3202309e-05, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0013456091, 0.00031660688, 8.7896463e-06, 1.239398e-05, 7.7072017e-06, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.0014544906, 0.00044995953, 2.6311244e-05, 1.8718417e-06, 2.2078777e-06, 1.8563273e-07, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.00018339524, 0.00018110928, 5.0219647e-06, 8.2965753e-06, 3.7900317e-07, 3.4385676e-08, 4.9054745e-08, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.00020406948, 1.3913689e-05, 8.4511225e-06, 3.1466922e-07, 4.9911698e-07, 1.8092486e-08, 5.8710168e-09, 1.5713514e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.000467012, 6.984714e-05, 4.052549e-06, 7.9218511e-07, 1.9313979e-07, 3.406965e-08, 1.4011289e-09, 1.3832559e-10, 1.3236655e-10, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.00054143689, 1.1630132e-05, 6.4835704e-06, 2.2378128e-07, 4.3471159e-08, 9.0300331e-09, 5.7498839e-10, 1.0203755e-10, 1.8234256e-13, 7.5191487e-12, 0.0, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.00084160025, 5.506615e-05, 3.4092571e-07, 1.8967153e-07, 1.0107761e-08, 2.0706063e-09, 5.8926272e-10, 5.1348358e-12, 1.8695112e-12, 1.885248e-13, 1.1409691e-13, 0.0, 0.0, 0.0, 0.0,
                    0.0, 0.00056343034, 4.4633814e-05, 3.3261046e-07, 2.0043823e-07, 1.4608797e-08, 1.935467e-11, 3.4474363e-10, 2.3030528e-11, 1.1705956e-12, 1.2995393e-13, 6.812978e-14, 1.514277e-14, 0.0, 0.0, 0.0,
                    0.0, 9.0551537e-05, 2.2609028e-07, 1.4402096e-07, 6.4197295e-09, 8.3422668e-09, 1.1672606e-09, 1.9050421e-10, 1.9903467e-11, 1.2764944e-12, 1.4391851e-14, 3.2500778e-15, 7.7976884e-16, 3.1218531e-17, 0.0, 0.0,
                    0.0, 0.00010450677, 6.4306586e-06, 2.1110565e-06, 1.553349e-07, 1.4467997e-08, 6.1351496e-10, 2.9127184e-12, 6.4443237e-12, 6.0953923e-13, 4.1457548e-14, 1.2466192e-15, 8.6235475e-18, 5.0954031e-17, 1.3625581e-17, 0.0,
                    0.0, 8.138564e-05, 2.2953268e-05, 1.1112276e-06, 7.8452563e-10, 3.3764876e-09, 9.4043359e-10, 6.7121496e-12, 1.1717653e-12, 4.985168e-14, 2.5320711e-15, 9.5636359e-16, 6.2521874e-17, 1.882182e-17, 4.7202048e-19, 2.5665136e-19;

            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize spherical harmonics surface distribution model");
    }
}

SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings::SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings(
        SecondDegreeZonalPeriodicSurfacePropertyDistributionModel model) :
        SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings(TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN)
{
    model_ = model;

    switch(model)
    {
        case SecondDegreeZonalPeriodicSurfacePropertyDistributionModel::albedo_knocke:
        {
            // Knocke Earth albedo model: Knocke, Philip et al. "Earth radiation pressure effects on satellites." Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988.
            a0 = 0.34;
            c0 = 0;
            c1 = 0.1;
            c2 = 0;
            a2 = 0.29;
            referenceEpoch = spice_interface::convertDateStringToEphemerisTime("1981 DEC 22");
            period = physical_constants::JULIAN_YEAR_IN_DAYS;

            break;
        }
        case SecondDegreeZonalPeriodicSurfacePropertyDistributionModel::emissivity_knocke:
        {
            // Knocke Earth emissivity model: Knocke, Philip et al. "Earth radiation pressure effects on satellites." Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988.
            a0 = 0.68;
            c0 = 0;
            c1 = -0.07;
            c2 = 0;
            a2 = -0.18;
            referenceEpoch = spice_interface::convertDateStringToEphemerisTime("1981 DEC 22");
            period = physical_constants::JULIAN_YEAR_IN_DAYS;

            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize second-degree zonal periodic surface distribution model");
    }
}

std::shared_ptr<electromagnetism::LuminosityModel> createLuminosityModel(
        const std::shared_ptr<LuminosityModelSettings>& modelSettings,
        const std::string &body)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<LuminosityModel> luminosityModel;

    switch(modelSettings->getLuminosityModelType())
    {
        case LuminosityModelType::constant_radiant_power:
        {
            auto constantLuminosityModelSettings =
                    std::dynamic_pointer_cast< ConstantLuminosityModelSettings >(modelSettings);
            if(constantLuminosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected constant luminosity model for body " + body );
            }

            luminosityModel = std::make_shared<ConstantLuminosityModel>(
                constantLuminosityModelSettings->getLuminosity()
            );
            break;
        }
        case LuminosityModelType::irradiance_based_radiant_power:
        {
            auto irradianceBasedLuminosityModelSettings =
                    std::dynamic_pointer_cast< IrradianceBasedLuminosityModelSettings >(modelSettings);
            if(irradianceBasedLuminosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected irradiance-based luminosity model for body " + body );
            }

            luminosityModel = std::make_shared<IrradianceBasedLuminosityModel>(
                irradianceBasedLuminosityModelSettings->getIrradianceAtDistanceFunction(),
                irradianceBasedLuminosityModelSettings->getDistance()
            );
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize luminosity model settings for " + body );
    }

    return luminosityModel;
}

std::shared_ptr<electromagnetism::SurfacePropertyDistribution> createSurfacePropertyDistribution(
        const std::shared_ptr<SurfacePropertyDistributionSettings>& distributionSettings,
        const std::string& body)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<SurfacePropertyDistribution> surfacePropertyDistribution;

    switch(distributionSettings->getSurfacePropertyDistributionType())
    {
        case SurfacePropertyDistributionType::constant:
        {
            auto constantSurfacePropertyDistributionSettings =
                    std::dynamic_pointer_cast<ConstantSurfacePropertyDistributionSettings>(distributionSettings);
            if(constantSurfacePropertyDistributionSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected constant surface property distribution for body " + body );
            }

            surfacePropertyDistribution = std::make_shared<ConstantSurfacePropertyDistribution>(
                    constantSurfacePropertyDistributionSettings->getConstantValue());
            break;
        }
        case SurfacePropertyDistributionType::spherical_harmonics:
        {
            auto sphericalHarmonicsSurfacePropertyDistributionSettings =
                    std::dynamic_pointer_cast<SphericalHarmonicsSurfacePropertyDistributionSettings>(distributionSettings);
            if(sphericalHarmonicsSurfacePropertyDistributionSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected spherical harmonics surface property distribution for body " + body );
            }

            surfacePropertyDistribution = std::make_shared<SphericalHarmonicsSurfacePropertyDistribution>(
                    sphericalHarmonicsSurfacePropertyDistributionSettings->getCosineCoefficients(),
                    sphericalHarmonicsSurfacePropertyDistributionSettings->getSineCoefficients());
            break;
        }
        case SurfacePropertyDistributionType::second_degree_zonal_periodic:
        {
            auto secondDegreeZonalPeriodicSurfacePropertyDistributionSettings =
                    std::dynamic_pointer_cast<SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings>(distributionSettings);
            if(secondDegreeZonalPeriodicSurfacePropertyDistributionSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected second-degree zonal periodic surface property distribution for body " + body );
            }

            surfacePropertyDistribution = std::make_shared<SecondDegreeZonalPeriodicSurfacePropertyDistribution>(
                    secondDegreeZonalPeriodicSurfacePropertyDistributionSettings->getA0(),
                    secondDegreeZonalPeriodicSurfacePropertyDistributionSettings->getC0(),
                    secondDegreeZonalPeriodicSurfacePropertyDistributionSettings->getC1(),
                    secondDegreeZonalPeriodicSurfacePropertyDistributionSettings->getC2(),
                    secondDegreeZonalPeriodicSurfacePropertyDistributionSettings->getA2(),
                    secondDegreeZonalPeriodicSurfacePropertyDistributionSettings->getReferenceEpoch(),
                    secondDegreeZonalPeriodicSurfacePropertyDistributionSettings->getPeriod());
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize surface property distribution settings for " + body );
    }

    return surfacePropertyDistribution;
}

std::unique_ptr<electromagnetism::SourcePanelRadiosityModel> createPanelRadiosityModel(
        const std::shared_ptr<PanelRadiosityModelSettings>& modelSettings,
        const std::string& body)
{
    using namespace tudat::electromagnetism;

    std::unique_ptr<SourcePanelRadiosityModel> panelRadiosityModel;

    switch(modelSettings->getPanelRadiosityModelType())
    {
        case PanelRadiosityModelType::constant:
        {
            auto constantPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< ConstantPanelRadiosityModelSettings >(modelSettings);
            if(constantPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected constant panel radiosity model for body " + body );
            }

            panelRadiosityModel = std::make_unique<ConstantSourcePanelRadiosityModel>(
                    constantPanelRadiosityModelSettings->getConstantRadiosity());
            break;
        }
        case PanelRadiosityModelType::albedo:
        {
            auto albedoPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< AlbedoPanelRadiosityModelSettings >(modelSettings);
            if(albedoPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected albedo panel radiosity model for body " + body );
            }

            panelRadiosityModel = std::make_unique<AlbedoSourcePanelRadiosityModel>(
                    createSurfacePropertyDistribution(
                            albedoPanelRadiosityModelSettings->getAlbedoDistribution(), body));
            break;
        }
        case PanelRadiosityModelType::thermal_delayed:
        {
            auto delayedThermalPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< DelayedThermalPanelRadiosityModelSettings >(modelSettings);
            if(delayedThermalPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected delayed thermal panel radiosity model for body " + body );
            }

            panelRadiosityModel = std::make_unique<DelayedThermalSourcePanelRadiosityModel>(
                    createSurfacePropertyDistribution(
                            delayedThermalPanelRadiosityModelSettings->getEmissivityDistribution(), body));
            break;
        }
        case PanelRadiosityModelType::thermal_angle_based:
        {
            auto angleBasedThermalPanelRadiosityModelSettings =
                    std::dynamic_pointer_cast< AngleBasedThermalPanelRadiosityModelSettings >(modelSettings);
            if(angleBasedThermalPanelRadiosityModelSettings == nullptr)
            {
                throw std::runtime_error(
                        "Error, expected angle-based thermal panel radiosity model for body " + body );
            }

            panelRadiosityModel = std::make_unique<AngleBasedThermalSourcePanelRadiosityModel>(
                    angleBasedThermalPanelRadiosityModelSettings->getMinTemperature(),
                    angleBasedThermalPanelRadiosityModelSettings->getMaxTemperature(),
                    createSurfacePropertyDistribution(
                            angleBasedThermalPanelRadiosityModelSettings->getEmissivityDistribution(), body));
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize panel radiosity model settings for " + body );
    }

    return panelRadiosityModel;
}

std::shared_ptr<electromagnetism::RadiationSourceModel> createRadiationSourceModel(
        const std::shared_ptr<RadiationSourceModelSettings>& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<RadiationSourceModel> radiationSourceModel;

    switch(modelSettings->getRadiationSourceModelType())
    {
    case RadiationSourceModelType::isotropic_point_source:
    {
        auto isotropicPointModelSettings =
                std::dynamic_pointer_cast< IsotropicPointRadiationSourceModelSettings >(modelSettings);

        if(isotropicPointModelSettings == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected isotropic point radiation source for body " + body );
        }
        if(isotropicPointModelSettings->getLuminosityModelSettings() == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected isotropic point radiation source to have a luminosity model for body " + body);
        }

        auto luminosityModel = createLuminosityModel(
                isotropicPointModelSettings->getLuminosityModelSettings(), body);

        radiationSourceModel = std::make_shared<IsotropicPointRadiationSourceModel>(luminosityModel);
        break;
    }
    case RadiationSourceModelType::statically_paneled_source:
    {
        auto paneledModelSettings =
                std::dynamic_pointer_cast< StaticallyPaneledRadiationSourceModelSettings >(modelSettings);

        if(paneledModelSettings == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected statically paneled radiation source for body " + body );
        }
        if( paneledModelSettings->getOriginalSourceName().empty())
        {
            throw std::runtime_error(
                    "Error, expected statically paneled radiation source to have an original source for body " + body);
        }
        if(paneledModelSettings->getNumberOfPanels() == 0)
        {
            throw std::runtime_error(
                    "Error, expected statically paneled radiation source to have at least one panel for body " + body);
        }
        if(paneledModelSettings->getPanelRadiosityModelSettings().empty())
        {
            throw std::runtime_error(
                    "Error, expected statically paneled radiation source to have at least one panel radiosity model for body " + body);
        }

        auto sourceBody = bodies.getBody(body);

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
        for (auto& radiosityModelSetting : paneledModelSettings->getPanelRadiosityModelSettings())
        {
            radiosityModels.push_back(createPanelRadiosityModel(radiosityModelSetting, body));
        }

        radiationSourceModel = std::make_shared<StaticallyPaneledRadiationSourceModel>(
                paneledModelSettings->getOriginalSourceName(),
                sourceBody->getShapeModel(),
                radiosityModels,
                paneledModelSettings->getNumberOfPanels(),
                paneledModelSettings->getOriginalSourceToSourceOccultingBodies());
        break;
    }
    case RadiationSourceModelType::dynamically_paneled_source:
    {
        auto paneledModelSettings =
                std::dynamic_pointer_cast< DynamicallyPaneledRadiationSourceModelSettings >(modelSettings);

        if(paneledModelSettings == nullptr)
        {
            throw std::runtime_error(
                    "Error, expected dynamically paneled radiation source for body " + body );
        }
        if( paneledModelSettings->getOriginalSourceName().empty())
        {
            throw std::runtime_error(
                    "Error, expected dynamically paneled radiation source to have an original source for body " + body);
        }
        if(paneledModelSettings->getPanelRadiosityModelSettings().empty())
        {
            throw std::runtime_error(
                    "Error, expected dynamically paneled radiation source to have at least one panel radiosity model for body " + body);
        }

        auto sourceBody = bodies.getBody(body);

        std::vector<std::unique_ptr<SourcePanelRadiosityModel>> radiosityModels;
        for (auto& radiosityModelSetting : paneledModelSettings->getPanelRadiosityModelSettings())
        {
            radiosityModels.push_back(createPanelRadiosityModel(radiosityModelSetting, body));
        }

        radiationSourceModel = std::make_shared<DynamicallyPaneledRadiationSourceModel>(
                paneledModelSettings->getOriginalSourceName(),
                sourceBody->getShapeModel(),
                radiosityModels,
                paneledModelSettings->getNumberOfPanelsPerRing(),
                paneledModelSettings->getOriginalSourceToSourceOccultingBodies());
        break;
    }
    default:
        throw std::runtime_error( "Error, do not recognize radiation source model settings for " + body );
    }

    return radiationSourceModel;
}

} // tudat
} // electromagnetism
