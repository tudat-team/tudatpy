/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createSurfacePropertyDistribution.h"

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
            // Original coefficients: https://github.com/DominikStiller/tudelft-hpb-project/blob/main/analysis/files/DLAM-1.txt
            // Initializer lists generated in https://github.com/DominikStiller/tudelft-hpb-project/blob/0a73ae8d03022852ff89874fbebf1eb7c58ea8e9/analysis/lunar_models.ipynb

            cosineCoefficients_ = Eigen::MatrixXd::Zero( 16, 16 );
            cosineCoefficients_ <<
                +1.9467246e-01, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +3.7619635e-03, -3.6562523e-02, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +2.6013759e-02, -2.3342472e-02, -1.5511737e-03, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +1.0791359e-02, +1.2984099e-02, -5.1947501e-04, +8.3037838e-04, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +1.5490099e-03, +2.4150199e-03, -4.5980599e-04, -5.8965344e-05, +1.0396812e-04, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                -4.7391417e-03, +9.7388352e-04, -6.6144592e-04, +3.6299540e-05, +3.9759108e-06, +2.2706280e-06, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +2.1741056e-03, +1.9854158e-03, +3.6293963e-05, -2.7638060e-05, -4.6429902e-06, -3.8219121e-06, +3.8641269e-07, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                -2.6214758e-03, +9.6586546e-04, -7.8726437e-05, +3.5085491e-05, +3.0497316e-06, -4.0231517e-08, -5.2535643e-08, -1.3545976e-08, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                -7.2392123e-03, -3.6002321e-05, -1.5504282e-04, -4.2175334e-06, +2.0927038e-06, -2.1213522e-07, +5.7936412e-08, -2.8380071e-09, +1.0183988e-09, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                -4.6530150e-03, +7.7363813e-04, +5.8165535e-05, -2.2643899e-05, -3.9981338e-07, -2.2805609e-08, -2.1918927e-08, +4.0524183e-09, +5.7278542e-10, +1.0169170e-10, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                -1.5925583e-02, -1.4417769e-05, +5.8118004e-05, -1.0071993e-06, +3.8769952e-07, -1.2667752e-08, +8.7088087e-09, +1.0255584e-09, +3.5006160e-10, +3.1475729e-11, +2.3275218e-12, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                -3.7225793e-03, -4.1507339e-04, -6.3009829e-05, -2.5513158e-06, +1.5175304e-07, -5.7092276e-08, +1.9616945e-09, +3.1055117e-10, +3.6901893e-11, +5.4058712e-12, +6.7325209e-14, +8.2735779e-14, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                -7.7462959e-03, +6.8495466e-04, +1.4207584e-04, -1.0655919e-06, -3.4816920e-07, +5.7062063e-09, +5.0926848e-11, +9.6192077e-11, +6.7701005e-12, +2.1793964e-13, -1.3931937e-13, -6.8635425e-14, -1.5758214e-14, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                -3.3546880e-03, -1.3212219e-04, +4.5390329e-05, -2.6017423e-07, +7.6728373e-08, +3.4665610e-08, -7.9035553e-10, +1.7610549e-11, -3.1363309e-13, -4.7416380e-13, +5.5748058e-14, -1.3900740e-14, -1.8879465e-15, -1.0161655e-16, +0.0000000e+00, +0.0000000e+00,
                -3.6441110e-04, -9.7545225e-05, +1.5509838e-05, +4.0179933e-06, +1.6956168e-07, +1.7457361e-09, +7.2904379e-10, +2.9305495e-11, +1.3563638e-13, -3.3491741e-13, +4.1828106e-14, +4.4231192e-15, -2.2486444e-17, -2.6053296e-17, +1.6362147e-18, +0.0000000e+00,
                -5.2703131e-04, -1.6753741e-04, -2.3173518e-06, +1.9601325e-06, +1.3155985e-07, +8.7049262e-11, -1.5396847e-10, +3.0639418e-11, +3.3787168e-12, +1.8243357e-15, +5.1578908e-15, +6.2635779e-16, +7.1664911e-17, +4.9580752e-17, -6.7885954e-19, +1.5503219e-19;
            // DLAM-1 was derived for 750 nm, for which the lunar albedo is 1.3x higher than for the average solar wavelength (Vasavada 2012)
            cosineCoefficients_ /= 1.3;

            sineCoefficients_ = Eigen::MatrixXd::Zero( 16, 16 );
            sineCoefficients_ <<
                +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, +5.6001365e-03, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, +8.1157719e-03, +3.8666951e-03, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, +4.4228069e-03, +8.4562088e-04, +5.6399842e-04, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, +4.0835265e-04, +1.3355072e-04, -3.6099797e-04, +2.3202309e-05, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, -1.3456091e-03, +3.1660688e-04, -8.7896463e-06, -1.2393980e-05, -7.7072017e-06, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, -1.4544906e-03, -4.4995953e-04, -2.6311244e-05, -1.8718417e-06, -2.2078777e-06, +1.8563273e-07, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, +1.8339524e-04, -1.8110928e-04, +5.0219647e-06, +8.2965753e-06, -3.7900317e-07, -3.4385676e-08, +4.9054745e-08, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, +2.0406948e-04, -1.3913689e-05, -8.4511225e-06, +3.1466922e-07, +4.9911698e-07, -1.8092486e-08, +5.8710168e-09, -1.5713514e-10, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, +4.6701200e-04, +6.9847140e-05, -4.0525490e-06, -7.9218511e-07, -1.9313979e-07, +3.4069650e-08, -1.4011289e-09, -1.3832559e-10, +1.3236655e-10, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, -5.4143689e-04, -1.1630132e-05, -6.4835704e-06, +2.2378128e-07, +4.3471159e-08, -9.0300331e-09, +5.7498839e-10, -1.0203755e-10, +1.8234256e-13, +7.5191487e-12, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, -8.4160025e-04, -5.5066150e-05, +3.4092571e-07, +1.8967153e-07, +1.0107761e-08, -2.0706063e-09, +5.8926272e-10, -5.1348358e-12, -1.8695112e-12, +1.8852480e-13, +1.1409691e-13, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, -5.6343034e-04, -4.4633814e-05, -3.3261046e-07, +2.0043823e-07, -1.4608797e-08, +1.9354670e-11, +3.4474363e-10, +2.3030528e-11, +1.1705956e-12, +1.2995393e-13, +6.8129780e-14, +1.5142770e-14, +0.0000000e+00, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, +9.0551537e-05, +2.2609028e-07, +1.4402096e-07, -6.4197295e-09, -8.3422668e-09, -1.1672606e-09, -1.9050421e-10, +1.9903467e-11, -1.2764944e-12, +1.4391851e-14, +3.2500778e-15, +7.7976884e-16, +3.1218531e-17, +0.0000000e+00, +0.0000000e+00,
                +0.0000000e+00, -1.0450677e-04, +6.4306586e-06, +2.1110565e-06, -1.5533490e-07, +1.4467997e-08, -6.1351496e-10, -2.9127184e-12, +6.4443237e-12, +6.0953923e-13, +4.1457548e-14, +1.2466192e-15, -8.6235475e-18, -5.0954031e-17, -1.3625581e-17, +0.0000000e+00,
                +0.0000000e+00, +8.1385640e-05, -2.2953268e-05, +1.1112276e-06, +7.8452563e-10, -3.3764876e-09, -9.4043359e-10, -6.7121496e-12, -1.1717653e-12, +4.9851680e-14, +2.5320711e-15, -9.5636359e-16, +6.2521874e-17, -1.8821820e-17, +4.7202048e-19, +2.5665136e-19;
            sineCoefficients_ /= 1.3;

            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize spherical harmonics surface distribution model");
    }
}

SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings::SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings(
        KnockeTypeSurfacePropertyDistributionModel model) :
        SecondDegreeZonalPeriodicSurfacePropertyDistributionSettings(TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN, TUDAT_NAN)
{
    model_ = model;

    switch(model)
    {
        case KnockeTypeSurfacePropertyDistributionModel::albedo_knocke:
        {
            // Knocke Earth albedo model: Knocke, Philip et al. "Earth radiation pressure effects on satellites." Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988.
            a0 = 0.34;
            c0 = 0;
            c1 = 0.10;
            c2 = 0;
            a2 = 0.29;
            referenceEpoch = spice_interface::convertDateStringToEphemerisTime("1981-12-22 00:00:00 UTC");
            period = physical_constants::JULIAN_YEAR_IN_DAYS;

            break;
        }
        case KnockeTypeSurfacePropertyDistributionModel::emissivity_knocke:
        {
            // Knocke Earth emissivity model: Knocke, Philip et al. "Earth radiation pressure effects on satellites." Astrodynamics Conference. American Institute of Aeronautics and Astronautics, 1988.
            a0 = 0.68;
            c0 = 0;
            c1 = -0.07;
            c2 = 0;
            a2 = -0.18;
            referenceEpoch = spice_interface::convertDateStringToEphemerisTime("1981-12-22 00:00:00 UTC");
            period = physical_constants::JULIAN_YEAR_IN_DAYS;

            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize second-degree zonal periodic surface distribution model");
    }
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
        case SurfacePropertyDistributionType::custom_surface_distribution:
        {
            auto customSurfacePropertyDistributionSettings =
                std::dynamic_pointer_cast<CustomSurfacePropertyDistributionSettings>(distributionSettings);
            if(customSurfacePropertyDistributionSettings == nullptr)
            {
                throw std::runtime_error(
                    "Error, expected custom surface property distribution for body " + body );
            }

            surfacePropertyDistribution = std::make_shared< CustomSurfacePropertyDistribution >(
                customSurfacePropertyDistributionSettings->getCustomFunction( ) );
            break;
        }
        default:
            throw std::runtime_error( "Error, do not recognize surface property distribution settings for " + body );
    }

    return surfacePropertyDistribution;
}


} // tudat
} // electromagnetism
