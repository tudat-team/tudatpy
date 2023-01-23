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
 *      Oliver Montenbruck, et al. Satellite Orbits. Springer Berlin Heidelberg, 2000.
 */

#ifndef TUDAT_LUMINOSITYMODEL_H
#define TUDAT_LUMINOSITYMODEL_H

#include <functional>
#include <memory>

#include "tudat/math/basic/mathematicalConstants.h"


namespace tudat
{
namespace electromagnetism
{

/*!
 * Class modeling the luminosity, or emitted electromagnetic power, of a radiation source.
 *
 * Only used by an isotropic point radiation source.
 */
class LuminosityModel
{
public:
    LuminosityModel() = default;

    virtual ~LuminosityModel() = default;

    /*!
     * Get the luminosity.
     *
     * @return Luminosity [W]
     */
    virtual double getLuminosity() const = 0;

    /*!
     * Update class members.
     *
     * @param currentTime Current simulation time
     */
    void updateMembers(double currentTime);

protected:
    virtual void updateMembers_(const double currentTime) {};

    double currentTime_{TUDAT_NAN};
};


/*!
 * Class modeling constant luminosity of a radiation source.
 */
class ConstantLuminosityModel : public LuminosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param luminosity Constant luminosity of the source [W]
     */
    explicit ConstantLuminosityModel(double luminosity) : luminosity_(luminosity) {}

    double getLuminosity() const override
    {
        return luminosity_;
    }

private:
    //! Constant luminosity of the source [W]
    double luminosity_;
};


/*!
 * Class modeling luminosity of a radiation source based on the irradiance at a given distance.
 *
 * This luminosity is useful to model variations in the Sun luminosity based on the solar irradiance measured
 * at Earth distance.
 *
 * The source is assumed to be isotropic such that the irradiance given would be found at any position with the
 * same distance from the source.
 */
class IrradianceBasedLuminosityModel : public LuminosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param irradianceAtDistanceFunction Function returning the irradiance [W/m²] at a given time
     * @param distance Distance from the source at which the irradiance was evaluated/measured
     */
    IrradianceBasedLuminosityModel(
            const std::function<double(double)>& irradianceAtDistanceFunction,
            double distance)
            : irradianceAtDistanceFunction_(irradianceAtDistanceFunction), distance_(distance) {}

    double getLuminosity() const override;

private:
    void updateMembers_(double currentTime) override;

    //! Function returning the irradiance [W/m²] at a given time
    std::function<double(double)> irradianceAtDistanceFunction_;

    //! Irradiance from the last evaluation of irradianceAtDistanceFunction_ [W/m²]
    double irradianceAtDistance_{};

    //! Distance from the source at which the irradiance was evaluated/measured
    double distance_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_LUMINOSITYMODEL_H
