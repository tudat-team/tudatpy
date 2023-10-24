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

double computeLuminosityFromIrradiance( const double irradiance, const double distance );

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

class VariableLuminosityModel : public LuminosityModel
{
public:
    /*!
     * Constructor.
     *
     * @param luminosity Constant luminosity of the source [W]
     */
    explicit VariableLuminosityModel( const std::function< double( const double ) > luminosityFunction ) :
        luminosityFunction_( luminosityFunction ), luminosity_( TUDAT_NAN ) {}

    double getLuminosity() const override
    {
        return luminosity_;
    }

private:

    virtual void updateMembers_(const double currentTime)  override
    {
        luminosity_ = luminosityFunction_( currentTime );
    }

    //! Constant luminosity of the source [W]
    const std::function< double( const double ) > luminosityFunction_;

    double luminosity_;
};


} // tudat
} // electromagnetism

#endif //TUDAT_LUMINOSITYMODEL_H
