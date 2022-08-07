#ifndef TUDAT_LUMINOSITYMODEL_H
#define TUDAT_LUMINOSITYMODEL_H

#include <functional>
#include <memory>

#include "tudat/math/basic/mathematicalConstants.h"

namespace tudat
{
namespace electromagnetism
{


class LuminosityModel
{
public:
    LuminosityModel() = default;

    virtual ~LuminosityModel() = default;

    virtual double getLuminosity() const = 0;

    virtual void updateMembers() {};
};


class ConstantLuminosityModel : public LuminosityModel
{
public:
    explicit ConstantLuminosityModel(double luminosity) : luminosity_(luminosity) {}

    double getLuminosity() const override
    {
        return luminosity_;
    }

private:
    double luminosity_;
};


class IrradianceBasedLuminosityModel : public LuminosityModel
{
public:
    IrradianceBasedLuminosityModel(
            const std::function<double()> irradianceAtDistanceFunction,
            double distance)
            : irradianceAtDistanceFunction_(irradianceAtDistanceFunction), distance_(distance) {}

    IrradianceBasedLuminosityModel(
            double irradianceAtDistance,
            double distance)
            : IrradianceBasedLuminosityModel([=] () { return irradianceAtDistance; }, distance) {}

    double getLuminosity() const override;

    void updateMembers() override;

private:
    std::function<double()> irradianceAtDistanceFunction_;
    double irradianceAtDistance_{};
    double distance_;
};

}
}

#endif //TUDAT_LUMINOSITYMODEL_H
