#ifndef TUDAT_RADIATIONSOURCEMODEL_H
#define TUDAT_RADIATIONSOURCEMODEL_H

#include <vector>
#include <tuple>
#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/electromagnetism/luminosityModel.h"

namespace tudat
{
namespace electromagnetism
{

typedef std::vector<std::tuple<double, Eigen::Vector3d>> IrradianceWithSourceList;


class RadiationSourceModel
{
public:
    explicit RadiationSourceModel(const std::function< Eigen::Vector3d( ) > sourcePositionFunction):
            sourcePositionFunction_(sourcePositionFunction) {}

    virtual ~RadiationSourceModel() = default;

    void updateMembers(double currentTime );

    virtual IrradianceWithSourceList evaluateIrradianceAtPosition(Eigen::Vector3d) const = 0;

protected:
    virtual void updateMembers_(const double currentTime) {};

    Eigen::Vector3d sourcePosition_;
    std::function< Eigen::Vector3d( ) > sourcePositionFunction_;
    double currentTime_{TUDAT_NAN};
};


class IsotropicPointRadiationSourceModel : public RadiationSourceModel
{
public:
    explicit IsotropicPointRadiationSourceModel(
            const std::function< Eigen::Vector3d( ) > sourcePositionFunction,
            const std::shared_ptr<LuminosityModel> luminosityModel):
            RadiationSourceModel(sourcePositionFunction),
            luminosityModel_(luminosityModel) {}

    IrradianceWithSourceList evaluateIrradianceAtPosition(Eigen::Vector3d targetPosition) const override;

    std::shared_ptr<LuminosityModel> getLuminosityModel() const
    {
        return luminosityModel_;
    }

private:
    void updateMembers_(double currentTime) override;

    std::shared_ptr<LuminosityModel> luminosityModel_;
};

} // tudat
} // electromagnetism

#endif //TUDAT_RADIATIONSOURCEMODEL_H
