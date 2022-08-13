#ifndef TUDAT_RADIATIONSOURCEMODEL_H
#define TUDAT_RADIATIONSOURCEMODEL_H

#include <vector>
#include <utility>
#include <functional>
#include <memory>

#include <Eigen/Core>

#include "tudat/astro/electromagnetism/luminosityModel.h"

namespace tudat
{
namespace electromagnetism
{

typedef std::vector<std::pair<double, Eigen::Vector3d>> IrradianceWithSourceList;


class RadiationSourceModel
{
public:
    explicit RadiationSourceModel() = default;

    virtual ~RadiationSourceModel() = default;

    void updateMembers(double currentTime );

    /*!
     * Evaluate the irradiance [W/m²] at a certain positition due to this source.
     * @param targetPosition Position where to evaluate the irradiance in local (i.e. source-fixed) coordinates
     * @return A list of irradiances and their source-fixed origin. Single element for point sources, multiple
     * elements for paneled sources.
     */
    virtual IrradianceWithSourceList evaluateIrradianceAtPosition(Eigen::Vector3d targetPosition) const = 0;

protected:
    virtual void updateMembers_(const double currentTime) {};

    double currentTime_{TUDAT_NAN};
};


class IsotropicPointRadiationSourceModel : public RadiationSourceModel
{
public:
    explicit IsotropicPointRadiationSourceModel(
            const std::shared_ptr<LuminosityModel> luminosityModel) :
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
