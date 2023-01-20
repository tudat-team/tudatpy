#include <iostream>
#include <cstdlib>

#include <Eigen/Core>

#include "tudat/astro/electromagnetism/occultationModel.h"
#include "tudat/astro/basic_astro/sphericalBodyShapeModel.h"

using namespace tudat::electromagnetism;

Eigen::Vector3d generatePointOnSource();
double estimateShadowFunction(OccultationModel& occultationModel);

const int nSamples = 1e7;
const double sourceRadius = 2.5;
const double occultingBodyRadius = 1;
const Eigen::Vector3d sourcePosition = Eigen::Vector3d::Zero();
const Eigen::Vector3d occultingBodyPosition = Eigen::Vector3d(0, 0, 4);
const Eigen::Vector3d targetPosition = Eigen::Vector3d(0, 0.2, 9);

/*
 * Estimate the shadow function by evaluating point-to-point visibility of evenly distributed
 * points on source. Ideally, the result should agree with the geometric shadow function
 * evaluation.
 */
int main()
{
    srand(time(nullptr));

    auto sourceShape = std::make_shared<tudat::basic_astrodynamics::SphericalBodyShapeModel>(sourceRadius);
    auto occultingBodyShape = std::make_shared<tudat::basic_astrodynamics::SphericalBodyShapeModel>(occultingBodyRadius);

    SingleOccultingBodyOccultationModel occultationModel(
            {}, [=]() { return occultingBodyPosition; }, occultingBodyShape);
    occultationModel.updateMembers(TUDAT_NAN);

    std::cout << "Actual shadow function: "
            << occultationModel.evaluateReceivedFractionFromExtendedSource(
                    sourcePosition,
                    sourceShape,
                    targetPosition) << std::endl;

    std::cout << "Estimated shadow function: "
            << estimateShadowFunction(occultationModel) << std::endl;

    return 0;
}

double estimateShadowFunction(OccultationModel& occultationModel)
{
    int nVisible = 0;
    for (int i = 0; i < nSamples; ++i)
    {
        auto pointOnSource = generatePointOnSource();
        if (occultationModel.evaluateReceivedFractionFromPointSource(pointOnSource, targetPosition))
        {
            nVisible++;
        }
    }
    return nVisible / static_cast<double>(nSamples);
}

Eigen::Vector3d generatePointOnSource()
{
    while (true)
    {
        double x = rand() / (RAND_MAX + 1.);
        double y = rand() / (RAND_MAX + 1.);
        if (x*x + y*y < 1)
        {
            return Eigen::Vector3d(x, y, 0) * sourceRadius;
        }
    }
}
