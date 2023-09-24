/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include "tudat/simulation/environment_setup/createOccultationModel.h"

#include "tudat/astro/electromagnetism/occultationModel.h"
#include "tudat/simulation/environment_setup/body.h"


namespace tudat
{
namespace simulation_setup
{

std::shared_ptr<tudat::electromagnetism::OccultationModel> createOccultationModel(
        const std::vector<std::string>& occultingBodies,
        const SystemOfBodies& bodies)
{
    using namespace tudat::electromagnetism;

    std::shared_ptr<OccultationModel> occultationModel;
    switch(occultingBodies.size())
    {
        case 0:
        {
            occultationModel = std::make_shared<NoOccultingBodyOccultationModel>();
            break;
        }
        case 1:
        {
            auto occultingBodyName = occultingBodies.front();
            auto occultingBody = bodies.at(occultingBodyName);
            occultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
                    occultingBodyName,
                    std::bind( &Body::getPosition, occultingBody),
                    occultingBody->getShapeModel());
            break;
        }
        default:
        {
            std::vector<std::function<Eigen::Vector3d()>> occultingBodyPositionFunctions{};
            std::vector<std::shared_ptr<basic_astrodynamics::BodyShapeModel>> occultingBodyShapeModels;
            for (const auto& occultingBodyName : occultingBodies)
            {
                auto occultingBody = bodies.at(occultingBodyName);
                occultingBodyPositionFunctions.emplace_back(std::bind( &Body::getPosition, occultingBody));
                occultingBodyShapeModels.push_back(occultingBody->getShapeModel());
            }
            occultationModel = std::make_shared<SimpleMultipleOccultingBodyOccultationModel>(
                    occultingBodies, occultingBodyPositionFunctions, occultingBodyShapeModels);
            break;
        }
    }

    return occultationModel;
}

} // tudat
} // electromagnetism


#include "tudat/simulation/environment_setup/createOccultationModel.h"
