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

namespace tudat
{
namespace simulation_setup
{

using namespace electromagnetism;

std::shared_ptr<electromagnetism::OccultationModel> createOccultationModel(
        const std::shared_ptr<LuminosityModelSettings>& modelSettings,
        const std::string &body,
        const SystemOfBodies& bodies)
{
    using namespace tudat::electromagnetism;

    auto occultingBodies = modelSettings->occultingBodies_;

    std::shared_ptr<OccultationModel> occultationModel;
    switch(occultingBodies.size())
    {
        case 0:
        {
            occultationModel = std::make_shared<NoOccultingBodyOccultationModel>(occultingBodies);
            break;
        }
        case 1:
        {
            auto occultingBodyName = occultingBodies_.front();
            auto occultingBody = bodies.at(occultingBodyName);
            occultationModel = std::make_shared<SingleOccultingBodyOccultationModel>(
                    occultingBodies,
                    std::bind( &Body::getPosition, occultingBody),
                    occultingBody->getShapeModel());
            break;
        }
        default:
        {
            throw std::runtime_error( "Error, only a single occulting body is supported for target body " + body);
        }
    }

    return occultationModel;
}


} // tudat
} // electromagnetism
