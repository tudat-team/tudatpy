/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_CREATEOCCULTATIONMODEL_H
#define TUDAT_CREATEOCCULTATIONMODEL_H

#include <memory>

#include "tudat/simulation/environment_setup/body.h"

namespace tudat
{
namespace simulation_setup
{

/*!
 * Settings for an occultation model. This model should be assigned to the observer (for 
 * radiation pressure, the target), as opposed to the source/occulted body.
 */
class OccultationModelSettings
{
public:
    /*!
     * Constructor.
     
     * @param occultingBodies Names of bodies to occult a source as seen from an observer 
     */
    explicit OccultationModelSettings(const std::vector<std::string>& occultingBodies) :
            occultingBodies_(occultingBodies) {}

    std::vector<std::string> getOccultingBodies() const
    {
        return occultingBodies_;
    }

private:
    std::vector<std::string> occultingBodies_
};

/*!
 * Create settings for an occultation model.
 *
 * @param occultingBodies Names of bodies to occult a source as seen from an observer/target
 * @return Shared pointer to settings for an occultation model
 */
inline std::shared_ptr<OccultationModelSettings>
occultationModelSettings(const std::vector<std::string>& occultingBodies = {})
{
    return std::make_shared<OccultationModelSettings>(occultingBodies);
}

/*!
 * Create occultation model from its settings.
 *
 * @param modelSettings Settings of the radiation source model
 * @param body Body to which the radiation source model belongs
 * @param bodies System of bodies
 * @return Shared pointer to radiation source model
 */
std::shared_ptr<electromagnetism::OccultationModel> createOccultationModel(
        const std::shared_ptr< OccultationModelSettings >& modelSettings,
        const std::string& body,
        const SystemOfBodies& bodies);


} // tudat
} // simulation_setup

#endif //TUDAT_CREATEOCCULTATIONMODEL_H
