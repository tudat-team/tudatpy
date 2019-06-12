//
// Created by ggarrett on 10-6-19.
//

#ifndef TUDATBUNDLE_TYPEDEFS_H
#define TUDATBUNDLE_TYPEDEFS_H

#include "Tudat/SimulationSetup/tudatSimulationHeader.h"

namespace tss = tudat::simulation_setup;

typedef std::shared_ptr<tss::BodySettings> stdBodySettingsPointer;
typedef boost::shared_ptr<tss::BodySettings> boostBodySettingsPointer;
typedef std::map<std::string, stdBodySettingsPointer> stdBodySettingsPointerMap;
typedef std::map<std::string, boostBodySettingsPointer> boostBodySettingsPointerMap;

#endif //TUDATBUNDLE_TYPEDEFS_H
