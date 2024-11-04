#ifndef TUDATBUNDLE_LIBLRO_H
#define TUDATBUNDLE_LIBLRO_H

#include <memory>
#include <string>

#include "tudat/simulation/propagation_setup/propagationResults.h"


void loadLROSpiceKernels();
void saveSimulationResults(
        const std::shared_ptr<tudat::propagators::SingleArcSimulationResults<>>& propagationResults,
        const std::string resultsFolder,
        const bool saveHistory = true);

#endif //TUDATBUNDLE_LIBLRO_H
