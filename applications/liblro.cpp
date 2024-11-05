#include "liblro.h"

#include <string>

#include <boost/filesystem.hpp>
#include <boost/algorithm/string/predicate.hpp>
#include <boost/range/adaptors.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/io/basicInputOutput.h"


// Loads SPICE kernels for LRO
void loadLROSpiceKernels()
{
    using namespace tudat::spice_interface;

    std::string basePath;
    if(const char* env_p = std::getenv("SPICE_BASE"))
    {
        basePath = std::string(env_p);
    }
    else
    {
        basePath = "./";
    }

    std::string path = basePath + "lro/data";

    // Leap seconds
    loadSpiceKernelInTudat(path + "/lsk/naif0012.tls");
    // Planetary orientation shapes
    loadSpiceKernelInTudat(path + "/pck/pck00010.tpc");
    // Planetary gravitational parameters
    loadSpiceKernelInTudat(path + "/pck/gm_de431.tpc");
    // Lunar frame
    loadSpiceKernelInTudat(path + "/fk/moon_080317.tf");
    loadSpiceKernelInTudat(path + "/pck/moon_pa_de421_1900_2050.bpc");

    // LRO spacecraft bus and instrument frames
    loadSpiceKernelInTudat(path + "/fk/lro_frames_2012255_v02.tf");
    // LRO spacecraft clock
    loadSpiceKernelInTudat(path + "/sclk/lro_clkcor_2022075_v00.tsc");

    // LRO ephemeris
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path + "/spk"), {})) {
        if (entry.path().extension() == ".bsp")
        {
            loadSpiceKernelInTudat(entry.path().string());
        }
    }

    // LRO orientation
    for(auto& entry : boost::make_iterator_range(boost::filesystem::directory_iterator(path + "/ck"), {})) {
        if (entry.path().extension() == ".bc")
        {
            loadSpiceKernelInTudat(entry.path().string());
        }
    }
}

void saveSimulationResults(
        const std::shared_ptr<tudat::propagators::SingleArcSimulationResults<>>& propagationResults,
        const std::string resultsFolder,
        const bool saveHistory)
{
    using namespace tudat::input_output;

    if (saveHistory)
    {
        auto stateHistory = propagationResults->getEquationsOfMotionNumericalSolution();
        auto dependentVariableHistory = propagationResults->getDependentVariableHistory();
        auto dependentVariableNames = propagationResults->getDependentVariableId();
        auto stateNames = propagationResults->getProcessedStateIds();

        writeDataMapToTextFile(stateHistory,
                               "state_history.csv",
                               resultsFolder,
                               "",
                               std::numeric_limits< double >::digits10,
                               std::numeric_limits< double >::digits10,
                               ",");

        writeDataMapToTextFile(dependentVariableHistory,
                               "dependent_variable_history.csv",
                               resultsFolder,
                               "",
                               std::numeric_limits< double >::digits10,
                               std::numeric_limits< double >::digits10,
                               ",");

        writeIdMapToTextFile(dependentVariableNames,
                             "dependent_variable_names.csv",
                             resultsFolder,
                             ";");

        writeIdMapToTextFile(stateNames,
                             "state_names.csv",
                             resultsFolder,
                             ";");
    }

    auto cpuTimeHistory = propagationResults->getCumulativeComputationTimeHistory();
    writeDataMapToTextFile(cpuTimeHistory,
                           "cpu_time.csv",
                           resultsFolder,
                           "",
                           std::numeric_limits< double >::digits10,
                           std::numeric_limits< double >::digits10,
                           ",");
}
