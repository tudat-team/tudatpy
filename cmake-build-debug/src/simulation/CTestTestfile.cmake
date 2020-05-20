# CMake generated Testfile for 
# Source directory: /home/ggarrett/Repositories/new/tudat_bundle/tudat/src/simulation
# Build directory: /home/ggarrett/Repositories/new/tudat_bundle/tudat/cmake-build-debug/src/simulation
# 
# This file includes the relevant testing commands required for 
# testing this directory and lists subdirectories to be tested as well.
add_test(test_simulation_EnvironmentModelSetup "/home/ggarrett/Repositories/new/tudat_bundle/tudat/cmake-build-debug/tests/test_simulation_EnvironmentModelSetup")
set_tests_properties(test_simulation_EnvironmentModelSetup PROPERTIES  _BACKTRACE_TRIPLES "/home/ggarrett/Repositories/new/tudat_bundle/tudat/CMakeLists.txt;323;add_test;/home/ggarrett/Repositories/new/tudat_bundle/tudat/src/simulation/CMakeLists.txt;21;TUDAT_ADD_TESTCASE;/home/ggarrett/Repositories/new/tudat_bundle/tudat/src/simulation/CMakeLists.txt;0;")
add_test(test_simulation_AccelerationModelSetup "/home/ggarrett/Repositories/new/tudat_bundle/tudat/cmake-build-debug/tests/test_simulation_AccelerationModelSetup")
set_tests_properties(test_simulation_AccelerationModelSetup PROPERTIES  _BACKTRACE_TRIPLES "/home/ggarrett/Repositories/new/tudat_bundle/tudat/CMakeLists.txt;323;add_test;/home/ggarrett/Repositories/new/tudat_bundle/tudat/src/simulation/CMakeLists.txt;25;TUDAT_ADD_TESTCASE;/home/ggarrett/Repositories/new/tudat_bundle/tudat/src/simulation/CMakeLists.txt;0;")
subdirs("environment")
subdirs("propagation")
