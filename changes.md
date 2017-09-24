##### Changes Tudat/tudat:development ... aleixpinardell/tudat:json

* Ignore .DS_Store files
      .gitignore
    
* Typedef for `std::map< std::string, std::vector< boost::shared_ptr< AccelerationSettings > > >` as `SingleSelectedAccelerationMap`
      Tudat/Astrodynamics/Aerodynamics/UnitTests/unitTestAerodynamicCoefficientsFromFile.cpp
      Tudat/Astrodynamics/Aerodynamics/UnitTests/unitTestAerodynamicMomentAndAerodynamicForce.cpp
      Tudat/Astrodynamics/Aerodynamics/UnitTests/unitTestControlSurfaceIncrements.cpp
      Tudat/Astrodynamics/Aerodynamics/UnitTests/unitTestTabulatedAerodynamicCoefficients.cpp
      Tudat/Astrodynamics/Aerodynamics/UnitTests/unitTestWindModel.cpp
      Tudat/Astrodynamics/BasicAstrodynamics/UnitTests/unitTestEmpiricalAcceleration.cpp
      Tudat/Astrodynamics/OrbitDetermination/UnitTests/orbitDeterminationTestCases.h
      Tudat/Astrodynamics/OrbitDetermination/UnitTests/unitTestMultiArcStateEstimation.cpp
      Tudat/Astrodynamics/OrbitDetermination/UnitTests/unitTestParameterInfluenceDetermination.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/uniTestPropagationTerminationCheckOnFinalStep.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestCowellStateDerivative.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestCustomStatePropagation.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestDependentVariableOutput.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestEnckeStateDerivative.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestGaussStateDerivative.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestMultiArcDynamics.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestMultiArcVariationalEquationPropagation.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestMultiTypeStatePropagation.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestPropagationTerminationReason.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestRotationalDynamicsPropagator.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestSequentialVariationalEquationIntegration.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestStoppingConditions.cpp
      Tudat/Astrodynamics/Propagators/UnitTests/unitTestVariationalEquationPropagation.cpp
      Tudat/Astrodynamics/Propulsion/UnitTests/unitTestThrustAcceleration.cpp
      Tudat/Astrodynamics/Relativity/UnitTests/unitTestRelativisticAccelerationCorrection.cpp
      Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h
      Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.cpp
      
* Function `getAtmosphereTablesPath`
      Tudat/Astrodynamics/Aerodynamics/UnitTests/unitTestNRLMSISE00Atmosphere.cpp
      Tudat/Astrodynamics/Aerodynamics/UnitTests/unitTestTabulatedAtmosphere.cpp
      Tudat/InputOutput/basicInputOutput.h
      Tudat/SimulationSetup/EnvironmentSetup/createAtmosphereModel.cpp
      Tudat/SimulationSetup/UnitTests/unitTestEnvironmentModelSetup.cpp
      
* Rename `central_gravity` to `point_mass_gravity`
      Tudat/Astrodynamics/BasicAstrodynamics/accelerationModelTypes.h
      
* Typedef for `std::map< std::string, std::vector< boost::shared_ptr< MassRateModel > > >` as `MassRateModelMap`
      Tudat/Astrodynamics/BasicAstrodynamics/massRateModel.h

* Provide absolute path for non-standard Spice kernels:
      Tudat/Astrodynamics/Gravitation/UnitTests/unitTestGravityFieldVariations.cpp
      Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/UnitTests/unitTestSphericalHarmonicPartials.cpp
      
* SH gravity fields coefficients in External directory
      Tudat/Astrodynamics/Gravitation/egm96_coefficients.dat
      Tudat/Astrodynamics/Gravitation/gglp_lpe200_sha.tab
      Tudat/Astrodynamics/Gravitation/jgmro_120d_sha.tab
      Tudat/Astrodynamics/OrbitDetermination/AccelerationPartials/UnitTests/unitTestSphericalHarmonicPartials.cpp
      Tudat/External/GravityModels/Earth/egm96.txt
      Tudat/External/GravityModels/Earth/ggm02c.txt
      Tudat/External/GravityModels/Earth/ggm02s.txt
      Tudat/External/GravityModels/Mars/jgmro120d.txt
      Tudat/External/GravityModels/Moon/glgm3150.txt
      Tudat/External/GravityModels/Moon/lpe200.txt
      Tudat/InputOutput/basicInputOutput.h
      
* Remove gomp
      Tudat/Astrodynamics/GroundStations/CMakeLists.txt
      
* Fix adding cpp in list of headers
      Tudat/Astrodynamics/Propagators/CMakeLists.txt
      
* Fix typo `body_transational_state_update`
      Tudat/Astrodynamics/Propagators/environmentUpdateTypes.h
      Tudat/Astrodynamics/Propagators/singleStateTypeDerivative.h
      
* CPU time
      Tudat/Astrodynamics/Propagators/integrateEquations.h
      Tudat/Mathematics/NumericalIntegrators/numericalIntegrator.h
      Tudat/Mathematics/NumericalIntegrators/rungeKutta4Integrator.h
      Tudat/Mathematics/NumericalIntegrators/rungeKuttaVariableStepSizeIntegrator.h
      Tudat/SimulationSetup/PropagationSetup/dynamicsSimulator.h
      Tudat/SimulationSetup/PropagationSetup/propagationTermination.cpp
      Tudat/SimulationSetup/PropagationSetup/propagationTermination.h
      Tudat/SimulationSetup/PropagationSetup/propagationTerminationSettings.h
      Tudat/SimulationSetup/PropagationSetup/variationalEquationsSolver.h
      
* Add `DataInterpolationSettings` class
      Tudat/Astrodynamics/Propulsion/UnitTests/unitTestThrustAcceleration.cpp
      Tudat/Mathematics/Interpolators/createInterpolator.h
      Tudat/SimulationSetup/PropagationSetup/accelerationSettings.h
      
* Fix typo `associatedThroustSource_`
      Tudat/Astrodynamics/Propulsion/thrustAccelerationModel.h
      Tudat/SimulationSetup/PropagationSetup/createMassRateModels.cpp
      
* Space weather data to External directory
      Tudat/Astrodynamics/Aerodynamics/sw19571001.txt
      Tudat/External/SpaceWeatherData/sw19571001.txt
      Tudat/InputOutput/basicInputOutput.h
      Tudat/SimulationSetup/UnitTests/unitTestEnvironmentModelSetup.cpp
      
* Fix typo `abberation`
      Tudat/External/SpiceInterface/UnitTests/unitTestSpiceInterface.cpp
      Tudat/External/SpiceInterface/spiceEphemeris.cpp
      Tudat/External/SpiceInterface/spiceEphemeris.h
      Tudat/External/SpiceInterface/spiceInterface.cpp
      Tudat/External/SpiceInterface/spiceInterface.h
      Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.cpp
      Tudat/SimulationSetup/EnvironmentSetup/createEphemeris.h
      
* Make `readMatrixFromFile` templated
      Tudat/InputOutput/CMakeLists.txt
      Tudat/InputOutput/matrixTextFileReader.h
      Tudat/InputOutput/matrixTextFileReader.cpp
      
* Add overload of `writeDataMapToTextFile`
      Tudat/InputOutput/basicInputOutput.h
      
* Add `readEigenVectorMapFromFile` Function
      Tudat/InputOutput/mapTextFileReader.h

* More info in the error message when reading tabulated aerodynamics coefficients file
      Tudat/InputOutput/multiDimensionalArrayReader.h
      
* Fix compile error
      Tudat/Mathematics/Interpolators/jumpDataLinearInterpolator.h
      
* Remove `const`
      Tudat/Mathematics/NumericalIntegrators/createNumericalIntegrator.h
      
* Add class `FromFileAerodynamicCoefficientSettings` and make `TabulatedAerodynamicCoefficientSettings<N>` and `TabulatedAerodynamicCoefficientSettings<1>` derive from it
      Tudat/SimulationSetup/EnvironmentSetup/createAerodynamicCoefficientInterface.h

* Add `constantMass` property for `BodySettings`
      Tudat/SimulationSetup/EnvironmentSetup/createBodies.cpp
      Tudat/SimulationSetup/EnvironmentSetup/createBodies.h

* Add class `SphericalHarmonicsModelGravityFieldSettings` to load SH gravity files using the name of a SH model (enumeration)
      Tudat/SimulationSetup/EnvironmentSetup/createGravityField.cpp
      Tudat/SimulationSetup/EnvironmentSetup/createGravityField.h
      Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.cpp
      
* Add optional argument to function `getDefaultBodySettings` to be able to customise interpolation step for ephemeris (was fixed to 300.0 s before)
      Tudat/SimulationSetup/EnvironmentSetup/defaultBodies.cpp

* Move `createAccelerationModelsMap` function to avoid recursive file-inclusion cycle
      Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.cpp
      Tudat/SimulationSetup/PropagationSetup/createAccelerationModels.h
      Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.cpp
      Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h

* Move `createTorqueModelsMap` function to avoid recursive file-inclusion cycle
      Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.cpp
      Tudat/SimulationSetup/PropagationSetup/createNumericalSimulator.h
      Tudat/SimulationSetup/PropagationSetup/createTorqueModel.cpp
      Tudat/SimulationSetup/PropagationSetup/createTorqueModel.h
      
* Add property `componentIndex_` to `SingleDependentVariableSaveSettings` to allow saving just one of the elements of the vector
      Tudat/SimulationSetup/PropagationSetup/propagationOutput.cpp
      Tudat/SimulationSetup/PropagationSetup/propagationOutput.h
      
* Move `getVectorDependentVariableFunction` above `getDoubleDependentVariableFunction` 
      Tudat/SimulationSetup/PropagationSetup/propagationOutput.h
      
* Create class `VariableSettings` to represent Tudat variables and make `SingleDependentVariableSaveSettings` derive from it
      Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.cpp
      Tudat/SimulationSetup/PropagationSetup/propagationOutputSettings.h
      
* Enable constructing `PropagatorSettings` using a map of `AccelerationSettings`
      Tudat/SimulationSetup/PropagationSetup/propagationSettings.h
      
* Typedef for `std::map< std::string, std::map< std::string, std::vector< boost::shared_ptr< TorqueSettings > > > >` as `SelectedTorqueMap`
      Tudat/SimulationSetup/PropagationSetup/torqueSettings.h
      
* Include `mapTextFileReader` in `tudatSimulationHeader`
      Tudat/SimulationSetup/tudatSimulationHeader.h