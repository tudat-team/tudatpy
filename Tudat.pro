# -------------------------------------------------
# Project created by QtCreator 2010-09-30T01:39:33
# -------------------------------------------------
QT -= gui
TARGET = Tudat
CONFIG += console
CONFIG -= app_bundle
TEMPLATE = app
SOURCES += Astrodynamics/Bodies/body.cpp \
    Astrodynamics/Bodies/CelestialBodies/celestialBody.cpp \
    Astrodynamics/Bodies/CelestialBodies/predefinedPlanets.cpp \
    Astrodynamics/Bodies/Vehicles/vehicle.cpp \
    Astrodynamics/Bodies/Vehicles/vehicleExternalModel.cpp \
    Astrodynamics/EnvironmentModels/unitTestSphericalHarmonicsGravityField.cpp \
    Astrodynamics/EnvironmentModels/sphericalHarmonicsGravityField.cpp \
    Astrodynamics/EnvironmentModels/predefinedGravityFieldModels.cpp \
    Astrodynamics/EnvironmentModels/gravityFieldModel.cpp \
    Astrodynamics/ForceModels/gravity.cpp \
    Astrodynamics/physicalConstants.cpp \
    Astrodynamics/unitTestPhysicalConstants.cpp \
    Astrodynamics/Propagators/propagator.cpp \
    Astrodynamics/Propagators/numericalPropagator.cpp \
    Astrodynamics/Propagators/bodyContainer.cpp \
    Astrodynamics/States/state.cpp \
    Astrodynamics/States/orbitalElements.cpp \
    Astrodynamics/States/orbitalElementConversions.cpp \
    Astrodynamics/States/keplerianElements.cpp \
    Astrodynamics/States/cartesianElements.cpp \
    Astrodynamics/States/unitTestOrbitalElementConversions.cpp \
    Basics/basicFunctions.cpp \
    Mathematics/basicMathematicsFunctions.cpp \
    Mathematics/randomNumberGenerator.cpp \
    Mathematics/unitTestRandomNumberGenerator.cpp \
    Mathematics/GeometricShapes/surfaceGeometry.cpp \
    Mathematics/GeometricShapes/sphereSegment.cpp \
    Mathematics/GeometricShapes/geometricShape.cpp \
    Mathematics/LinearAlgebra/linearAlgebra.cpp \
    Mathematics/NumericalIntegration/rungeKutta4thOrderFixedStepsize.cpp \
    Mathematics/NumericalIntegration/integrator.cpp \
    Mathematics/NumericalIntegration/euler.cpp \
    Mathematics/unitTestUnitConversions.cpp \
    Output/writingOutputToFile.cpp
HEADERS += Astrodynamics/Bodies/body.h \
    Astrodynamics/Bodies/CelestialBodies/predefinedPlanets.h \
    Astrodynamics/Bodies/CelestialBodies/celestialBody.h \
    Astrodynamics/Bodies/Vehicles/vehicle.h \
    Astrodynamics/Bodies/Vehicles/vehicleExternalModel.h \
    Astrodynamics/EnvironmentModels/predefinedGravityFieldParameters.h \
    Astrodynamics/EnvironmentModels/gravityFieldParameters.h \
    Astrodynamics/ForceModels/gravity.h \
    Astrodynamics/ForceModels/forceModel.h \
    Astrodynamics/physicalConstants.h \
    Astrodynamics/unitTestPhysicalConstants.h \
    Astrodynamics/Propagators/propagator.h \
    Astrodynamics/Propagators/numericalPropagator.h \
    Astrodynamics/Propagators/bodyContainer.h \
    Astrodynamics/States/state.h \
    Astrodynamics/States/orbitalElements.h \
    Astrodynamics/States/orbitalElementConversions.h \
    Astrodynamics/States/keplerianElements.h \
    Astrodynamics/States/cartesianElements.h \
    Astrodynamics/States/unitTestOrbitalElementConversions.h \
    Basics/basicFunctions.h \
    Mathematics/unitConversions.h \
    Mathematics/unitTestUnitConversions.h \
    Mathematics/basicMathematicsFunctions.h \
    Mathematics/randomNumberGenerator.h \
    Mathematics/unitTestRandomNumberGenerator.h \
    Mathematics/GeometricShapes/surfaceGeometry.h \
    Mathematics/GeometricShapes/sphereSegment.h \
    Mathematics/GeometricShapes/geometricShape.h \
    Mathematics/LinearAlgebra/linearAlgebra.h \
    Mathematics/NumericalIntegration/singleStepIntegrationMethods.h \
    Mathematics/NumericalIntegration/rungeKutta4thOrderFixedStepsize.h \
    Mathematics/NumericalIntegration/integrator.h \
    Mathematics/NumericalIntegration/euler.h \
    Output/outputHandling.h \
    Output/writingOutputToFile.h
INCLUDEPATH += External/Eigen-2.0.15 \
    Astrodynamics \
    Astrodynamics/Bodies \
    Astrodynamics/Bodies/CelestialBodies \
    Astrodynamics/Bodies/Vehicles \
    Astrodynamics/EnvironmentModels \
    Astrodynamics/ForceModels \
    Astrodynamics/Propagators \
    Astrodynamics/States \
    Basics \
    Mathematics \
    Mathematics/GeometricShapes \
    Mathematics/LinearAlgebra \
    Mathematics/NumericalIntegration \
    Output
