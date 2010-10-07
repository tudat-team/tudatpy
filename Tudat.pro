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
    Astrodynamics/Bodies/Vehicles/vehicle.cpp \
    Astrodynamics/Bodies/Vehicles/vehicleExternalModel.cpp \
    Astrodynamics/EnvironmentModels/gravityFieldParameters.cpp \
    Astrodynamics/ForceModels/gravity.cpp \
    Astrodynamics/Propagators/propagator.cpp \
    Astrodynamics/Propagators/numericalPropagator.cpp \
    Astrodynamics/Propagators/bodyContainer.cpp \
    Basics/basicFunctions.cpp \
    Mathematics/basicMathematicsFunctions.cpp \
    Mathematics/GeometricShapes/surfaceGeometry.cpp \
    Mathematics/GeometricShapes/sphereSegment.cpp \
    Mathematics/GeometricShapes/geometricShape.cpp \
    Mathematics/LinearAlgebra/linearAlgebra.cpp \
    Mathematics/NumericalIntegration/rungeKutta4thOrderFixedStepsize.cpp \
    Mathematics/NumericalIntegration/integrator.cpp \
    Mathematics/NumericalIntegration/euler.cpp \
    Output/writingOutputToFile.cpp \
    earthOrbitExample.cpp
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
    Astrodynamics/Propagators/propagator.h \
    Astrodynamics/Propagators/numericalPropagator.h \
    Astrodynamics/Propagators/bodyContainer.h \
    Basics/basicFunctions.h \
    Mathematics/unitConversions.h \
    Mathematics/basicMathematicsFunctions.h \
    Mathematics/GeometricShapes/surfaceGeometry.h \
    Mathematics/GeometricShapes/sphereSegment.h \
    Mathematics/GeometricShapes/geometricShape.h \
    Mathematics/LinearAlgebra/linearAlgebra.h \
    Mathematics/NumericalIntegration/singleStepIntegrationMethods.h \
    Mathematics/NumericalIntegration/rungeKutta4thOrderFixedStepsize.h \
    Mathematics/NumericalIntegration/integrator.h \
    Mathematics/NumericalIntegration/euler.h \
    Output/writingOutputToFile.h \
    Output/outputHandling.h
INCLUDEPATH += External/Eigen-2.0.15 \
    Astrodynamics \
    Astrodynamics/Bodies \
    Astrodynamics/Bodies/CelestialBodies \
    Astrodynamics/Bodies/Vehicles \
    Astrodynamics/EnvironmentModels \
    Astrodynamics/ForceModels \
    Astrodynamics/Propagators \
    Basics \
    Mathematics \
    Mathematics/GeometricShapes \
    Mathematics/LinearAlgebra \
    Mathematics/NumericalIntegration \
    Output
