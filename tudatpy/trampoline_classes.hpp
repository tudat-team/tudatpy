//
// Created by ggarrett on 24-04-20.
//

#ifndef TUDATBUNDLE_TRAMPOLINE_CLASSES_HPP
#define TUDATBUNDLE_TRAMPOLINE_CLASSES_HPP

#include "Tudat/Astrodynamics/Ephemerides/ephemeris.h"
#include "Tudat/Basics/timeType.h"

#include "Tudat/Astrodynamics/Ephemerides/tabulatedEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/constantEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/keplerEphemeris.h"
#include "Tudat/Astrodynamics/Ephemerides/customEphemeris.h"

#include "Tudat/Astrodynamics/BasicAstrodynamics/accelerationModel.h" // For AccelerationModel

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>


namespace te = tudat::ephemerides;
namespace tba = tudat::basic_astrodynamics;

namespace tudatpy {

    namespace basic_astrodynamics {

        template<class AccelerationModelBase = tba::AccelerationModel<Eigen::Vector3d>>
        class PyAccelerationModel
                : public AccelerationModelBase {
        public :
            using AccelerationModelBase::AccelerationModelBase;

            Eigen::Vector3d getAcceleration() override {
                PYBIND11_OVERLOAD_PURE(Eigen::Vector3d, AccelerationModelBase, getAcceleration,)

            }
        };
    };

    namespace ephemerides {

        template<class EphemerisBase = te::Ephemeris>
        class PyEphemeris
                : public EphemerisBase {
        public:
            using EphemerisBase::EphemerisBase;

            Eigen::Vector6d getCartesianState(const double secondsSinceEpoch) override {
                PYBIND11_OVERLOAD_PURE(Eigen::Vector6d, EphemerisBase, getCartesianState, secondsSinceEpoch)
            }

            Eigen::Matrix<long double, 6, 1> getCartesianLongState(const double secondsSinceEpoch) override {
                PYBIND11_OVERLOAD_PURE(Eigen::Vector3d, EphemerisBase, getCartesianLongState, secondsSinceEpoch)
            }

            Eigen::Vector6d getCartesianStateFromExtendedTime(const tudat::Time &currentTime) override {
                PYBIND11_OVERLOAD_PURE(Eigen::Vector3d, EphemerisBase, getCartesianStateFromExtendedTime, currentTime)
            }

            Eigen::Matrix<long double, 6, 1>
            getCartesianLongStateFromExtendedTime(const tudat::Time &currentTime) override {
                PYBIND11_OVERLOAD_PURE(Eigen::Vector3d, EphemerisBase, getCartesianLongStateFromExtendedTime,
                                       currentTime)
            }


        };

//        template<class TabulatedCartesianEphemerisBase = te::TabulatedCartesianEphemeris<>>
//        class PyTabulatedCartesianEphemeris
//                : public PyEphemeris<TabulatedCartesianEphemerisBase> {
//
//            // Complete PyTabulatedCartesianEphemeris overrides and implement.
//        };
//
//        template<class ConstantEphemerisBase = te::ConstantEphemeris>
//        class PyConstantEphemeris
//                : public PyEphemeris<ConstantEphemerisBase> {
//        public:
//            using PyEphemeris<ConstantEphemerisBase>::PyEphemeris;
//
//            Eigen::Vector6d getCartesianState(const double secondsSinceEpoch) override {
//                PYBIND11_OVERLOAD_PURE(Eigen::Vector6d, ConstantEphemerisBase, getCartesianState, secondsSinceEpoch)
//            }
//
//            Eigen::Matrix<long double, 6, 1> getCartesianLongState(const double secondsSinceEpoch) override {
//                PYBIND11_OVERLOAD(Eigen::Vector3d, ConstantEphemerisBase, getCartesianLongState, secondsSinceEpoch);
//            }
//
//            Eigen::Vector6d getCartesianStateFromExtendedTime(const tudat::Time &currentTime) override {
//                PYBIND11_OVERLOAD(Eigen::Vector3d, ConstantEphemerisBase, getCartesianStateFromExtendedTime,
//                                  currentTime);
//            }
//
//            Eigen::Matrix<long double, 6, 1>
//            getCartesianLongStateFromExtendedTime(const tudat::Time &currentTime) override {
//                PYBIND11_OVERLOAD(Eigen::Vector3d, ConstantEphemerisBase, getCartesianLongStateFromExtendedTime,
//                                  currentTime);
//            }
//
////            void updateConstantState(const Eigen::Vector6d &newState) override {
////                PYBIND11_OVERLOAD(void, ConstantEphemerisBase, updateConstantState,
////                                  newState);
////            }
//
//
//        };
//
//        template<class KeplerEphemerisBase = te::KeplerEphemeris>
//        class PyKeplerEphemeris
//                : public PyEphemeris<KeplerEphemerisBase> {
//
//            // TODO: Complete PyKeplerEphemeris overrides and implement.
//        };
//
//        template<class CustomEphemerisBase = te::CustomEphemeris>
//        class PyCustomEphemeris
//                : public PyEphemeris<CustomEphemerisBase> {
//
//            // TODO: Complete PyCustomEphemeris overrides and implement.
//        };


    }

    namespace simulation_setup {


    }

}
#endif //TUDATBUNDLE_TRAMPOLINE_CLASSES_HPP
