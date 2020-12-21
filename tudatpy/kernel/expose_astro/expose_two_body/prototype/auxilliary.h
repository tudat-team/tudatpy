#ifndef TUDATPY_AUXILLIARY_H
#define TUDATPY_AUXILLIARY_H

namespace tudat {
namespace prototype {

double calculateSphereOfInfluence(const double distanceToCentralBody,
                                  const double ratioOfOrbitingToCentralBodyMass);

double calculateSphereOfInfluence(const double distanceToCentralBody,
                                  const double massOrbitingBody,
                                  const double massCentralBody);

}// namespace prototype
}// namespace tudat

#endif//TUDATPY_AUXILLIARY_H
