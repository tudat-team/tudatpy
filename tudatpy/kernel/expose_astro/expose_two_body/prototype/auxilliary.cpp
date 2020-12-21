
#include <math.h>

namespace tudat {
namespace prototype {

double calculateSphereOfInfluence(const double distanceToCentralBody,
                                  const double ratioOfOrbitingToCentralBodyMass) {
  // Return the radius of the sphere of influence.
  return distanceToCentralBody * std::pow(ratioOfOrbitingToCentralBodyMass, 0.4);
}

double calculateSphereOfInfluence(const double distanceToCentralBody,
                                  const double massOrbitingBody,
                                  const double massCentralBody) {
  // Return the radius of the sphere of influence.
  return calculateSphereOfInfluence(distanceToCentralBody, massOrbitingBody / massCentralBody);
}

}// namespace prototype
}// namespace tudat