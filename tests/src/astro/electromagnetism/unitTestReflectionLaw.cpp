/*    Copyright (c) 2010-2022, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MAIN

#include <iostream>

#include <boost/test/unit_test.hpp>

#include "tudat/basics/testMacros.h"
#include "tudat/astro/electromagnetism/reflectionLaw.h"
#include "tudat/math/basic/mathematicalConstants.h"


namespace tudat
{
namespace unit_tests
{

using namespace tudat::electromagnetism;
using Vector = Eigen::Vector3d;

BOOST_AUTO_TEST_SUITE(test_reflection_law)

BOOST_AUTO_TEST_CASE( testMirrorlikeReflection )
{
    const auto normal = Vector(0, 0, 1);

    TUDAT_CHECK_MATRIX_CLOSE(Vector(0, 0, 1), computeMirrorlikeReflection(Vector(0, 0, -1), normal), 1e-15);
    TUDAT_CHECK_MATRIX_CLOSE(Vector(-1, -1, 1), computeMirrorlikeReflection(Vector(-1, -1, -1), normal), 1e-15);
    TUDAT_CHECK_MATRIX_CLOSE(Vector(0, 0, 0), computeMirrorlikeReflection(Vector(-1, 0, 0), normal), 1e-15);
}



BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_ReflectedFraction_PointingAwayFromSource )
{
    const auto expectedReflectedFraction = 0;

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 1, 1).normalized();
    const auto observerDirection = Vector(1, 1, 1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(0.2, 0.4, 0.4, true);
    const auto actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_ReflectedFraction_PointingAwayFromObserver )
{
    const auto expectedReflectedFraction = 0;

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 1, -1).normalized();
    const auto observerDirection = Vector(1, 1, -1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(0.2, 0.4, 0.4, true);
    const auto actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_ReflectedFraction_PurelyAbsorbing )
{
    const auto expectedReflectedFraction = 0;

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 1, -1).normalized();
    const auto observerDirection = Vector(1, 1, 1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(1, 0, 0);
    const auto actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_ReflectedFraction_PurelyAbsorbing_InstantaneousReradiation )
{
    const auto expectedReflectedFraction = 1 / mathematical_constants::PI;

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 1, -1).normalized();
    const auto observerDirection = Vector(1, 1, 1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(1, 0, 0, true);
    const auto actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_ReflectedFraction_LambertianReflection )
{
    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 1, -1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(0.5, 0, 0.5);

    const auto expectedReflectedFraction = 0.5 / mathematical_constants::PI;
    const auto observerDirection = Vector(0, 0, 1).normalized();
    const auto actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_ReflectedFraction_MixedReflection_NotInMirroredPath )
{
    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 1, -1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(0, 0.5, 0.5);

    const auto expectedReflectedFraction = 0.5 / mathematical_constants::PI;
    const auto observerDirection = Vector(0, 0, 1).normalized();
    const auto actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_ReflectedFraction_MixedReflection_InMirroredPath )
{
    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 0, -1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(0, 0.5, 0.5);

    const auto expectedReflectedFraction = 0.5 / mathematical_constants::PI + 0.5 * sqrt(2);
    const auto observerDirection = Vector(1, 0, 1).normalized();
    const auto actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-10);
}



BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_Reaction_PointingAway )
{
    const auto expectedReaction = Vector::Zero();

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 1, 1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(0.2, 0.4, 0.4, true);
    const auto actualReaction = reflectionLaw.evaluateReactionVector(normal, incomingDirection);

    TUDAT_CHECK_MATRIX_CLOSE(actualReaction, expectedReaction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_Reaction_PurelyAbsorbing )
{
    const auto expectedReaction = Vector(-1, -1, -1).normalized();

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(-1, -1, -1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(1, 0, 0);
    const auto actualReaction = reflectionLaw.evaluateReactionVector(normal, incomingDirection);

    TUDAT_CHECK_MATRIX_CLOSE(actualReaction, expectedReaction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_Reaction_PurelyAbsorbing_InstantaneousReradiation )
{
    const auto expectedReaction = Vector(0, 0, -(1 + 2./3));

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(0, 0, -1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(1, 0, 0, true);
    const auto actualReaction = reflectionLaw.evaluateReactionVector(normal, incomingDirection);

    TUDAT_CHECK_MATRIX_CLOSE(actualReaction, expectedReaction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_Reaction_PurelySpecularReflection )
{
    const auto expectedReaction = Vector(0, 0, -2);

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(0, 0, -1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(0, 1, 0);
    const auto actualReaction = reflectionLaw.evaluateReactionVector(normal, incomingDirection);

    TUDAT_CHECK_MATRIX_CLOSE(actualReaction, expectedReaction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testSpecularDiffuseMix_Reaction_PurelyDiffuseReflection )
{
    const auto expectedReaction = Vector(0, 0, -(1 + 2./3));

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(0, 0, -1).normalized();
    SpecularDiffuseMixReflectionLaw reflectionLaw(0, 0, 1);
    const auto actualReaction = reflectionLaw.evaluateReactionVector(normal, incomingDirection);

    TUDAT_CHECK_MATRIX_CLOSE(actualReaction, expectedReaction, 1e-15);
}



BOOST_AUTO_TEST_CASE( testLambertian_ReflectedFraction )
{
    // Similar to testSpecularDiffuseMix_ReflectedFraction_LambertianReflection
    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(1, 1, -1).normalized();
    LambertianReflectionLaw reflectionLaw(0.5);

    const auto observerDirection = Vector(0, 0, 1).normalized();

    auto expectedReflectedFraction = 0.5 / mathematical_constants::PI;
    auto actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-15);

    reflectionLaw.setDiffuseReflectivity(0.8);

    expectedReflectedFraction = 0.8 / mathematical_constants::PI;
    actualReflectedFraction =
            reflectionLaw.evaluateReflectedFraction(normal, incomingDirection, observerDirection);

    BOOST_CHECK_CLOSE(actualReflectedFraction, expectedReflectedFraction, 1e-15);
}

BOOST_AUTO_TEST_CASE( testLambertian_Reaction )
{
    // Similar to testSpecularDiffuseMix_Reaction_PurelyDiffuseReflection
    const auto expectedReaction = Vector(0, 0, -(1 + 2./3));

    const auto normal = Vector(0, 0, 1);
    const auto incomingDirection = Vector(0, 0, -1).normalized();
    LambertianReflectionLaw reflectionLaw(1);
    const auto actualReaction = reflectionLaw.evaluateReactionVector(normal, incomingDirection);

    TUDAT_CHECK_MATRIX_CLOSE(actualReaction, expectedReaction, 1e-15);
}


BOOST_AUTO_TEST_SUITE_END()

} // namespace unit_tests
} // namespace tudat
