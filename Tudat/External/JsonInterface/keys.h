/*    Copyright (c) 2010-2017, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDAT_JSONINTERFACE_KEYS_H
#define TUDAT_JSONINTERFACE_KEYS_H

#include <string>
#include <vector>

namespace tudat
{

namespace json_interface
{

//! Keys recognised by json_interface.
struct Keys
{
    static const std::string simulation;
    struct Simulation
    {
        static const std::string startEpoch;
        static const std::string endEpoch;
        static const std::string globalFrameOrigin;
        static const std::string globalFrameOrientation;
        static const std::string spiceKernels;
        static const std::string preloadSpiceData;
    };

    static const std::string bodies;
    struct Body
    {
        static const std::string useDefaultSettings;
        static const std::string mass;
        static const std::string referenceArea;

        static const std::string aerodynamics;
        struct Aerodynamics
        {
            static const std::string type;
            static const std::string referenceArea;
            static const std::string dragCoefficient;
            static const std::string forceCoefficients;
            static const std::string momentCoefficients;
            static const std::string areCoefficientsInAerodynamicFrame;
            static const std::string areCoefficientsInNegativeAxisDirection;
        };

        static const std::string atmosphere;
        struct Atmosphere
        {
            static const std::string type;
            static const std::string densityScaleHeight;
            static const std::string constantTemperature;
            static const std::string densityAtZeroAltitude;
            static const std::string specificGasConstant;
            static const std::string atmosphereFile;
            static const std::string spaceWeatherFile;
        };

        static const std::string ephemeris;
        struct Ephemeris
        {
            static const std::string type;
            static const std::string frameOrigin;
            static const std::string frameOrientation;
            static const std::string makeMultiArc;
            static const std::string correctForStellarAbberation;
            static const std::string correctForLightTimeAbberation;
            static const std::string convergeLighTimeAbberation;
            static const std::string initialTime;
            static const std::string finalTime;
            static const std::string timeStep;
            static const std::string interpolator;
            static const std::string useLongDoubleStates;
            static const std::string bodyIdentifier;
            static const std::string useCircularCoplanarApproximation;
            static const std::string constantState;
            // static const std::string customStateFunction;
            static const std::string initialStateInKeplerianElements;
            static const std::string epochOfInitialState;
            static const std::string centralBodyGravitationalParameter;
            static const std::string rootFinderAbsoluteTolerance;
            static const std::string rootFinderMaximumNumberOfIterations;
            static const std::string bodyStateHistory;
        };

        static const std::string gravityField;
        struct GravityField
        {
            static const std::string type;
            static const std::string gravitationalParameter;
            static const std::string referenceRadius;
            static const std::string cosineCoefficients;
            static const std::string sineCoefficients;
            static const std::string associatedReferenceFrame;
            static const std::string sphericalHarmonicsFile;
            static const std::string maximumDegree;
            static const std::string maximumOrder;
            static const std::string gravitationalParameterIndex;
            static const std::string referenceRadiusIndex;
        };

        static const std::string radiationPressure;
        struct RadiationPressure
        {
            static const std::string type;
            static const std::string referenceArea;
            static const std::string radiationPressureCoefficient;
            static const std::string ocultingBodies;
        };
    };

    static const std::string propagators;
    struct Propagator
    {
        static const std::string integratedStateType;
        static const std::string type;
        static const std::string centralBodies;
        static const std::string bodiesToPropagate;
        static const std::string initialStates;
        static const std::string initialStateTypes;

        static const std::string termination;
        struct Termination
        {

        };

        static const std::string accelerations;
        struct Acceleration
        {

        };
    };

    static const std::string integrator;
    struct Integrator
    {
        static const std::string type;
        static const std::string initialTime;
        static const std::string initialTimeStep;
        static const std::string saveFrequency;
        static const std::string rungeKuttaCoefficientSet;
        static const std::string minimumStepSize;
        static const std::string maximumStepSize;
        static const std::string relativeErrorTolerance;
        static const std::string absoluteErrorTolerance;
        static const std::string safetyFactorForNextStepSize;
        static const std::string maximumFactorIncreaseForNextStepSize;
        static const std::string minimumFactorDecreaseForNextStepSize;
    };

    static const std::string output;
};

// FIXME: what about arrays?
//! Class for specifying a key tree (key.subkey.subsubkey ...) used to access data from `json` objects.
class KeyTree : public std::vector< std::string >
{
public:
    //! Inherit constructors.
    using vector< std::string >::vector;

    //! Constructor with a single string key.
    /*!
     * Constructor with a single string key.
     * \param key The key to be accessed.
     */
    KeyTree( const std::string& key ) : KeyTree( std::vector< std::string >( { key } ) ) { }

    //! Constructor with a single char key.
    /*!
     * Constructor with a single char key.
     * \param key The key to be accessed.
     */
    //! Constructor with a single char key.
    KeyTree( const char* key ) : KeyTree( std::string( key ) ) { }

    //! -DOC
    KeyTree& operator+ ( const KeyTree& subkeys ) const
    {
        KeyTree* compoundKeyTree = new KeyTree( *this );
        for ( std::string subkey : subkeys )
        {
            compoundKeyTree->push_back( subkey );
        }
        return *compoundKeyTree;
    }
};

//! String representation for `KeyTree`, as key.subkey.subsubkey ...
inline std::ostream& operator<< ( std::ostream & stringRepresentation, KeyTree const & keyTree )
{
    for ( unsigned int i = 0; i < keyTree.size(); ++i )
    {
        stringRepresentation << keyTree.at( i );
        if ( i < keyTree.size() - 1 )
        {
            stringRepresentation << ".";
        }
    }
    return stringRepresentation;
}

//! Key trees recognised by `json_interface`.
struct KeyTrees
{
    struct Simulation
    {
        static const KeyTree startEpoch;
        static const KeyTree endEpoch;
        static const KeyTree globalFrameOrigin;
        static const KeyTree globalFrameOrientation;
        static const KeyTree spiceKernels;
        static const KeyTree preloadSpiceData;
    };

    /*
    struct Integrator
    {
        static const KeyTree initialTime;
    };
    */
};

} // namespace json_interface

} // namespace tudat

#endif // TUDAT_JSONINTERFACE_KEYS_H
