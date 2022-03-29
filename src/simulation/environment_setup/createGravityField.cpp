/*    Copyright (c) 2010-2019, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#include <boost/make_shared.hpp>

#include "tudat/interface/spice/spiceInterface.h"
#include "tudat/astro/gravitation/timeDependentSphericalHarmonicsGravityField.h"
#include "tudat/astro/gravitation/triAxialEllipsoidGravity.h"
#include "tudat/simulation/environment_setup/createGravityField.h"
#include "tudat/io/basicInputOutput.h"

namespace tudat
{

namespace simulation_setup
{

//! Get the path of the SH file for a SH model.
std::string getPathForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel )
{
    switch ( sphericalHarmonicsModel )
    {
    case egm96:
        return paths::getGravityModelsPath( ) + "/Earth/egm96.txt";
    case ggm02c:
        return paths::getGravityModelsPath( ) + "/Earth/ggm02c.txt";
    case ggm02s:
        return paths::getGravityModelsPath( ) + "/Earth/ggm02s.txt";
    case glgm3150:
        return paths::getGravityModelsPath( ) + "/Moon/glgm3150.txt";
    case lpe200:
        return paths::getGravityModelsPath( ) + "/Moon/lpe200.txt";
    case jgmro120d:
        return paths::getGravityModelsPath( ) + "/Mars/jgmro120d.txt";
    default:
        std::cerr << "No path known for Spherical Harmonics Model " << sphericalHarmonicsModel << std::endl;
        throw;
    }
}

//! Get the associated reference frame for a SH model.
std::string getReferenceFrameForSphericalHarmonicsModel( const SphericalHarmonicsModel sphericalHarmonicsModel )
{
    switch ( sphericalHarmonicsModel )
    {
    case egm96:
    case ggm02c:
    case ggm02s:
        return "IAU_Earth";
    case glgm3150:
    case lpe200:
        return "IAU_Moon";
    case jgmro120d:
        return "IAU_Mars";
    default:
        std::cerr << "No reference frame known for Spherical Harmonics Model " << sphericalHarmonicsModel << std::endl;
        throw;
    }
}

//! Constructor with custom model.
FromFileSphericalHarmonicsGravityFieldSettings::FromFileSphericalHarmonicsGravityFieldSettings(
        const std::string& filePath, const std::string& associatedReferenceFrame,
        const int maximumDegree, const int maximumOrder,
        const int gravitationalParameterIndex, const int referenceRadiusIndex,
        const double gravitationalParameter, const double referenceRadius ) :
    SphericalHarmonicsGravityFieldSettings( gravitationalParameter, referenceRadius, Eigen::MatrixXd( ),
                                            Eigen::MatrixXd( ), associatedReferenceFrame ),
    filePath_( filePath ),
    maximumDegree_( maximumDegree ),
    maximumOrder_( maximumOrder ),
    gravitationalParameterIndex_( gravitationalParameterIndex ),
    referenceRadiusIndex_( referenceRadiusIndex )
{
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients;
    std::pair< double, double > referenceData =
            readGravityFieldFile( filePath, maximumDegree, maximumOrder, coefficients,
                                  gravitationalParameterIndex, referenceRadiusIndex );
    gravitationalParameter_ = gravitationalParameterIndex >= 0 ? referenceData.first : gravitationalParameter;
    referenceRadius_ = referenceRadiusIndex >= 0 ? referenceData.second : referenceRadius;
    cosineCoefficients_ = coefficients.first;
    sineCoefficients_ = coefficients.second;
}

//! Constructor with model included in Tudat.
FromFileSphericalHarmonicsGravityFieldSettings::FromFileSphericalHarmonicsGravityFieldSettings(
        const SphericalHarmonicsModel sphericalHarmonicsModel ) :
    FromFileSphericalHarmonicsGravityFieldSettings( getPathForSphericalHarmonicsModel( sphericalHarmonicsModel ),
                                                    getReferenceFrameForSphericalHarmonicsModel( sphericalHarmonicsModel ),
                                                    50, 50, 0, 1 )
{
    sphericalHarmonicsModel_ = sphericalHarmonicsModel;
}


//! Function to read a gravity field file
std::pair< double, double  > readGravityFieldFile(
        const std::string& fileName, const int maximumDegree, const int maximumOrder,
        std::pair< Eigen::MatrixXd, Eigen::MatrixXd >& coefficients,
        const int gravitationalParameterIndex, const int referenceRadiusIndex )
{
    // Attempt to open gravity file.
    std::fstream stream( fileName.c_str( ), std::ios::in );
    if( stream.fail( ) )
    {
        throw std::runtime_error( "Pds gravity field data file could not be opened: " + fileName );
    }

    // Declare variables for reading file.
    std::vector< std::string > vectorOfIndividualStrings;
    vectorOfIndividualStrings.resize( 4 );
    std::string line;


    double gravitationalParameter = TUDAT_NAN;
    double referenceRadius = TUDAT_NAN;

    if( ( gravitationalParameterIndex >= 0 ) &&
            ( referenceRadiusIndex >= 0 ) )
    {
        // Get first line of file.
        std::getline( stream, line );

        // Get reference radius and gravitational parameter from first line of file.
        boost::algorithm::trim( line );
        boost::algorithm::split( vectorOfIndividualStrings,
                                 line,
                                 boost::algorithm::is_any_of( "\t, " ),
                                 boost::algorithm::token_compress_on );
        if( gravitationalParameterIndex >= static_cast< int >( vectorOfIndividualStrings.size( ) ) ||
                referenceRadiusIndex >= static_cast< int >( vectorOfIndividualStrings.size( ) ) )
        {
            throw std::runtime_error( "Error when reading gravity field file, requested header index exceeds file contents" );
        }

        gravitationalParameter = std::stod( vectorOfIndividualStrings[ gravitationalParameterIndex ] );
        referenceRadius = std::stod( vectorOfIndividualStrings[ referenceRadiusIndex ] );
    }
    else if( ( !( gravitationalParameterIndex >= 0 ) &&
               ( referenceRadiusIndex >= 0 ) ) ||
             ( ( gravitationalParameterIndex >= 0 ) &&
               !( referenceRadiusIndex >= 0 ) ) )
    {
        throw std::runtime_error( "Error when reading gravity field file, must retrieve either both or neither of Re and mu" );
    }


    // Declare variables for reading in cosine and sine coefficients.
    int currentDegree = 0, currentOrder = 0;
    Eigen::MatrixXd cosineCoefficients = Eigen::MatrixXd( maximumDegree + 1, maximumOrder + 1 );
    cosineCoefficients.setZero( );
    Eigen::MatrixXd sineCoefficients = Eigen::MatrixXd( maximumDegree + 1, maximumOrder + 1 );
    sineCoefficients.setZero( );

    // Read coefficients up to required maximum degree and order.
    while ( !stream.fail( ) && !stream.eof( ) &&
            ( currentDegree <= maximumDegree && currentOrder < maximumOrder )  )
    {
        // Read current line
        std::getline( stream, line );

        // Trim input string (removes all leading and trailing whitespaces).
        boost::algorithm::trim( line );

        // Split string into multiple strings, each containing one element from a line from the
        // data file.
        boost::algorithm::split( vectorOfIndividualStrings,
                                 line,
                                 boost::algorithm::is_any_of( ", " ),
                                 boost::algorithm::token_compress_on );

        // Check current line for consistency
        if( vectorOfIndividualStrings.size( ) != 0 )
        {
            if( vectorOfIndividualStrings.size( ) < 4 )
            {
                std::string errorMessage = "Error when reading pds gravity field file, number of fields is " +
                        std::to_string( vectorOfIndividualStrings.size( ) );
                throw std::runtime_error( errorMessage );
            }
            else
            {
                // Read current degree and orde from line.
                currentDegree = std::stoi( vectorOfIndividualStrings[ 0 ] );
                currentOrder = std::stoi( vectorOfIndividualStrings[ 1 ] );

                // Set cosine and sine coefficients for current degree and order.
                if( currentDegree <= maximumDegree && currentOrder <= maximumOrder )
                {
                    cosineCoefficients( currentDegree, currentOrder ) =
                            std::stod( vectorOfIndividualStrings[ 2 ] );
                    sineCoefficients( currentDegree, currentOrder ) =
                            std::stod( vectorOfIndividualStrings[ 3 ] );
                }
            }
        }
    }

    // Set cosine coefficient at (0,0) to 1.
    cosineCoefficients( 0, 0 ) = 1.0;
    coefficients = std::make_pair( cosineCoefficients, sineCoefficients );

    return std::make_pair( gravitationalParameter, referenceRadius );
}

//! Function to create a gravity field model.
std::shared_ptr< gravitation::GravityFieldModel > createGravityFieldModel(
        const std::shared_ptr< GravityFieldSettings > gravityFieldSettings,
        const std::string& body,
        const SystemOfBodies& bodies,
        const std::vector< std::shared_ptr< GravityFieldVariationSettings > >& gravityFieldVariationSettings )
{
    using namespace tudat::gravitation;

    // Declare return object.
    std::shared_ptr< GravityFieldModel > gravityFieldModel;

    // Check which type of gravity field model is to be created.
    switch( gravityFieldSettings->getGravityFieldType( ) )
    {
    case central:
    {
        // Check whether settings for point mass gravity field model are consistent with its type.
        std::shared_ptr< CentralGravityFieldSettings > centralFieldSettings =
                std::dynamic_pointer_cast< CentralGravityFieldSettings >( gravityFieldSettings );
        if( centralFieldSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected central field settings when making gravity field model for body " +
                        body);
        }
        else if( gravityFieldVariationSettings.size( ) != 0 )
        {
            throw std::runtime_error( "Error, requested central gravity field, but field variations settings are not empty." );
        }
        else
        {
            // Create and initialize point mass gravity field model.
            gravityFieldModel = std::make_shared< GravityFieldModel >(
                        centralFieldSettings->getGravitationalParameter( ) );
        }
        break;
    }
    case central_spice:
    {
        if( gravityFieldVariationSettings.size( ) != 0 )
        {
            throw std::runtime_error( "Error, requested central gravity field, but field variations settings are not empty." );
        }
        else
        {
            // Create and initialize point mass gravity field model from Spice.
            gravityFieldModel = std::make_shared< GravityFieldModel >(
                        spice_interface::getBodyGravitationalParameter( body ) );
        }

        break;
    }
    case spherical_harmonic:
    {
        // Check whether settings for spherical harmonic gravity field model are consistent with
        // its type.
        std::shared_ptr< SphericalHarmonicsGravityFieldSettings > sphericalHarmonicFieldSettings =
                std::dynamic_pointer_cast< SphericalHarmonicsGravityFieldSettings >(
                    gravityFieldSettings );

        if( sphericalHarmonicFieldSettings == nullptr )
        {
            throw std::runtime_error(
                        "Error, expected spherical harmonic field settings when making gravity field model of "
                        + body );
        }
        else
        {
            std::function< void( ) > inertiaTensorUpdateFunction;
            if( bodies.count( body ) == 0 )
            {
                inertiaTensorUpdateFunction = std::function< void( ) >( );
            }
            else
            {
                inertiaTensorUpdateFunction =
                    std::bind( &Body::setBodyInertiaTensorFromGravityFieldAndExistingMeanMoment, bodies.at( body ), true );
                if( sphericalHarmonicFieldSettings->getScaledMeanMomentOfInertia( ) == sphericalHarmonicFieldSettings->getScaledMeanMomentOfInertia( ) )
                {
                    bodies.at( body )->setBodyInertiaTensor(
                                sphericalHarmonicFieldSettings->getInertiaTensor( ),
                                sphericalHarmonicFieldSettings->getScaledMeanMomentOfInertia( )) ;

                }
            }

            // Check consistency of cosine and sine coefficients.
            if( ( sphericalHarmonicFieldSettings->getCosineCoefficients( ).rows( ) !=
                  sphericalHarmonicFieldSettings->getSineCoefficients( ).rows( ) ) ||
                    ( sphericalHarmonicFieldSettings->getCosineCoefficients( ).cols( ) !=
                      sphericalHarmonicFieldSettings->getSineCoefficients( ).cols( ) ) )
            {
                throw std::runtime_error(
                            std::string( "Error when making spherical harmonic field, sine and " ) +
                            std::string( "cosine matrix  sizes are not equal for body " ) + body );
            }
            else
            {

                if( gravityFieldVariationSettings.size( ) == 0 &&
                        sphericalHarmonicFieldSettings->getCreateTimeDependentField( ) == 0 )
                {
                    // Create and initialize spherical harmonic gravity field model.
                    gravityFieldModel = std::make_shared< SphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ),
                                sphericalHarmonicFieldSettings->getAssociatedReferenceFrame( ),
                                inertiaTensorUpdateFunction );
                }
                else
                {
                    if( bodies.at( body )->getGravityFieldModel( ) != nullptr )
                    {
                        std::string errorMessage = "Warning when making time-dependent gravity field model for body " + body +
                                " existing gravity field is not empty but overwritten in Body! ";
                        throw std::runtime_error( errorMessage );
                    }

                    // Create preliminary TimeDependentSphericalHarmonicsGravityField, without actual variation settings.
                    gravityFieldModel = std::make_shared< TimeDependentSphericalHarmonicsGravityField >(
                                sphericalHarmonicFieldSettings->getGravitationalParameter( ),
                                sphericalHarmonicFieldSettings->getReferenceRadius( ),
                                sphericalHarmonicFieldSettings->getCosineCoefficients( ),
                                sphericalHarmonicFieldSettings->getSineCoefficients( ),
                                sphericalHarmonicFieldSettings->getAssociatedReferenceFrame( ),
                                inertiaTensorUpdateFunction );
                }


            }
        }
        break;
    }
    default:
        throw std::runtime_error(
                    "Error, did not recognize gravity field model settings type " +
                    std::to_string(
                        gravityFieldSettings->getGravityFieldType( ) ) );
    }

    return gravityFieldModel;
}

//! Function to create gravity field settings for a homogeneous triaxial ellipsoid
std::shared_ptr< SphericalHarmonicsGravityFieldSettings > createHomogeneousTriAxialEllipsoidGravitySettings(
        const double axisA, const double axisB, const double axisC, const double ellipsoidDensity,
        const int maximumDegree, const int maximumOrder,
        const std::string& associatedReferenceFrame,
        const double gravitationalConstant)
{
    // Compute reference quantities
    double ellipsoidGravitationalParameter = gravitation::calculateTriAxialEllipsoidVolume(
                axisA, axisB, axisC ) * ellipsoidDensity * gravitationalConstant;
    double ellipsoidReferenceRadius = gravitation::calculateTriAxialEllipsoidReferenceRadius(
                axisA, axisB, axisC );

    // Compute coefficients
    std::pair< Eigen::MatrixXd, Eigen::MatrixXd > coefficients =
            gravitation::createTriAxialEllipsoidNormalizedSphericalHarmonicCoefficients(
                axisA, axisB, axisC, maximumDegree, maximumOrder );

    return std::make_shared< SphericalHarmonicsGravityFieldSettings >(
                ellipsoidGravitationalParameter, ellipsoidReferenceRadius, coefficients.first,
                coefficients.second, associatedReferenceFrame );
}

PolyhedronGravityFieldSettings::PolyhedronGravityFieldSettings (
        const double gravitationalConstant,
        const double density,
        const Eigen::MatrixXd& verticesCoordinates,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const std::string& associatedReferenceFrame):
    GravityFieldSettings( polyhedron ),
    gravitationalConstantTimesDensity_( gravitationalConstant * density ),
    gravitationalParameter_( TUDAT_NAN ),
    volume_( TUDAT_NAN ),
    verticesCoordinates_( verticesCoordinates ),
    verticesDefiningEachFacet_( verticesDefiningEachFacet ),
    associatedReferenceFrame_( associatedReferenceFrame )
{
    const unsigned int numberOfVertices = verticesCoordinates_.rows();
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();

    if ( numberOfFacets != 2 * ( numberOfVertices - 2) )
    {
        throw std::runtime_error( "Number of polyhedron facets and vertices not consistent." );
    }

    // Compute edges in polyhedron
    computeVerticesAndFacetsDefiningEachEdge();

    // Compute facet dyads
    computeFacetNormalsAndDyads();

    // Compute edge dyads
    computeEdgeDyads();
}

PolyhedronGravityFieldSettings::PolyhedronGravityFieldSettings (
        const double gravitationalParameter,
        const Eigen::MatrixXd& verticesCoordinates,
        const Eigen::MatrixXi& verticesDefiningEachFacet,
        const std::string& associatedReferenceFrame):
        GravityFieldSettings( polyhedron ),
        gravitationalParameter_( gravitationalParameter ),
        verticesCoordinates_( verticesCoordinates ),
        verticesDefiningEachFacet_( verticesDefiningEachFacet ),
        associatedReferenceFrame_( associatedReferenceFrame )
{
    const unsigned int numberOfVertices = verticesCoordinates_.rows();
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();

    if ( numberOfFacets != 2 * ( numberOfVertices - 2) )
    {
        throw std::runtime_error( "Number of polyhedron facets and vertices not consistent." );
    }

    // Compute edges in polyhedron
    computeVerticesAndFacetsDefiningEachEdge();

    // Compute facet dyads
    computeFacetNormalsAndDyads();

    // Compute edge dyads
    computeEdgeDyads();

    // Compute volume
    computeVolume();

    // Define value of gravitational constant times the density
    gravitationalConstantTimesDensity_ = gravitationalParameter_ / volume_;
}

void PolyhedronGravityFieldSettings::computeVerticesAndFacetsDefiningEachEdge ( )
{
    const unsigned int numberOfVertices = verticesCoordinates_.rows();
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();
    const unsigned int numberOfEdges = 3 * ( numberOfVertices - 2 );

    verticesDefiningEachEdge_ = Eigen::MatrixXi::Constant( numberOfEdges, 2, TUDAT_NAN );
    facetsDefiningEachEdge_ = Eigen::MatrixXi::Constant( numberOfEdges, 2, TUDAT_NAN );

    unsigned int numberOfInsertedEdges = 0;
    for ( unsigned int facet = 0; facet < numberOfFacets; ++facet )
    {
        // Extract edges from face
        std::vector< std::vector< int > > edgesToInsert;
        const int vertex0 = verticesDefiningEachFacet_(facet,0);
        const int vertex1 = verticesDefiningEachFacet_(facet,1);
        const int vertex2 = verticesDefiningEachFacet_(facet,2);
        edgesToInsert.push_back( { vertex0, vertex1 } );
        edgesToInsert.push_back( { vertex1, vertex2 } );
        edgesToInsert.push_back( { vertex2, vertex0 } );

        // Check if all edges of the facet have been included in verticesDefiningEachEdge. If so, remove the edge from
        // the edgesToInsert vector and add the facet to facetsDefiningEachEdge
        for ( unsigned int edge = 0; edge < numberOfInsertedEdges and !edgesToInsert.empty(); ++edge )
        {
            // Loop over edges of current facet still to be inserted
            for ( unsigned int i = 0; i < edgesToInsert.size(); ++i)
            {
                if ( ( edgesToInsert.at(i).at(0) == verticesDefiningEachEdge_(edge, 0) && edgesToInsert.at(i).at(1) == verticesDefiningEachEdge_(edge, 1) ) ||
                        ( edgesToInsert.at(i).at(0) == verticesDefiningEachEdge_(edge, 1) && edgesToInsert.at(i).at(1) == verticesDefiningEachEdge_(edge, 0) ) )
                {
                    edgesToInsert.erase(edgesToInsert.begin() + i);
                    facetsDefiningEachEdge_(edge, 1) = facet;
                }
            }
        }

        // Check if any of the facet's edges still needs to be added to verticesDefiningEachEdge, and add it/them if so
        for ( unsigned int i = 0; i < edgesToInsert.size(); ++i)
        {
            verticesDefiningEachEdge_(numberOfInsertedEdges,0) = edgesToInsert.at(i).at(0);
            verticesDefiningEachEdge_(numberOfInsertedEdges,1) = edgesToInsert.at(i).at(1);
            facetsDefiningEachEdge_(numberOfInsertedEdges, 0) = facet;
            ++numberOfInsertedEdges;
        }
    }

    // Sanity checks
    if ( numberOfInsertedEdges != numberOfEdges )
    {
        throw std::runtime_error( "Extracted number of polyhedron edges not correct." );
    }
    for ( unsigned int i = 0; i < numberOfEdges; ++i )
    {
        for (unsigned int j : {0,1} )
        {
            if ( verticesDefiningEachEdge_(i,j) != verticesDefiningEachEdge_(i,j) )
            {
                throw std::runtime_error( "The vertices defining some edge were not selected." );
            }
            else if ( facetsDefiningEachEdge_(i,j) != facetsDefiningEachEdge_(i,j) )
            {
                throw std::runtime_error( "The facets defining some edge were not selected." );
            }
        }
    }
}

void PolyhedronGravityFieldSettings::computeFacetNormalsAndDyads ( )
{
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();

    facetNormalVectors_.resize(numberOfFacets);
    facetDyads_.resize(numberOfFacets);

    // Loop over facets, and for each facet compute the facet normal and the facet dyad
    for (unsigned int facet = 0; facet < numberOfFacets; ++facet)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,2),0);

        // Compute outward-pointing vector normal to facet
        facetNormalVectors_.at(facet) = (vertex1 - vertex0).cross(vertex2 - vertex1);
        facetNormalVectors_.at(facet).normalize();

        // Compute facet dyad (outer product)
        facetDyads_.at(facet) = facetNormalVectors_.at(facet) * facetNormalVectors_.at(facet).transpose();
    }

}

void PolyhedronGravityFieldSettings::computeEdgeDyads ( )
{
    const unsigned int numberOfEdges = verticesDefiningEachEdge_.rows();

    edgeDyads_.resize(numberOfEdges);

    // Loop over edges, and for each edge compute the edge dyad
    for (unsigned int edge = 0; edge < numberOfEdges; ++edge)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates_.block<1,3>(verticesDefiningEachEdge_(edge,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates_.block<1,3>(verticesDefiningEachEdge_(edge,1),0);

        // Compute edge normals, with arbitrary direction
        Eigen::Vector3d facetNormalFacetA = facetNormalVectors_.at(facetsDefiningEachEdge_(edge, 0) );
        Eigen::Vector3d edgeNormalFacetA = ( vertex0 - vertex1).cross(facetNormalFacetA );
        edgeNormalFacetA.normalize();
        Eigen::Vector3d facetNormalFacetB = facetNormalVectors_.at(facetsDefiningEachEdge_(edge, 1) );
        Eigen::Vector3d edgeNormalFacetB = ( vertex0 - vertex1).cross(facetNormalFacetB );
        edgeNormalFacetB.normalize();

        // Check order of vertices for the facet associated with edgeNormalA, and correct the edge normal directions based on that
        unsigned int vertex0FacetAIndex = TUDAT_NAN, vertex1FacetAIndex = TUDAT_NAN;
        for (unsigned int facetVertex = 0; facetVertex < 3; ++facetVertex )
        {
            if ( verticesDefiningEachFacet_(facetsDefiningEachEdge_(edge,0), facetVertex) == verticesDefiningEachEdge_(edge, 0) )
            {
                vertex0FacetAIndex = facetVertex;
            }
            else if ( verticesDefiningEachFacet_(facetsDefiningEachEdge_(edge,0), facetVertex) == verticesDefiningEachEdge_(edge, 1) )
            {
                vertex1FacetAIndex = facetVertex;
            }
        }

        if (( vertex0FacetAIndex == vertex1FacetAIndex + 1) || ( vertex0FacetAIndex == 0 && vertex1FacetAIndex == 2 ) )
        {
            edgeNormalFacetB = - edgeNormalFacetB;
        }
        else
        {
            edgeNormalFacetA = - edgeNormalFacetA;
        }

        // Compute edge dyads: outer product
        edgeDyads_.at(edge) = facetNormalFacetA * edgeNormalFacetA.transpose() + facetNormalFacetB * edgeNormalFacetB.transpose();

    }
}

void PolyhedronGravityFieldSettings::computeVolume ( )
{
    const unsigned int numberOfFacets = verticesDefiningEachFacet_.rows();

    // Initialize volume
    volume_ = 0;

    for (unsigned int facet = 0; facet < numberOfFacets; ++facet)
    {
        Eigen::Vector3d vertex0 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,0),0);
        Eigen::Vector3d vertex1 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,1),0);
        Eigen::Vector3d vertex2 = verticesCoordinates_.block<1,3>(verticesDefiningEachFacet_(facet,2),0);

        volume_ += 1.0/6.0 * vertex0.dot( vertex1.cross( vertex2 ) );
    }

}


} // namespace simulation_setup

} // namespace tudat
