/*! \file writingOutputToFile.h
 *    Header file that defines the class containing all funtionality pertaining
 *    to writing output to file included in Tudat.
 *
 *    Path              : /Output/
 *    Version           : 5
 *    Check status      : Checked
 *
 *    Author            : K. Kumar
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : K.Kumar@tudelft.nl
 *
 *    Author            : D. Dirkx
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : d.dirkx@tudelft.nl
 *
 *    Checker           : B. Romgens
 *    Affiliation       : Delft University of Technology
 *    E-mail address    : bart.romgens@gmail.com
 *
 *    Date created      : 12 August, 2010
 *    Last modified     : 2 February, 2010
 *
 *    References
 *
 *    Notes
 *
 *    Copyright (c) 2010-2011 Delft University of Technology.
 *
 *    This software is protected by national and international copyright.
 *    Any unauthorized use, reproduction or modification is unlawful and
 *    will be prosecuted. Commercial and non-private application of the
 *    software in any form is strictly prohibited unless otherwise granted
 *    by the authors.
 *
 *    The code is provided without any warranty; without even the implied
 *    warranty of merchantibility or fitness for a particular purpose.
 *
 *    Changelog
 *      YYMMDD    Author              Comment
 *      100914    K. Kumar            File created.
 *      100928    K. Kumar            Completed missing comments, changed
 *                                    writeIntegrationHistoryToFile( ) to
 *                                    writePropagationHistoryToFile( ).
 *      100929    B. Romgens          Spelling mistakes corrected.
 *      100929    K. Kumar            Added checked code written by D. Dirkx.
 *      110202    K. Kumar            Updated writePropagationHistoryToFile()
 *                                    to work with State*.
 */

#ifndef WRITINGOUTPUTTOFILE_H
#define WRITINGOUTPUTTOFILE_H

// Include statements.
#include <map>
#include "Astrodynamics/States/state.h"
#include "Mathematics/GeometricShapes/compositeSurfaceGeometry.h"
#include "Mathematics/GeometricShapes/singleSurfaceGeometry.h"

//! Writing output to file class.
/*!
 * Class containing all functionality pertaining to writing output to file.
 */
class WritingOutputToFile
{
public:

    //! Default destructor.
    /*!
     * Default destructor.
     */
    virtual ~WritingOutputToFile( ) { }

    //! Write propagation history to file.
    /*!
     * Writes propagation history stored in a map to file. The data is stored
     * in State objects.
     * \param propagationHistory Map of propagation history.
     * \param outputFilename Output filename.
     */
    void writePropagationHistoryToFile( std::map< double, State >& propagationHistory,
                                        const std::string& outputFilename );

    //! Write single surface geometry to a file.
    /*!
     * Writes the points on a SingleSurfaceGeometry object to a file.
     * Each row contains the x, y and z coordinate of the point, each next row
     * defines a single new point.
     * \param pointerToSingleSurfaceGeometry Geometry which is to be written.
     * \param numberOfLines Defines how many points are taken over the 1st
     *          independent variable.
     * \param numberOfPoints Defines how many points are taken over the 2nd
     *          independent variable.
     * \param filename Name of the file to which the points are written.
     * \param writeType Defines whether to append or write to the file given by
     *          filename, should be "a" for append and "w" for write.
     * \param isIndependentVariableInverted Boolean flag which if set to true
     *          inverts which independent variable is treated as 1st and which
     *          as 2nd.
     */
    void writeSingleSurfaceGeometryPointsToFile(
            SingleSurfaceGeometry* pointerToSingleSurfaceGeometry,
            const int& numberOfLines, const int& numberOfPoints,
            const std::string& filename, const int& writeType,
            const bool& isIndependentVariableInverted );

    //! Write composite surface geometry to a file.
    /*!
     *  Writes the single surface geometries in a composite surface geometry to
     *  a file. The writeSingleGeometryPointsToFile() function is called for
     *  each surface geometry.
     *  \param pointerToCompositeSurfaceGeometry Geometry from which there
     *          is to be written.
     *  \param arrayOfNumberOfLines Array of how many points to take over the 1st
     *          independent variables of single surface geometries.
     *  \param arrayOfNumberOfPoints Array of how many points to take over the 2nd
     *          independent variables of single surface geometries.
     *  \param filename Name of the file to which the points are written.
     *  \param writeType Defines whether to append or write to the file given
     *          by filename,  should be "a" for append and "w" for write.
     *  \param isIndependentVariableInvertedArray Array of booleans which if
     *          set to true invert which independent variable is treated as 1st
     *          and which as 2nd for each single surface geometry.
     */
    void writeCompositeSurfaceGeometryPointsToFile(
            CompositeSurfaceGeometry* pointerToCompositeSurfaceGeometry,
            int* arrayOfNumberOfLines, int* arrayOfNumberOfPoints,
            const std::string& filename, const int& writeType,
            bool* isIndependentVariableInvertedArray );

protected:

private:

    //! Output file.
    /*!
     * Output file.
     */
    static std::ofstream outputFile_;
};

#endif // WRITINGOUTPUTTOFILE_H

// End of file.
