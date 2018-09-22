/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rigths reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 *
 */

#ifndef TUDAT_LINEAR_FIELD_TRANSFORM_H
#define TUDAT_LINEAR_FIELD_TRANSFORM_H

#include "Tudat/InputOutput/fieldTransform.h"

namespace tudat
{
namespace input_output
{

//! Linear field transform class.
/*!
 * This class can be used to linearly transform an input field (string). The linear transformation
 * is of the form: result = slope * input + intercept.
 * This class is derived from the FieldTransform abstract base class.
 * \sa FieldTransform
 */
class LinearFieldTransform : public FieldTransform
{
public:

    //! Constructor of the linear field transform.
    /*!
     * Constructor of the linear field transformation. Linear transformation is of the form:
     * \f$y=a*x+b\f$, where \f$a\f$ is the slope and b is the intercept.
     * \param aSlope Slope of the linear field transform.
     * \param anIntercept Intercept of the linear field transform.
     */
    LinearFieldTransform( const double aSlope, const double anIntercept )
        : slope( aSlope ), intercept( anIntercept )
    { }

    //! Default destructor.
    ~LinearFieldTransform( ) { }

    //! Transform input string.
    /*!
     * Returns a transformed string, according to the linear transformation:
     * result = slope * input + intercept.
     * \param input Input string.
     * \return Shared-pointer to transformed string.
     */
    std::shared_ptr< std::string > transform( const std::string& input );

protected:

    //! Slope of the transformation.
    const double slope;

    //! Intercept of the transformation.
    const double intercept;

private:
};

//! Typedef for shared-pointer to LinearFieldTransform object.
typedef std::shared_ptr< LinearFieldTransform > LinearFieldTransformPointer;

} // namespace input_output
} // namespace tudat

#endif // TUDAT_LINEAR_FIELD_TRANSFORM_H
