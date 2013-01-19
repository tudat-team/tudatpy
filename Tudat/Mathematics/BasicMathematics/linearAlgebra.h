/*    Copyright (c) 2010-2013, Delft University of Technology
 *    All rights reserved.
 *
 *    Redistribution and use in source and binary forms, with or without modification, are
 *    permitted provided that the following conditions are met:
 *      - Redistributions of source code must retain the above copyright notice, this list of
 *        conditions and the following disclaimer.
 *      - Redistributions in binary form must reproduce the above copyright notice, this list of
 *        conditions and the following disclaimer in the documentation and/or other materials
 *        provided with the distribution.
 *      - Neither the name of the Delft University of Technology nor the names of its contributors
 *        may be used to endorse or promote products derived from this software without specific
 *        prior written permission.
 *
 *    THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS
 *    OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF
 *    MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE
 *    COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 *    EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE
 *    GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED
 *    AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 *    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED
 *    OF THE POSSIBILITY OF SUCH DAMAGE.
 *
 *    Changelog
 *      YYMMDD    Author            Comment
 *      130116    K. Kumar          File created.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_LINEAR_ALGEBRA_H
#define TUDAT_LINEAR_ALGEBRA_H

#include <Eigen/Core>

namespace tudat
{
namespace basic_mathematics
{
namespace linear_algebra
{

//! Flip matrix rows.
/*!
 * Flips all rows of an Eigen-matrix, i.e., order of rows is reversed.
 * \param matrixToFlip Matrix that should be flipped. The flipping is done in place, i.e. the input
 *          matrix is modified directly (no copies are made).
 */
static inline void flipMatrixRows( Eigen::MatrixXd& matrixToFlip )
{
    // Loop through each row.
    for ( int i = 0; i < static_cast< int >( matrixToFlip.rows( ) / 2 ); i++ )
    {
        // Save top row in a temporary matrix.
        const Eigen::MatrixXd temporaryRow = matrixToFlip.row( i );

        // Set top row to bottom row.
        matrixToFlip.row( i ) = matrixToFlip.row( matrixToFlip.rows( ) - 1 - i );

        // Set bottom row to temporarily stored matrix.
        matrixToFlip.row( matrixToFlip.rows( ) - 1 - i ) = temporaryRow;
    }
}

} // linear_algebra
} // basic_mathematics
} // tudat

#endif // TUDAT_LINEAR_ALGEBRA_H
