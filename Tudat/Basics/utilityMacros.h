/*    Copyright (c) 2010-2015, Delft University of Technology
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
 *      120202    S. Billemont      File created.
 *      121212    K. Kumar          Added header include guard.
 *
 *    References
 *
 *    Notes
 *
 */

#ifndef TUDAT_UTILITY_MACROS_H
#define TUDAT_UTILITY_MACROS_H

//! Suppress compile time unused parameter warnings.
/*!
 * This macro suppresses the warning that a parameter passed to a function/method is not used, but
 * still passed in. This commonly occurs when overriding virtual methods.
 * 
 * Example:
 * int foo( int bar) {
 *     TUDAT_UNUSED_PARAMETER(bar);
 *     return -1;
 * }
 *
 * \param unusedParameter Parameter that is not being used in a function.
 */
#define TUDAT_UNUSED_PARAMETER( unusedParameter ) { ( void ) unusedParameter; }

//! Mark a declaration as deprecated.
/*!
 * When a declaration is marked as deprecated, it should be avoided, as most likely it performs 
 * something unwanted, superseded by another functionality or marked for removal in a future a 
 * release. The user should check where the deprecation warning comes from in order to get
 * information from the developer on how to resolve the warning. Commonly one should use another,
 * recommended function or alternative code block.
 *
 * Snippets:
 * TUDAT_DEPRECATED("Deprecated function!", void f() );
 * TUDAT_DEPRECATED("Deprecated type!",     typedef int my_int);
 * TUDAT_DEPRECATED("Deprecated variable!", int bar);
 *
 * Example:
 * //! Do something. 
 * // Deprecated; use Foo::NewFunc(int, float) instead.
 * // \see Foo::NewFunc(int, float)
 * TUDAT_DEPRECATED("Deprecated, you should use Foo::NewFunc(int, float) instead.",
 *                   void OldFunc(int a, float b));
 *
 * //! Do something better
 * void NewFunc(int a, double b);
 *
 * \param message Warning message to show when compiling. One should mention "Deprecated", and why
 *                  the function was deprecated in a short one liner.
 * \param expression The expression which has to be deprecated.
 */
#ifdef __GNUC__
#define TUDAT_DEPRECATED( message, expression ) expression __attribute__ ( ( deprecated( message ) ) )
#elif defined( _MSC_VER )
#define TUDAT_DEPRECATED( message, expression ) __declspec( deprecated( message ) ) expression
#else
#pragma message( "WARNING: You need to implement DEPRECATED for this compiler" )
#define TUDAT_DEPRECATED( message, expression ) expression
#endif

#endif // TUDAT_UTILITY_MACROS_H
