//
// Created by ggarrett on 19-05-20.
//

#ifndef TUDAT_UTIL_H
#define TUDAT_UTIL_H


namespace tudat {
namespace input_output {
//! Function to compare if two lists of aerodynamic coefficient independent variables are equal
/*!
 * Function to compare if two lists of aerodynamic coefficient independent variables (vector of vector of doubles) are equal
 * \param list1 First list that is to be compared.
 * \param list2 Second list that is to be compared.
 * \return True of the two lists are completely equal in size and contents, false otherwise.
 */
bool compareIndependentVariables( const std::vector< std::vector< double > >& list1,
                                  const std::vector< std::vector< double > >& list2 );

}
}
#endif //TUDAT_UTIL_H

