
def get_current_year():
    """
    gets current year
    :return:
    """
    return datetime.datetime.now().year

COPYRIGHT = (
f"""
/*    Copyright (c) 2010-{get_current_year()}, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */
"""
)

def get_all_files(path):
    """
    gets all files recursively
    :param path:
    :return:
    """
    for root, dirs, files in os.walk(path):
        for file in files:
            yield os.path.join(root, file)


def