''' 
Copyright (c) 2010-2020, Delft University of Technology
All rigths reserved

This file is part of the Tudat. Redistribution and use in source and 
binary forms, with or without modification, are permitted exclusively
under the terms of the Modified BSD license. You should have received
a copy of the license with this file. If not, please or visit:
http://tudat.tudelft.nl/LICENSE.
'''

import time
import datetime


def elapsed_time_since_J2000_to_date(elapsed_time: float) -> str:

    # Transform the reference date to UNIX time
    reference_UNIX_time = time.mktime(datetime.datetime.strptime("2020-01-01", "%Y-%m-%d").timetuple())

    # Add elapsed time to reference date
    UNIX_time = int(reference_UNIX_time + elapsed_time)

    # Transform UNIX time to date
    date = datetime.datetime.utcfromtimestamp(UNIX_time).strftime('%Y-%m-%d')

    return date
