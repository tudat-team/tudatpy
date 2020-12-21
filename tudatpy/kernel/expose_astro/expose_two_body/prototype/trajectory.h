/*    Copyright (c) 2010-2018, Delft University of Technology
 *    All rights reserved
 *
 *    This file is part of the Tudat. Redistribution and use in source and
 *    binary forms, with or without modification, are permitted exclusively
 *    under the terms of the Modified BSD license. You should have received
 *    a copy of the license with this file. If not, please or visit:
 *    http://tudat.tudelft.nl/LICENSE.
 */

#ifndef TUDATPY_TRAJECTORY_H
#define TUDATPY_TRAJECTORY_H

class Trajectory {

  Trajectory(
      std::shared_ptr<BaseState> referenceState,
      std::shared_ptr<ReferenceFrame> referenceFrame,
      double referenceEpoch = 0
      );

  sample(int n_points)

};

#endif//TUDATPY_TRAJECTORY_H
