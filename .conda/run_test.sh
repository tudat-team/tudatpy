#!/bin/bash

python3 -c "import tudatpy.constants as c; print(c.SPEED_OF_LIGHT)"
python3 -c "import tudatpy.spice_interface as spice; spice.load_standard_spice_kernels()"
python3 -c "import tudatpy.io as io; print(io.get_tudat_data_path())"