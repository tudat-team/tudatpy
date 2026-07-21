"""Validate the packaged acceleration contract against Tudatpy."""

from tudatpy.json_interface import ACCELERATION_CONTRACT_PATH
from tudatpy.json_interface.validate_contract import main

if __name__ == "__main__":
    main(
        [
            str(ACCELERATION_CONTRACT_PATH),
            "tudatpy.dynamics.propagation_setup.acceleration",
            "--types",
            "tudatpy.dynamics.environment_setup",
            "tudatpy.dynamics.propagation_setup",
        ]
    )
