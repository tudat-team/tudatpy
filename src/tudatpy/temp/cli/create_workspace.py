if __name__ == "__main__":
    import argparse

    parser = argparse.ArgumentParser(
        description=("Configure a workspace given "
                     "a tudat-space.yml file.")
    )
    parser.add_argument(
        "forge_file_directory",
        help=(
            "the directory containing the conda-forge.yml file "
            "used to configure the feedstock"
        ),
    )

    args = parser.parse_args()
    main(args.forge_file_directory)
