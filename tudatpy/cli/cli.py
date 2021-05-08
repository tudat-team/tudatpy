import argparse
import os
import logging
import sys
# from .._version import __version__
import create_project


class Subcommand(object):
    #: The name of the subcommand
    subcommand = None
    aliases = []

    def __init__(self, parser, help=None):
        subcommand_parser = parser.add_parser(
            self.subcommand, help=help, description=help, aliases=self.aliases
        )
        subcommand_parser.set_defaults(subcommand_func=self)
        self.subcommand_parser = subcommand_parser

    def __call__(self, args):
        pass


class CreateProject(Subcommand):
    subcommand = "project"
    aliases = ["create"]

    def __init__(self, parser):
        super(CreateProject, self).__init__(
            parser,
            "Create a tudat project.",
        )
        scp = self.subcommand_parser
        scp.add_argument(
            "-n",
            "--name",
            help="The name of the project to be created.",
        )
        scp.add_argument(
            "-d",
            "--directory",
            default=os.getcwd(),
            help="The directory of the project to be created in. (default = .)",
        )
        scp.add_argument(
            "-t",
            "--type",
            default="hybrid",
            help="The type of project to be created. (default = hybrid)",
        )
        scp.add_argument(
            "-c",
            "--config",
            default=None,
            help="The config of the project to be created. (default = None)",
        )

    def __call__(self, args):
        create_project.main(
            args.name,
            args.type,
            args.directory,
            args.config
        )


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        prog="tudat",
        description="a tool for performing research and education in astrodynamics.",
    )
    subparser = parser.add_subparsers()
    # TODO: Consider allowing plugins/extensions using entry_points.
    # https://reinout.vanrees.org/weblog/2010/01/06/zest-releaser-entry-points.html
    for subcommand in Subcommand.__subclasses__():
        subcommand(subparser)
    # And the alias for rerender
    parser.add_argument(
        "--version",
        action="version",
        # version=__version__,
        help="Show conda-smithy's version, and exit.",
    )

    if not sys.argv[1:]:
        args = parser.parse_args(["--help"])
    else:
        args = parser.parse_args()

    args.subcommand_func(args)


if __name__ == "__main__":
    main()
