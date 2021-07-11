import argparse
import os
import logging
import sys
# from .._version import __version__
import document
import build
import sphinx
import configparser


def load_config(path):
    config = configparser.ConfigParser()
    config.read(path)
    config["COMMON"]["configpath"] = path
    return config


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


class Document(Subcommand):
    subcommand = "document"
    aliases = ["d", "doc"]

    def __init__(self, parser):
        super(Document, self).__init__(
            parser,
            "tudat-bundle command line interface.",
        )
        scp = self.subcommand_parser
        scp.add_argument(
            "-p",
            "--project",
            # action="append",
            default=["tudat", "tudatpy"],
            help="Generate documented source, build projects and build documentation.",
        )
        scp.add_argument(
            "-c",
            "--config",
            default=".multidoc.cfg",
            help="Generate documented source, build projects and build documentation.",
        )

    def __call__(self, args):
        if not os.path.isfile(args.config):
            raise FileNotFoundError(".multidoc.cfg not found in current"
                                    " directory, or was not provided correctly as "
                                    "an argument with -c/--config")
        config = load_config(args.config)
        document.main(
            args.project,
            config
        )


class Build(Subcommand):
    subcommand = "build"
    aliases = ["b"]

    def __init__(self, parser):
        super(Build, self).__init__(
            parser,
            "tudat-bundle build command.",
        )
        scp = self.subcommand_parser
        scp.add_argument(
            "-p",
            "--project",
            # action="append",
            default=["tudat", "tudatpy"],
            help="Generate documented source, build projects and build documentation.",
        )
        scp.add_argument(
            "-c",
            "--config",
            default=".multidoc.cfg",
            help="Generate documented source, build projects and build documentation.",
        )

    def __call__(self, args):
        if not os.path.isfile(args.config):
            raise FileNotFoundError(".multidoc.cfg not found in current"
                                    " directory, or was not provided correctly as "
                                    "an argument with -c/--config")
        config = load_config(args.config)

        build.main(
            args.project,
            config
        )


class Sphinx(Subcommand):
    subcommand = "sphinx"
    aliases = ["s"]

    def __init__(self, parser):
        super(Sphinx, self).__init__(
            parser,
            "tudat-bundle sphinx command.",
        )
        scp = self.subcommand_parser
        scp.add_argument(
            "-p",
            "--project",
            # action="append",
            default=["tudat", "tudatpy"],
            help="Select project.",
        )
        scp.add_argument(
            "-c",
            "--config",
            default=".multidoc.cfg",
            help="Generate documented source, build projects and build documentation.",
        )

    def __call__(self, args):
        if not os.path.isfile(args.config):
            raise FileNotFoundError(".multidoc.cfg not found in current"
                                    " directory, or was not provided correctly as "
                                    "an argument with -c/--config")
        config = load_config(args.config)

        sphinx.main(
            args.project,
            config
        )
class All(Subcommand):
    subcommand = "all"
    aliases = ["a"]

    def __init__(self, parser):
        super(All, self).__init__(
            parser,
            "tudat-bundle all command.",
        )
        scp = self.subcommand_parser
        scp.add_argument(
            "-p",
            "--project",
            # action="append",
            default=["tudat", "tudatpy"],
            help="Select project.",
        )
        scp.add_argument(
            "-c",
            "--config",
            default=".multidoc.cfg",
            help="Provide .multidoc.cfg path.",
        )

    def __call__(self, args):
        if not os.path.isfile(args.config):
            raise FileNotFoundError(".multidoc.cfg not found in current"
                                    " directory, or was not provided correctly as "
                                    "an argument with -c/--config")
        config = load_config(args.config)
        document.main(
            args.project,
            config
        )
        config = load_config(args.config)
        build.main(
            args.project,
            config
        )
        config = load_config(args.config)
        sphinx.main(
            args.project,
            config
        )

# class BuildBundle(Subcommand):
#     subcommand = "build"
#     aliases = []
#
#     def __init__(self, parser):
#         super(BuildBundle, self).__init__(
#             parser,
#             "tudat-bundle build command.",
#         )
#         scp = self.subcommand_parser
#         scp.add_argument(
#             "-t",
#             "--tests",
#             default=False,
#             help="Build tudat-bundle.",
#         )
#
#     def __call__(self, args):
#         build_bundle.main(
#             args.tests
#         )


def main():
    logging.basicConfig(level=logging.INFO)

    parser = argparse.ArgumentParser(
        prog="multidoc",
        description="a tool for working with the tudat-bundle.",
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
        help="Show version, and exit.",
    )

    if not sys.argv[1:]:
        args = parser.parse_args(["--help"])
    else:
        args = parser.parse_args()

    args.subcommand_func(args)


if __name__ == "__main__":
    main()
