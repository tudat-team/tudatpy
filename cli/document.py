import sys
import os

mp = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "tudat-multidoc", "multidoc")
sys.path.insert(0, mp)
import pathlib

from multidoc.generate import generate_pybind_documented, generate_cpp_documented
from multidoc.config import process_config
from multidoc.parsing import guess_project_type


class ProjectSourcMissing(Exception):
    pass


class ProjectDocstringsMissing(Exception):
    pass


def main(projects, config):
    projects_found, common_configs, _ = process_config(config)
    projects = [projects] if type(projects) == str else projects


    # for project in [projects[1]]: # tudatpy only!
    for project in projects: # all projects
        # confirm project exists in .multidoc.cfg
        if project not in projects_found:
            raise KeyError(f"{project} key not found in config file.")

        # superimpose project config onto common config.
        project_kwargs = common_configs.copy()
        project_kwargs.update(config[project])

        # determine then confirm determined project_src exists.
        project_src = project_kwargs.pop("project_src", f"./{project}")
        if not os.path.isdir(project_src):
            raise ProjectSourcMissing(
                f"project_api={project_src} was not found. Please confirm that the "
                f".multidoc.cfg file points to the correct source directory "
                f"for {project} (default: project_api='./{project}'")

        # determine then confirm determined project_api exists.
        project_api = project_kwargs.pop("project_api", f"./docstrings")
        if not os.path.isdir(project_api):
            raise ProjectDocstringsMissing(
                f"{project_api} was not found. Please confirm that the "
                f".multidoc.cfg file points to the correct api directory "
                f"for {project} (default='./docstrings'")

        # determine otherwise guess project type.

        project_type = project_kwargs.pop("project_type", None)
        if project_type is None:
            project_type = guess_project_type(f'./{project}')
        documented_output = project_kwargs.pop("documented_output", None)
        config_path = project_kwargs.pop("configpath", None)
        dest = os.path.join(os.path.dirname(os.path.abspath(config_path)), documented_output)

        if project_type == "py":
            generate_pybind_documented(
                api_prefix=project_api,
                target_src=project_src,
                dest=dest)
        #elif project_type == "cpp":
        #    generate_cpp_documented(
        #        api_prefix=project_api,
        #        target_src=project_src,
        #        dest=dest)
