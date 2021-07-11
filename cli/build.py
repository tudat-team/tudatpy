import subprocess
import os
import sys

mp = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "tudat-multidoc", "multidoc")
sys.path.insert(0, mp)

from multidoc.config import process_config


def main(projects, config):
    projects_found, common_configs, env = process_config(config)
    bash_env = os.environ.copy()
    bash_env.update(env)
    build_scripts = []
    # p = os.path.dirname(os.path.realpath(os.path.join(__file__, "..", "..")))
    projects = [projects] if type(projects) == str else projects
    for project in projects:
        print(project)

        # superimpose project config onto common config.
        project_kwargs = common_configs.copy()
        project_kwargs.update(config[project])
        if "build_script" in project_kwargs.keys():
            build_scripts.append(
                os.path.join(os.path.dirname(os.path.abspath(project_kwargs["configpath"])),
                             project_kwargs["build_script"]))
    for bs in set(build_scripts):
        subprocess.call(bs, env=bash_env, cwd=os.path.dirname(os.path.abspath(common_configs["configpath"])))
