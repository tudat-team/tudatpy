import subprocess
import os
import sys

mp = os.path.join(os.path.abspath(os.path.dirname(__file__)), "..", "tudat-multidoc", "multidoc")
sys.path.insert(0, mp)

from multidoc.config import process_config


def main(projects, config):
    projects_found, common_configs, env = process_config(config)
    cwd = os.path.dirname(os.path.abspath(common_configs["configpath"]))
    bash_env = os.environ.copy()
    bash_env.update(env)
    projects = [projects] if type(projects) == str else projects

    for project in [projects[1]]:# tudatpy only!

        # superimpose project config onto common config.
        project_kwargs = common_configs.copy()
        project_kwargs.update(config[project])
        output_sphinx = project_kwargs.pop("sphinx_output_directory", f".docs-output/{project}/sphinx")
        documented_output = project_kwargs.pop("documented_output", f".{project}-documented")

        # if cpp, build doxygen first
        if project_kwargs["project_type"] == "cpp":
            project_doxygen_source = project_kwargs.pop("project_doxygen_source", "docs/doxygen")
            doxygen_conf_in = os.path.join(cwd, documented_output, project_doxygen_source, "Doxyfile.in")
            doxygen_conf_tmp = os.path.join(cwd, documented_output, project_doxygen_source, "Doxyfile")
            with open(doxygen_conf_in, "r") as f:
                s = f.read()
                relative_doxygen_output_directory = project_kwargs.pop("project_kwargs", f".docs-output/{project}/doxygen")
                doxygen_output_directory = os.path.abspath(os.path.join(cwd, relative_doxygen_output_directory))
                if not os.path.exists(doxygen_output_directory):
                    os.makedirs(doxygen_output_directory)
                s = s.replace("@DOXYGEN_OUTPUT_DIRECTORY", doxygen_output_directory)

            with open(doxygen_conf_tmp, "w") as f:
                f.write(s)

            bash_env.update({"DOXYGEN_OUTPUT_DIRECTORY": doxygen_output_directory})

            command_doxygen = "doxygen Doxyfile"
            # TODO: This is wrong, musn't be executed on source of dox.
            process = subprocess.Popen(command_doxygen.split(), env=bash_env, cwd=os.path.join(cwd, documented_output, project_doxygen_source))
            stdout, stderr = process.communicate()
            # os.remove(doxygen_conf_tmp)

        # build sphinx docs
        if project_kwargs["project_type"] == "cpp":
            pass
            #    project_sphinx_source = os.path.join(cwd, documented_output, project_kwargs.pop("project_sphinx_source", "docs/sphinx/source"))
        elif project_kwargs["project_type"] == "py":
            build_directory = project_kwargs.pop("build_directory", "cmake-build-release")
            project_sphinx_source = os.path.join(build_directory, documented_output, project_kwargs.pop("project_sphinx_source", "docs/sphinx/source"))
            command = f"sphinx-build -b html {project_sphinx_source} {output_sphinx}"
            process = subprocess.Popen(command.split(), env=bash_env, cwd=cwd)
            stdout, stderr = process.communicate()
        else:
            raise KeyError("Project type can only be py or cpp")

        if not os.path.exists(output_sphinx):
            os.makedirs(output_sphinx)
        #command = f"sphinx-build -b html {project_sphinx_source} {output_sphinx}"
#
        #process = subprocess.Popen(command.split(), env=bash_env, cwd=cwd)
        #stdout, stderr = process.communicate()
