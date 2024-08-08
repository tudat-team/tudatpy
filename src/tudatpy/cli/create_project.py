import os
import jinja2
import shutil
import yaml

CMAKELISTS_FILENAME = "CMakeLists.txt"
CONFIG_FILENAME = "tudat-space.yml"
MAIN_PY_FILENAME = "main.py"
KERNEL_FILENAME = "pybind_module.cpp"
MAIN_CPP_FILENAME = "main.cpp"


def yaml2dict(path):
    with open(path) as file:
        dict_ = yaml.load(file, Loader=yaml.FullLoader)
        return dict_


def render_templates(templates, path, **template_kwargs):
    for template in templates:
        file_path = template.filename.replace(".template", "")
        file_name = os.path.split(file_path)[-1]
        with open(os.path.join(path, file_name), "w") as f:
            f.write(template.render(template_kwargs))


# def render_simulation_template(simulation)


def create_project(project_name,
                   project_type,
                   project_path=".",
                   template_path="./templates",
                   **kwargs):
    """

    Parameters
    ----------
    project_name: str
    project_type: str {"hybrid", "cpp", "python"}
    project_git: str [optional]

    Returns
    -------

    """

    support_folder = kwargs.get("support_folder", "support")
    kernel_folder = kwargs.get("kernel_folder", "kernel")

    template_loader = jinja2.FileSystemLoader(searchpath=template_path)
    template_env = jinja2.Environment(loader=template_loader)

    abs_path_project = os.path.join(project_path, project_name)
    abs_path_kernel = os.path.join(abs_path_project, kernel_folder)
    abs_support_path = os.path.join(abs_path_project, support_folder)

    os.makedirs(abs_path_project)
    os.makedirs(abs_support_path)

    kernel_templates = []
    support_templates = []
    root_templates = []

    config_template = template_env.get_template(CONFIG_FILENAME)
    root_templates.append(config_template)

    if (project_type == "hybrid") | (project_type == "python"):
        main_py_template = template_env.get_template(MAIN_PY_FILENAME)
        root_templates.append(main_py_template)

    if (project_type == "hybrid") | (project_type == "cpp"):
        cmakelists_template = template_env.get_template(CMAKELISTS_FILENAME)
        root_templates.append(cmakelists_template)

        if project_type == "hybrid":
            os.makedirs(abs_path_kernel)
            kernel_cpp_template = template_env.get_template(KERNEL_FILENAME)
            kernel_templates.append(kernel_cpp_template)

        elif project_type == "cpp":
            main_cpp_template = template_env.get_template(MAIN_CPP_FILENAME)
            root_templates.append(main_cpp_template)

    elif project_type == "python":
        pass
        # main_py_template = template_env.get_template("main.py.template")
        # root_templates.append(main_py_template)

    template_kwargs = {"project_name": project_name,
                       "project_type": project_type,
                       **kwargs}

    render_templates(root_templates, abs_path_project, **template_kwargs)
    render_templates(kernel_templates, abs_path_kernel, **template_kwargs)
    render_templates(support_templates, abs_support_path, **template_kwargs)


def test():
    try:
        create_project(
            project_name="test_hybrid",
            project_type="hybrid",
            project_path="./test_space",
            **yaml2dict("config/base.yml")
        )
        create_project(
            project_name="test_cpp",
            project_type="cpp",
            project_path="./test_space",
            **yaml2dict("config/base.yml")
        )
        create_project(
            project_name="test_python",
            project_type="python",
            project_path="./test_space",
            **yaml2dict("config/base.yml")
        )
    except FileExistsError:
        shutil.rmtree("./test_space")
    finally:
        pass


def main(
        project_name,
        project_type="python",
        project_path=".",
        project_config_directory=None
):
    kwargs = yaml2dict("config/base.yml")
    if project_config_directory:
        kwargs.update(yaml2dict(project_config_directory))
    create_project(
        project_name=project_name,
        project_type=project_type,
        project_path=project_path,
        **kwargs
    )


if __name__ == "__main__":
    test()
    # import argparse
    #
    # parser = argparse.ArgumentParser(
    #     description=("Configure a project given a conda-forge.yml file.")
    # )
    # parser.add_argument(
    #     "project_name",
    #     help=(
    #         "the name of the project that is to be created"
    #     ),
    # )
    # args = parser.parse_args()
    # main(args.project_name)
