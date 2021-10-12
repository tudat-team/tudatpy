import os

from . import logger
import yaml
from yaml.parser import ParserError, ScannerError
from multidoc.regex import p_api_tag


def yaml2dict(path, include_name_error=False, **kwargs):
    """Yaml file parser.

    Parameters
    ----------
    path : os.PathLike or str
        Path of the ``yaml`` file to be loaded.
    _locals : dict[str, Any]
        List of definitions to parse in the yaml files. See examples for how
        they affect yaml loading.
    include_name_error : bool, default=True
        Include tag evaluations that return a NameError

    Examples
    ========
    Given the following example ``yaml`` file:

    .. code-block:: yaml
        :caption: example.yaml

        package:
          name: name-cpp    # [cpp]
          name: name-py     # [py]

        modules:
          - module
          - module-py       # [py]
          - module-cpp      # [cpp]
          - module-not-cpp  # [not cpp]
          - module-not-py   # [not py]
          - module-both     # [py or cpp]
          - module-both     # [py and cpp]

    >>> yaml2dict("../tests/example.yaml")
    {'package': None, 'modules': ['module']}

    >>> yaml2dict("../tests/example.yaml", cpp=True)
    {'package': {'name': 'name-cpp'}, 'modules': ['module', 'module-cpp']}

    >>> yaml2dict("../tests/example.yaml", py=True)
    {'package': {'name': 'name-py'}, 'modules': ['module', 'module-py']}

    >>> yaml2dict("../tests/example.yaml", cpp=True, py=False)
    {'package': {'name': 'name-cpp'}, 'modules': ['module', 'module-cpp', 'module-not-py']}


    Returns
    -------
    dict
        ``yaml`` loaded into memory as Python dict, parsed for definitions.

    """
    # update the local variable space for eval of tag
    locals().update(kwargs) if kwargs else None
    logger.info(f"Parsing yaml file: {path} with kwargs: {kwargs}")
    # open the yaml file!
    with open(path) as file:
        # logger.info(f"Parsing yaml file :{file}")

        # read the raw lines
        raw_lines = file.readlines()

        # create list for processed lines
        processed_lines = []

        # iterate through lines
        for line in raw_lines:
            match = p_api_tag.match(line)
            if match:  # there's an expr on this line
                try:
                    # if the expr is True, add to lines, else ignore
                    if eval(match.group("expr")):
                        processed_lines += line
                except NameError:
                    # if the expr raises a NameError (e.g. undefined var in expr)
                    if include_name_error:
                        processed_lines += line
                    else:
                        pass
            else:
                processed_lines += line
        # apply line check truth to all raw lines, retrieving those that return true
        try:
            return yaml.load("".join(processed_lines), yaml.Loader)
        except (ParserError, ScannerError) as e:
            broken_path = f"BROKEN-{kwargs}-{os.path.basename(path)}"
            with open(broken_path, "w") as f:
                f.write("".join(processed_lines))
            logger.error(f"Broken .yaml file {path} dumped as {broken_path}.")
            raise e


if __name__ == "__main__":
    d = yaml2dict("../testing/example-1.yaml", cpp=True, py=False)
    print(d)
