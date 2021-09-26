from pydantic import BaseModel
from typing import List, Optional
from multidoc.parsing.io import yaml2dict


#
# class DocStyle(BaseModel):
#     pass
#
# #
# class NumpyDoc(DocStyle):
#     name: str
#     short_summary: Optional[str]
#     deprecation_warning: Optional[str]
#     extended_summary: Optional[str]
#     parameters: Optional[List[Parameter]]
#     returns: Optional[Returns] or Optional[List[Returns]]
#     yields: Optional[List[Yields] or Yields]
#     other_parameters: Optional[List[Parameter]]
#     raises: Optional[List[Raises] or Raises]
#     warns: Optional[List[Raises] or Raises]
#     warnings: Optional[str]
#     see_also: Optional[str]
#     notes: Optional[str]
#     references: Optional[str]
#     examples: Optional[str]


# class APIElement(BaseModel):
#     name: str

class APIDeclaration:
    pass


class FileBased(BaseModel):
    """FileBased declaration ``pydantic.BaseModel`` data structure.

    Attributes
    ----------


    """

    @classmethod
    def parse_yaml(cls, path, **kwargs):
        """

        Parameters
        ----------
        path
        local

        Returns
        -------

        """
        return cls.parse_obj(yaml2dict(path, **kwargs))


class DirBased(FileBased):

    @classmethod
    def parse_yaml(cls, path, local: dict = None):
        """

        Parameters
        ----------
        path
        local

        Returns
        -------

        """
        local = local if local is not None else {}
        return cls.parse_obj(yaml2dict(path, local))


class Parameter(BaseModel):
    """Parameter docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Parameter.name : str
    Parameter.type : Optional[str]
    Parameter.description : Optional[str]

    Examples
    --------


    .. code-block:: yaml
      :caption: Used in :class:`~multidoc.parsing.Function` parameters.

      functions:
       - name: foo_function
         parameters:
          - name: bar_parameter
            type: Any
            description: "The bar parameter"

    .. code-block:: yaml
       :caption: Used in :class:`~multidoc.parsing.Class` method parameters.

       classes:
        - name: FooClass
          methods:
           - name: bar_method
             parameters:
              - name: bar_parameter
                type: Any
                description: "The bar parameter"

    """
    name: str
    type: Optional[str]
    description: Optional[str]


class Returns(BaseModel):
    """Returns docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Returns.name : Optional[str]
    Returns.type : Optional[str]
    Returns.description : str

    """

    name: Optional[str]
    type: Optional[str]
    description: Optional[str]


class Yields(BaseModel):
    """Yields docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Yields.name: Optional[str]
    Yields.type: Optional[str]
    Yields.description: Optional[str]

    """
    name: Optional[str]
    type: Optional[str]
    description: Optional[str]


class Raises(BaseModel):
    """Raises docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Raises.name: str
    Raises.type: Optional[str]
    Raises.description: Optional[str]

    """
    name: str
    type: Optional[str]
    description: Optional[str]


class Function(BaseModel):
    """Function docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Function.name : str
    Function.short_summary : Optional[str]  # test
    Function.deprecation_warning : Optional[str]
    Function.extended_summary : Optional[str]
    Function.parameters : Optional[List[Parameter]]
    Function.returns : Optional[Returns] or Optional[List[Returns]]
    Function.yields : Optional[List[Yields] or Yields]
    Function.other_parameters : Optional[List[Parameter]]
    Function.raises : Optional[List[Raises] or Raises]
    Function.warns : Optional[List[Raises] or Raises]
    Function.warnings : Optional[str]
    Function.see_also : Optional[str]
    Function.notes : Optional[str]
    Function.references : Optional[str]
    Function.examples : Optional[str]

    """
    name: str
    short_summary: Optional[str]
    deprecation_warning: Optional[str]
    extended_summary: Optional[str]
    parameters: Optional[List[Parameter]]
    returns: Optional[Returns] or Optional[List[Returns]]
    yields: Optional[List[Yields] or Yields]
    other_parameters: Optional[List[Parameter]]
    raises: Optional[List[Raises] or Raises]
    warns: Optional[List[Raises] or Raises]
    warnings: Optional[str]
    see_also: Optional[str]
    notes: Optional[str]
    references: Optional[str]
    examples: Optional[str]


class Class(BaseModel):
    """Class docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Class.name : str
    Class.short_summary : Optional[str]  # test
    Class.deprecation_warning : Optional[str]
    Class.extended_summary : Optional[str]
    Class.parameters : Optional[List[Parameter]]
    Class.returns : Optional[Returns] or Optional[List[Returns]]
    Class.yields : Optional[List[Yields] or Yields]
    Class.other_parameters : Optional[List[Parameter]]
    Class.raises : Optional[List[Raises] or Raises]
    Class.warns : Optional[List[Raises] or Raises]
    Class.warnings : Optional[str]
    Class.see_also : Optional[str]
    Class.notes : Optional[str]
    Class.references : Optional[str]
    Class.examples : Optional[str]
    Class.methods : Optional[str]

    """

    name: str
    short_summary: Optional[str]
    deprecation_warning: Optional[str]
    extended_summary: Optional[str]
    parameters: Optional[List[Parameter]]
    attributes: Optional[List[Parameter]]
    yields: Optional[List[Yields] or Yields]
    other_parameters: Optional[List[Parameter] or Parameter]
    raises: Optional[List[Raises] or Raises]
    warns: Optional[List[Raises] or Raises]
    warnings: Optional[str]
    see_also: Optional[str]
    notes: Optional[str]
    references: Optional[str]
    examples: Optional[str]
    methods: Optional[List[Function]]


class Constant(BaseModel):
    """Constant docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Constant.summary: str
    Constant.extended_summary: Optional[str]
    Constant.see_also: Optional[str]
    Constant.references: Optional[str]
    Constant.examples: Optional[str]

    """
    summary: str
    extended_summary: Optional[str]
    see_also: Optional[str]
    references: Optional[str]
    examples: Optional[str]


class Config(BaseModel):
    """Multidoc module config ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Config.name: Optional[str]
    Config.version: Optional[str]

    """
    name: Optional[str]
    version: Optional[str]


class Module(FileBased):
    """Module docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Module.config: Optional[Config]
    Module.summary: Optional[str]
    Module.extended_summary: Optional[str]
    Module.routine_listings: Optional[str]
    Module.see_also: Optional[str]
    Module.notes: Optional[str]
    Module.references: Optional[str]
    Module.examples: Optional[str]
    Module.classes: Optional[List[Class]]
    Module.functions: Optional[List[Function]]
    Module.constants: Optional[List[Constant]]

    """
    config: Optional[Config]

    summary: Optional[str]
    extended_summary: Optional[str]
    routine_listings: Optional[str]
    see_also: Optional[str]
    notes: Optional[str]
    references: Optional[str]
    examples: Optional[str]

    # MODULE structure
    classes: Optional[List[Class]]
    functions: Optional[List[Function]]
    constants: Optional[List[Constant]]


class Package(Module):
    """Module docstring ``pydantic.BaseModel`` data structure.

    Attributes
    ----------
    Package.config: Optional[Config]
    Package.summary: Optional[str]
    Package.extended_summary: Optional[str]
    Package.routine_listings: Optional[str]
    Package.see_also: Optional[str]
    Package.notes: Optional[str]
    Package.references: Optional[str]
    Package.examples: Optional[str]
    Package.classes: Optional[List[Class]]
    Package.functions: Optional[List[Function]]
    Package.constants: Optional[List[Constant]]
    Package.modules: Optional[List[str]]

    """
    modules: Optional[List[str]]

    # def __init__(self):
    #     super().__init__()
    #


if __name__ == "__main__":
    import json

    # r = Package.schema()
    # print(json.dumps(r, indent=4))
    # p = Package()
    # print(p)
