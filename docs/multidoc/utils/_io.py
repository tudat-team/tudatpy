import yaml


def yaml2dict(path):
    with open(path) as file:
        dict_ = yaml.load(file, Loader=yaml.FullLoader)
        return dict_
