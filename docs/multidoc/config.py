def process_config(config):
    config = config._sections
    projects_found = list(config.keys())
    projects_found.remove("COMMON") if "COMMON" in config.keys() else False
    projects_found.remove("ENV") if "ENV" in config.keys() else False
    common_configs = config.pop("COMMON", {})
    environment = config.pop("ENV", {})
    return projects_found, common_configs, environment
