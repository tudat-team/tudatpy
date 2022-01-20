import jinja2

template_loader = jinja2.FileSystemLoader(searchpath=template_path)
template_env = jinja2.Environment(loader=template_loader)

