import os

def parts(path):
    p, f = os.path.split(path)
    return parts(p) + [f] if f else [p]
