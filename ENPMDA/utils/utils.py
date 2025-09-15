import re
import os

def atoi(text):
    return int(text) if text.isdigit() else text


def natural_keys(text):
    """
    alist.sort(key=natural_keys) sorts in human order
    http://nedbatchelder.com/blog/200712/human_sorting.html
    """
    return [atoi(c) for c in re.split(r"(\d+)", text)]

def normalize_user_path(p: str) -> str:
    """Treat './<abs>' like '<abs>'; otherwise return abspath(p)."""
    dot_slash = '.' + os.sep
    if p.startswith(dot_slash) and os.path.isabs(p[len(dot_slash):]):
        p = p[len(dot_slash):]
    return os.path.abspath(p)