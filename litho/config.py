"""
loads a configuration from 3 files, high priority overwrites low priority:

1. A config.yml found in the current working directory (high priority)
2. ~/.dimmilitho.yml specific for the machine
3. the default config is in this file

"""

__all__ = ["CONFIG"]

import logging
import pathlib

import hiyapyco

default_config = """
keySample: valueSample
"""

home = pathlib.Path.home()
cwd = pathlib.Path.cwd()
cwd_config = cwd / "config.yml"
home_config = home / ".dimmilitho.yml"
module_path = pathlib.Path(__file__).parent.absolute()
repo_path = module_path.parent

CONFIG = hiyapyco.load(
    str(default_config),
    str(cwd_config),
    str(home_config),
    failonmissingfiles=False,
    loglevelmissingfiles=logging.DEBUG,
)

CONFIG["module_path"] = module_path
CONFIG["repo_path"] = repo_path
CONFIG["gdslib"] = repo_path / "gdslib"
CONFIG["samples"] = module_path / "samples"

CONFIG["scalerate"] = 45 / 70

if __name__ == "__main__":
    print(CONFIG["repo"])
