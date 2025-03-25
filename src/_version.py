"""
New release versions are made using [Semantic Versioning](https://semver.org/).

> Given a version number MAJOR.MINOR.PATCH, increment the:

    1. MAJOR version when you make incompatible API changes,
    2. MINOR version when you add functionality in a backwards compatible manner, and
    3. PATCH version when you make backwards compatible bug fixes.

    Additional labels for pre-release and build metadata are available as extensions to the MAJOR.MINOR.PATCH format.

Examples:
"0.0.0.dev0"
"1.2.34.dev0"
"1.2.34a0"
"1.2.34"
"""

import os
import sys

__version__ = "0.0.0"  # Default value


def _read_version_from_pyproject():
    try:
        if sys.version_info >= (3, 11):
            import tomllib
        else:
            import tomli as tomllib
    except ImportError:
        return "0.0.0"

    parent_path = os.path.dirname(os.path.dirname(os.path.realpath(__file__)))
    pyproject_path = os.path.join(parent_path, "pyproject.toml")

    if not os.path.exists(pyproject_path):
        return "0.0.0"

    with open(pyproject_path, "rb") as f:
        data = tomllib.load(f)
        return data.get("project", {}).get("version", "0.0.0")


try:
    from importlib.metadata import PackageNotFoundError, version
except ImportError:
    # Python < 3.8, but should not be the case.
    version = None
    PackageNotFoundError = Exception

try:
    __version__ = version("pyhim")
except (PackageNotFoundError, TypeError):
    __version__ = _read_version_from_pyproject()
