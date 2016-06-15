#!/usr/bin/env python
"""Installs extra dependencies.

Run all the install scripts in the proper order. After each install, the
corresponding environment vars are sourced to make sure the following dependencies
can be built properly.
"""

import json
import os


def main():
    """Main program."""
    # Load the dependencies data
    with open('dependencies.json') as f:
        dependencies = json.load(f)

    # Install each with an install_command line.
    for d in dependencies:
        install_command = d.get('install_command')
        if install_command is not None:
            print install_command
            os.system(install_command)


if __name__ == '__main__':
    main()
