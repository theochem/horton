#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2016 The HORTON Development Team
#
# This file is part of HORTON.
#
# HORTON is free software; you can redistribute it and/or
# modify it under the terms of the GNU General Public License
# as published by the Free Software Foundation; either version 3
# of the License, or (at your option) any later version.
#
# HORTON is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License
# along with this program; if not, see <http://www.gnu.org/licenses/>
#
# --
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
