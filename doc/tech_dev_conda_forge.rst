..
    : HORTON: Helpful Open-source Research TOol for N-fermion systems.
    : Copyright (C) 2011-2022 The HORTON Development Team
    :
    : This file is part of HORTON.
    :
    : HORTON is free software; you can redistribute it and/or
    : modify it under the terms of the GNU General Public License
    : as published by the Free Software Foundation; either version 3
    : of the License, or (at your option) any later version.
    :
    : HORTON is distributed in the hope that it will be useful,
    : but WITHOUT ANY WARRANTY; without even the implied warranty of
    : MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    : GNU General Public License for more details.
    :
    : You should have received a copy of the GNU General Public License
    : along with this program; if not, see <http://www.gnu.org/licenses/>
    :
    : --

Making a conda-forge release
############################

Conda-forge is the recommended way to make a HORTON release.

If you have a commit you wish to release, please take the following steps:

1. Tag commit and push to github. This should build a tarball release.
2. Make a fork of the conda-forge feedstock: https://github.com/conda-forge/horton-feedstock.  
3. Update `recipe/meta.yml` to reflect the new version. Drop the build number back to 0. Be sure to update 
   the sha256 hash and version as well.
4. PR against the official conda-forge feedstock. Monitor the CI for failures.
5. If all checks pass, notify the maintainers. 
