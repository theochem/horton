#!/usr/bin/env python
# -*- coding: utf-8 -*-
# HORTON: Helpful Open-source Research TOol for N-fermion systems.
# Copyright (C) 2011-2017 The HORTON Development Team
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

import json
from string import Template
import textwrap


def mywrap(s):
    lines = textwrap.wrap(s, 60, initial_indent='', subsequent_indent=' '*16,
                          break_on_hyphens=False)
    return ' \\\n'.join(lines)


def format_packages(key, dependencies, result):
    names = []
    names_doc = []
    names_dev = []
    for d in dependencies:
        name = d.get(key)
        if name is not None:
            if d['purpose'] == 'doc':
                names_doc.append(name)
            elif d['purpose'] == 'dev':
                names_dev.append(name)
            else:
                names.append(name)
    names = ' '.join(names)
    names_doc = ' '.join(names_doc)
    names_dev = ' '.join(names_dev)

    if len(names) > 0:
        result['dependencies_%s' % key] = mywrap(names)
    if len(names_doc) > 0:
        result['doc_dependencies_%s' % key] = mywrap(names_doc)
    if len(names_dev) > 0:
        result['dev_dependencies_%s' % key] = mywrap(names_dev)


def format_list(suffix, line_formatter, dependencies, result):
    lines = []
    lines_doc = []
    lines_dev = []
    for d in dependencies:
        line = line_formatter(d)
        if line is not None:
            if d['purpose'] == 'doc':
                lines_doc.append(line)
            elif d['purpose'] == 'dev':
                lines_dev.append(line)
            else:
                lines.append(line)
    result['dependencies_%s' % suffix] = '\n'.join(lines)
    result['doc_dependencies_%s' % suffix] = '\n'.join(lines_doc)
    result['dev_dependencies_%s' % suffix] = '\n'.join(lines_dev)


def formatter_rst(d):
    line = '* %s %s: %s' % (d['title'], d['version_requirement'], d['homepage'])
    comment = d.get('comment')
    if comment is not None:
        line += ' (%s)' % comment
    return line


def formatter_macports_rst(d):
    macports_name = d.get('macports_name')
    if macports_name is not None:
        macports_url = d['macports_url']
        return '* ``%s``: %s' % (macports_name, macports_url)


def formatter_macports_command(d):
    macports_name = d.get('macports_name')
    if macports_name is not None:
        macports_command = d['macports_command']
        return '    %s' % (macports_command)


def dict_raise_on_duplicates(ordered_pairs):
    """Reject duplicate keys when loading JSON data."""
    d = {}
    for k, v in ordered_pairs:
        if k in d:
           raise ValueError("duplicate key: %r" % (k,))
        else:
           d[k] = v
    return d


def collect_fields():
    # Load the dependencies data
    with open('../dependencies.json') as f:
        dependencies = json.load(f, object_pairs_hook=dict_raise_on_duplicates)

    result = {}

    # dependencies_rst and dep_dependencies_rst
    format_list('rst', formatter_rst, dependencies, result)

    # dependencies_* and doc_dependencies_*
    format_packages('fedora_25_rpm', dependencies, result)
    format_packages('fedora_24_rpm', dependencies, result)
    format_packages('fedora_24_pip', dependencies, result)
    format_packages('fedora_22_rpm', dependencies, result)
    format_packages('fedora_22_pip', dependencies, result)
    format_packages('ubuntu_16_deb', dependencies, result)
    format_packages('ubuntu_16_pip', dependencies, result)
    format_packages('ubuntu_15_deb', dependencies, result)
    format_packages('ubuntu_15_pip', dependencies, result)
    format_packages('macports_pip', dependencies, result)

    # macports stuff
    format_list('macports_rst', formatter_macports_rst, dependencies, result)
    format_list('macports_command', formatter_macports_command, dependencies, result)


    # turn dependencies into a dictonary to get a few final things
    dependencies = dict((d['name'], d) for d in dependencies)
    result['custom_install_libxc'] = dependencies['libxc']['install_command']
    result['custom_install_libint'] = dependencies['libint']['install_command']

    return result


big_fat_warning = """\
..
    : THIS FILE IS AUTOMATICALLY GENERATED. CHANGES TO THIS FILE WILL BE OVERWRITTEN
    : WHEN REBUILDING THE DOCUMENTATION. MAKE CHANGES IN
    :     %s
    : OR
    :     ../dependencies.json
    : INSTEAD.

"""


def substitute_file(fields, fn):
    assert fn.endswith('.template')
    fields['big_fat_warning'] = big_fat_warning % fn
    with open(fn) as f:
        t = Template(f.read())
    with open(fn[:-9], 'w') as f:
        f.write(t.safe_substitute(fields))


def main():
    fields = collect_fields()
    substitute_file(fields, 'user_download_and_install_linux.rst.template')
    substitute_file(fields, 'user_download_and_install_mac.rst.template')
    substitute_file(fields, 'user_download_and_install_windows.rst.template')
    substitute_file(fields, 'tech_dev_git.rst.template')


if __name__ == '__main__':
    main()
