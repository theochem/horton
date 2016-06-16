#!/usr/bin/env python
"""Generates documentation and ammends it to the gh-pages branch.

This script can manage different documentation versions and generates an index.html with
the directory of different documentation versions. The style of the index is similar to
the read_the_docs theme in Sphinx.
"""

from distutils.version import StrictVersion
import os
from glob import glob
import json
import shutil
import subprocess


css_template = r"""
body {
    font-family: "Lato",sans-serif;
    background-color: #edf0f2;
    font-size: 16px;
}
a {
    color: #2980B9;
    text-decoration: none;
    cursor: pointer;
}
a:visited {
    color: #9B59B6;
}
h1 {
    font-size: 28px;
    font-weight: 700;
    font-family: "Roboto Slab",sans-serif;
}
h2 {
    font-size: 24px;
    font-weight: 700;
    font-family: "Roboto Slab",sans-serif;
}
span.fineprint {
    font-size: 12px;
    color: #999999;
}
"""


html_template = r'''<!doctype html>

<html lang="en">
<head>
  <meta charset="utf-8">
  <title>HORTON Documentation Directory</title>
  <meta name="author" content="Toon Verstraelen">
  <link href='https://fonts.googleapis.com/css?family=Roboto+Slab:400,700' rel='stylesheet' type='text/css'>
  <link href='https://fonts.googleapis.com/css?family=Lato' rel='stylesheet' type='text/css'>
<style>
{css}
</style>
</head>
<body>
<h1>HORTON Documentation Directory</h1>
<h2>Stable releases</h2>
<ul>{stable}</ul>
<h2>Beta releases (feature freeze)</h2>
<ul>{beta}</ul>
<h2>Alpha releases (features being added)</h2>
<ul>{alpha}</ul>
</body>
</html>
'''


release_template = r'''
<li><a href="{version}/index.html">HORTON {version}</a> released on {release_date}.<br \>
<span class='fineprint'>Documentation built on {doc_build_date} from HORTON {describe} ({doc_release_date}) .</span>
</li>
'''

redirect_template = r'''<!doctype html>

<html lang="en">
<head>
  <meta charset="utf-8">
  <title>Latest HORTON Documentation (redirect)</title>
  <meta http-equiv="Refresh" content="2; url=../{version}/index.html">
  <meta name="author" content="Toon Verstraelen">
  <link href='https://fonts.googleapis.com/css?family=Roboto+Slab:400,700' rel='stylesheet' type='text/css'>
  <link href='https://fonts.googleapis.com/css?family=Lato' rel='stylesheet' type='text/css'>
<style>
{css}
</style>
</head>
<body>
<h2>Redirecting to the latest stable version of the HORTON documentation ...</h2>
<p>... or you can click <a href="../{version}/index.html">here</a>
</html>
'''


def version_sort_key(releaseinfo):
    """Return a sort key for versioninfo dictionaries."""
    return StrictVersion(releaseinfo['version'])


def format_releases(releases):
    """Format a list of releases into HTML bullet points.

    Parameters
    ----------
    versions : list of str
               The list of versions.
    """
    version_lines = [release_template.format(**releaseinfo) for releaseinfo in releases]
    return ''.join(version_lines)


def update_index():
    """Update index.html (the documentation directory)."""
    releases = {
        'stable': [],
        'beta': [],
        'alpha': [],
    }

    for fn_ri in sorted(glob('*.*.*/releaseinfo.json')):
        with open(fn_ri) as f:
            releaseinfo = json.load(f)
        print 'Found', releaseinfo['version']
        if 'b' in releaseinfo['version']:
            releases['beta'].append(releaseinfo)
        elif 'a' in releaseinfo['version']:
            releases['alpha'].append(releaseinfo)
        else:
            releases['stable'].append(releaseinfo)
    for values in releases.itervalues():
        values.sort(key=version_sort_key, reverse=True)

    # Write the index file
    subs = dict((key, format_releases(versions))
                for key, versions in releases.iteritems())
    subs['css'] = css_template
    with open('index.html', 'w') as f:
        f.write(html_template.format(**subs))

    # Write the redirect file.
    if len(releases['stable']) > 0:
        subs = releases['stable'][0].copy()
        subs['css'] = css_template
        if not os.path.isdir('latest'):
            os.mkdir('latest')
        with open('latest/index.html', 'w') as f:
            f.write(redirect_template.format(**subs))


def call(command, cwd=None):
    """Call a command and print it on screen.

    Parameters
    ----------
    command : list
              A list of command line arguments, the first being the executable.
    cwd : str
          When given, the script is executed in the directory `cwd`.
    """
    command_str = ' '.join(('\'%s\'' % arg) if (' ' in arg) else arg for arg in command)
    if cwd is None:
        print command_str
    else:
        print command_str, '(in %s)' % cwd
    return subprocess.check_output(command, cwd=cwd)


def main():
    """Main program."""
    # Get the version string from the git tag
    version = call(['git', 'describe', '--tags'])
    version = version[:version.find('-')]

    # de-cruft
    call(['./cleanfiles.sh'])

    # Build the docs
    call(['make', 'html'], cwd='doc')

    # Copy releas info file
    shutil.copy('doc/releaseinfo.json', 'doc/_build/html/releaseinfo.json')

    # Copy the docs to the gh-pages branch, removing the old copy for the current version
    call(['git', 'checkout', 'gh-pages'])
    if os.path.exists(version):
        call(['git', 'rm', '-r', version])
    os.mkdir(version)
    call(['rsync', '-av', 'doc/_build/html/', version])

    # Update the documentation directory
    update_index()

    # Add new files to commit
    call(['git', 'add', version, 'latest', 'index.html'])

    # Commit
    call(['git', 'commit', '-a', '--amend', '-m', 'Automatic documentation update', '--no-verify'])


if __name__ == '__main__':
    main()
