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
"""Simulation tool for the trapdoor workflow on Travis.

It may be worth switching to pygit2 in future. (It seems better but does not work properly
on fedora 23.)
"""

import argparse
import os
import shutil
import subprocess

import git


class RepoError(Exception):
    """Raised when interaction with GIT repo fails."""

    pass


class Log(object):
    """A tiny logger object."""

    def __init__(self, name):
        """Initialize Log instance.

        Parameters
        ----------
        name : str
               The program name.
        """
        self._name = name

    def __call__(self, message, indent=4):
        """Print a message on screen."""
        print '/%s/ %s%s' % (self._name, ' '*indent, message)

    def section(self, name):
        """Dectorator to add section output to function.

        Parameters
        ----------
        name : str
               The function name.
        """
        def decorator(fn):
            """A function that returns the wrapped version of the original function."""
            def wrapper(*args, **kwargs):
                """The wrapper with additional printing."""
                self('BEGIN %s.' % name, indent=2)
                result = fn(*args, **kwargs)
                self('END %s.' % name, indent=2)
                return result
            return wrapper
        return decorator


log = Log('SIMULATE_TRAPDOOR_PR')


def main():
    """Main program."""
    # A few general things.
    args = parse_args()
    repo = git.Repo('.')

    # Get the QAWORKDIR. Create if it does not exist yet.
    qaworkdir = os.getenv('QAWORKDIR', 'qaworkdir')
    if not os.path.isdir(qaworkdir):
        os.makedirs(qaworkdir)

    # Pre-flight checks.
    orig_head_name, merge_head_name = run_pre_flight_checks(repo)

    try:
        make_temporary_merge(repo, merge_head_name)
        trapdoor_workflow(repo, args.script, qaworkdir, args.skip_ancestor)
    finally:
        roll_back(repo, orig_head_name, merge_head_name)


def parse_args():
    """Parse command-line arguments."""
    parser = argparse.ArgumentParser(description='Simulate trapdoor test locally.')
    parser.add_argument('script', type=str, metavar='trapdoor', help='Path to trapdoor script.')
    parser.add_argument('-s', '--skip-ancestor', default=False, action='store_true',
                        help='Do not run the trapdoor on master and re-use result for '
                             'ancestor from previous run.')
    return parser.parse_args()


@log.section('pre-flight checks')
def run_pre_flight_checks(repo):
    """Run some initial checks before doing anything.

    Parameters
    ----------
    repo : git.Repo
           A repository object from GitPython.
    """
    log('Check whether all changes are committed.')
    # (Can be made better in future version by temporarily auto-committing.)
    if repo.is_dirty():
        raise RepoError('Not all changes are committed.')

    log('Check whether master is up to date with origin/master.')
    remote_refs = git_ls_remote('origin')
    if remote_refs['refs/heads/master'] != repo.heads.master.object.hexsha:
        raise RepoError('Master is not up to date.')

    log('Check whether master is currently _not_ checked out.')
    if repo.active_branch == repo.heads.master:
        raise RepoError('The master branch is checked out.')

    log('Check whether temporary branch name does not exist yet.')
    orig_head_name = repo.active_branch.name
    merge_head_name = '%s-trapdoor-tmp-merge' % orig_head_name
    if merge_head_name in repo.heads:
        raise RepoError('The branch %s exists, probably due to earlier failures of this program.')

    return orig_head_name, merge_head_name


def git_ls_remote(url):
    """Run git ls-remote with GitPython and parse output.

    Parameters
    ----------
    url : str
          The remote repository URL.

    Returns
    -------
    refs : dict
           A dictionary with (ref, shahex) combinations.
    """
    remote_refs = {}
    git_interface = git.cmd.Git()
    for ref in git_interface.ls_remote(url).split('\n'):
        hash_ref_list = ref.split('\t')
        remote_refs[hash_ref_list[1]] = hash_ref_list[0]
    return remote_refs


@log.section('temporary merge')
def make_temporary_merge(repo, merge_head_name):
    """Make a temporary merge to run trapdoor_*.py feature.

    Parameters
    ----------
    repo : git.Repo
           A repository object from GitPython.
    merge_head_name : str
                      The name of the branch in which the temporary merge will be made.
                      This branches of from the original branch when the simulator was
                      called.
    """
    log('Make a new branch: %s.' % merge_head_name)
    merge_head = repo.create_head(merge_head_name)
    log('Checkout new branch: %s.' % merge_head_name)
    merge_head.checkout()
    log('Merge with master.')
    repo.index.merge_tree(repo.heads.master)
    repo.index.merge_tree(merge_head)
    log('Check if merge went well.')
    unmerged_blobs = repo.index.unmerged_blobs()
    for _path, list_of_blobs in unmerged_blobs.iteritems():
        for (stage, _blob) in list_of_blobs:
            # Now we can check each stage to see whether there were any conflicts
            if stage != 0:
                raise RepoError('Merge with master failed.')
    log('Commit merge.')
    repo.index.commit('Temporary merge for trapdoor simulation', [
        repo.heads.master.commit,
        merge_head.commit
    ])


@log.section('trapdoor workflow')
def trapdoor_workflow(repo, script, qaworkdir, skip_ancestor):
    """Run the trapdoor scripts in the right order.

    Parameters
    ----------
    repo : git.Repo
           A repository object from GitPython.
    script : str
             The relative path to the trapdoor script.
    qaworkdir : str
                The location of the QA work directory.
    skip_ancestor : bool
                    If True, the trapdoor script is not executed in the ancestor.
    """
    subprocess.call([script, 'feature'])
    if skip_ancestor:
        subprocess.call([script, 'report'])
    else:
        copied_script = os.path.join(qaworkdir, os.path.basename(script))
        shutil.copy(script, copied_script)
        shutil.copy('tools/qa/trapdoor.py', os.path.join(qaworkdir, 'trapdoor.py'))
        # Check out the master branch. (We should be constructing the ancestor etc. but
        # that should come down to the same thing for a PR.)
        repo.heads.master.checkout()
        subprocess.call([copied_script, 'ancestor'])
        subprocess.call([copied_script, 'report'])


@log.section('roll back')
def roll_back(repo, orig_head_name, merge_head_name):
    """Revert to the original head and clean up.

    Parameters
    ----------
    repo : git.Repo
           A repository object from GitPython.
    orig_head_name : str
                     The name of the original feature branch.
    merge_head_name : str
                      The name of the branch in which the temporary merge will be made.
                      This branches of from the original branch when the simulator was
                      called.

    """
    log('Remove lingering changes.')
    repo.head.reset(index=True, working_tree=True)
    log('Checkout the original head.')
    repo.heads[orig_head_name].checkout()
    # Check if temporary branch exsists
    if merge_head_name in repo.heads:
        log('Remove the temporary merge (with -D).')
        repo.delete_head(merge_head_name, '-D')
    else:
        log('No temporary merge found.')


if __name__ == '__main__':
    main()
