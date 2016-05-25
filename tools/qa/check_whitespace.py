#!/usr/bin/env python
"""Script that checks for white space errors in every commit from master to feature.

This script ignores and command-line arguments. Just don't provide any.
"""

import subprocess
import sys


def main():
    """Run ``git diff --check`` on every relevant commit."""
    # Get the common ancestor with the master branch
    command = ['git', 'merge-base', 'master', 'HEAD']
    ancestor = subprocess.check_output(command).strip()

    # Get the list of commit ids and descriptions between the current and master branch.
    # The latest one is printed first and the HEAD of the master branch is included
    command = ['git', 'log', '%s^..HEAD' % ancestor, '--pretty=oneline', '--color=never']
    commits_str = subprocess.check_output(command)

    # Parse the output
    commits = []
    for line in commits_str.split('\n')[:-1]:
        pos = line.find(' ')
        commit_id = line[:pos]
        message = line[pos + 1:]
        commits.append((commit_id, message))

    # Loop over all commits and check the diffs
    error_count = 0
    for icommit in xrange(len(commits) - 1):
        print 'Checking whitespace in %s %s' % commits[icommit]
        command = ['git', '-c', 'core.whitespace=tab-in-indent', 'diff', commits[icommit + 1][0], commits[icommit][0], '--check']
        proc = subprocess.Popen(command)
        retcode = proc.wait()
        if retcode != 0:
            error_count += 1

    # return code
    if error_count > 0:
        sys.exit(1)


if __name__ == '__main__':
    main()
