#!/usr/bin/env python
#!/usr/bin/python2.7
# vim:ft=python
#
# Simple program to run CTest on a branch and submit the results to 
# CDash.
#

import os, sys, string
import subprocess
from argparse import ArgumentParser

dashboard_dir='/home/suzanka/Shared/PhD/gimic/cdash'
cdash_script='CDashTestBranch.cmake'
project_repository="git@repo.ctcc.no:gimic.git"
project_name="GIMIC"

parser = ArgumentParser(description="Setup build configurations.")
parser.add_argument('--repository',
        action='store',
        default=project_repository,
        help='set project repository [default: %(default)s]')
parser.add_argument('--dashboard-dir',
        action='store',
        default=dashboard_dir,
        help='directory where the CTest scripts are located [default: %(default)s]')
parser.add_argument('--project-name',
        action='store',
        default=project_name,
        help='set project name [default: %(default)s]')
parser.add_argument('--build',
        action='store',
        default='Debug',
        help='set build configuration [default: %(default)s]')
parser.add_argument('--model',
        action='store',
        default='Experimental',
        help='set testing model [default: %(default)s]')
parser.add_argument('--branch',
        action='store',
        default='master',
        help='set branch [default: %(default)s]')
parser.add_argument('--site',
        action='store',
        default=None,
        help='set site name')
parser.add_argument('--cc',
        action='store',
        default=None,
        help='set C compiler name')
parser.add_argument('--cxx',
        action='store',
        default=None,
        help='set C++ compiler name')
parser.add_argument('--fc',
        action='store',
        default=None,
        help='set Fortran compiler name')
parser.add_argument('--toolchain-name',
        action='store',
        default=None,
        help='set toolchain [default: %(default)s]')
parser.add_argument('--mpi',
        action='store',
        type=int,
        default=0,
        metavar='N',
        help='run tests using N MPI procs [default: %(default)s]')
parser.add_argument('--omp',
        action='store',
        type=int,
        default=0,
        metavar='N',
        help='run tests using N OpenMP threads [default: %(default)s]')
parser.add_argument('--coverage',
        action='store_true',
        default=False,
        help='Switch on code coverage [default: %(default)s]')
parser.add_argument('--memcheck',
        action='store_true',
        default=False,
        help='Run valgrind [default: %(default)s]')
parser.add_argument('--dry-run',
        action='store_true',
        default=False,
        help='Only show commands without executing [default: %(default)s]')

args = parser.parse_args()

def main():
    script = "{0}/{1}".format(dashboard_dir, cdash_script)

    if not os.environ.has_key('PROJECT_NAME') or args.project_name:
        os.environ['PROJECT_NAME'] = project_name
    if not os.environ.has_key('PROJECT_REPOSITORY') or args.project_repository:
        os.environ['PROJECT_REPOSITORY'] = project_repository
    if not os.environ.has_key('DASHBOARD_DIR') or args.dashboard_dir:
        os.environ['DASHBOARD_DIR'] = dashboard_dir

    ctest = "ctest -V " 
    ctest += "-C {0} ".format(args.build)
    ctest += "-S {0}".format(script)
    ctest += ",model={0}".format(args.model)
    ctest += ",branch={0}".format(args.branch)
    ctest += ",coverage={0}".format(args.coverage)
    ctest += ",memcheck={0}".format(args.memcheck)
    if args.site:
        ctest += ",site={0}".format(args.site)
    if args.toolchain_name:
        ctest += ",toolchain_name={0}".format(args.toolchain_name)
    if args.mpi > 0:
        ctest += ",mpi={0}".format(args.mpi)
    if args.omp > 0:
        ctest += ",omp={0}".format(args.omp)
    if args.cc:
        ctest += ",cc={0}".format(args.cc)
    if args.cxx:
        ctest += ",cxx={0}".format(args.cxx)
    if args.fc:
        ctest += ",fc={0}".format(args.fc)
    
    if not args.cc and os.environ.has_key('CC'):
        ctest += ",cc={0}".format(os.environ['CC'])
    if not args.cxx and os.environ.has_key('CXX'):
        ctest += ",cxx={0}".format(os.environ['CXX'])
    if not args.fc and os.environ.has_key('FC'):
        ctest += ",fc={0}".format(os.environ['FC'])

    print ctest
    s = 0
    if not args.dry_run:
        p = subprocess.Popen(ctest,
                shell=True,
                stdin=subprocess.PIPE,
                stdout=subprocess.PIPE)
        s = p.communicate()[0]
    return s

if __name__ == '__main__':
    main()
