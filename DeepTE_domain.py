#!/usr/bin/env python

##DeepTE_domain is to detect TE domains used for improving performance of DeepTE

##BUILT-IN MODULES
import argparse
import sys
import subprocess
from distutils.spawn import find_executable

from scripts import DeepTE_detect_domain as domain_detection


def get_parsed_args():

    parser = argparse.ArgumentParser(description="DeepTE_domain detect TE domains")

    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory to store intermediate files of "
                                                                     "each step. Default: ./ ")

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files. "
                                                                    "Default: ./ ")

    parser.add_argument("-i", dest='ipt_seq', help="Input sequences that are unknown TEs or DNA sequences")

    parser.add_argument("-s", dest='supfile_dir', help="provide supplementary dir that contains required files.")

    parser.add_argument("--hmmscan", dest="hmmscan", default="/usr/bin/hmmscan", help="File path to hmmscan executable"
                                                                                      "Default: /usr/bin/hmmscan")


    ##parse of parameters
    args = parser.parse_args()
    return args


def main(argv=None):

    if argv is None:
        argv = sys.argv
    args = get_parsed_args()

    ##Check whether the files are provided
    if args.ipt_seq is None:
        print('Cannot find input sequence, please provide the file!')
        return  ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.ipt_seq, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return

    if args.supfile_dir is None:
        print ('Cannot find supplementary dir, please provide that')
        return

    hmmscan_exe = ''
    if find_executable(args.hmmscan) is not None:
        hmmscan_exe = args.hmmscan
        print('Hmmscan executable can be found')
    else:
        print("Cannot find Hmmscan executable, please check if it has been installed.")
        return

    ###########################################
    ##create the working and output directories
    working_dir = args.working_dir
    if not working_dir.endswith('/'):
        working_dir = working_dir + '/'
    else:
        working_dir = working_dir

    output_dir = args.output_dir
    if not output_dir.endswith('/'):
        output_dir = output_dir + '/'
    else:
        output_dir = output_dir


    ########################
    ##set the required files
    ipt_seq = args.ipt_seq
    supfile_dir = args.supfile_dir

    ##cp the minifam to the current output dir
    cmd = 'cp -r ' + supfile_dir + '/minifam* ' + working_dir
    subprocess.call(cmd, shell=True)

    ##detect domain
    te_domain_pattern_dic = domain_detection.collect_domain_from_six_rdframe(ipt_seq, hmmscan_exe, working_dir)
    ##write out results
    with open(output_dir + '/opt_te_domain_pattern.txt', 'w+') as opt_dm_pt:
        for eachnm in te_domain_pattern_dic:
            opt_dm_pt.write(eachnm + '\t' + te_domain_pattern_dic[eachnm] + '\n')



if __name__ == "__main__":
    main()