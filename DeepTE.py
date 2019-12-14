#!/usr/bin/env python

##updation 12.14 add argument to download model automatically

##updation 9.26 distinguish the TE and no TE function

##DeepTE is to classify unknown TEs into different TE orders and families

##BUILT-IN MODULES
import os
import argparse
import sys
import time
import subprocess

from scripts import DeepTE_pipeline_no_modification as pipeline_no_m
from scripts import DeepTE_pipeline_yes_modification as pipeline_yes_m
from scripts import DeepTE_pipeline_unknown_sequence as pipeline_uns
from scripts import DeepTE_generate_CNN_dataset as generate_dataset
from scripts import DeepTE_combine_opt as combine_opt

def get_parsed_args():

    parser = argparse.ArgumentParser(description="DeepTE classify TEs via neural network")

    ##require files
    parser.add_argument("-d", dest='working_dir', default="./", help="Working directory to store intermediate files of "
                                                                     "each step. Default: ./ ")

    parser.add_argument("-o", dest='output_dir', default="./", help="Output directory to store the output files. "
                                                                    "Default: ./ ")

    parser.add_argument("-i", dest='ipt_seq', help="Input sequences that are unknown TEs or DNA sequences")



    parser.add_argument("-sp", dest='sp_type',help="Provide species type of the input sequence")

    ##optional files

    ##updation 12.14 add an argument to allow users to download model dir
    parser.add_argument("-m", dest='model_name', help="Provide one of model names: -m P or -m M or -m F or -m O or -m U."
                                                      "This argument will directly download the model dir,"
                                                      "So users do not need to initiate -m_dir.")

    parser.add_argument("-m_dir", dest='model_dir',help="Provide model_dir that could be downloaded from website"
                                                        "If users already download the model dir and can directly use this argument to indicate the model directory")


    parser.add_argument("-fam", dest='te_fam',help="Provide TE family name for the input te sequence"
                                                   "Default: All"
                                                   "All: the input sequence is unknown TEs"
                                                   "ClassI: the input sequence is ClassI TEs"
                                                   "ClassII: the input sequence is ClassII_sub1 TEs that do not contain Helitron TEs"
                                                   "LTR: the input sequence is LTR TEs"
                                                   "nLTR: the input sequence is nLTR TEs"
                                                   "LINE: the input sequence is LINE TEs"
                                                   "SINE: the input sequence is SINE TEs"
                                                   "Domain: the input sequence is Class II TEs with specified super families")

    parser.add_argument("-modify", dest='domain_file', help="If users set this argument, "
                                                            "users need to provide domain file"
                                                            " generated from another script")
    ##updation 9.26 distinguish TEs and no TE function
    parser.add_argument('-UNS', dest='yes', help="If users set this argument, "
                                                        "users need change the -i to the the DNA sequences."
                                                        "This function will classify the sequences into TEs, CDS, or Intergenic sequences,"
                                                        "-sp and -fam do not need to provide" )


    ##parse of parameters
    args = parser.parse_args()
    return args

##updation 12.14
##define a function to indicate where to download model
def download_model (working_dir,model_name,google_drive_path):

    ##create a dir in the working_dir to store the model dir
    download_model_dir = working_dir + '/download_' + model_name + '_model_dir'
    if not os.path.exists(download_model_dir):
        os.makedirs(download_model_dir)

    cmd = "wget --load-cookies /tmp/cookies.txt \"https://docs.google.com/uc?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate " + \
          "\'https://docs.google.com/uc?export=download&id=" + google_drive_path + "\' -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\\1\\n/p')&id=" + google_drive_path + "\" -O " + \
          download_model_dir + "/" + model_name + "_model.tar.gz " + "&& rm -rf /tmp/cookies.txt"
    subprocess.call(cmd, shell=True)

    cmd = 'tar -xvzf ' +  download_model_dir + '/' + model_name + '_model.tar.gz -C ' + download_model_dir + '/'
    subprocess.call(cmd,shell=True)

    target_download_model_dir = download_model_dir + '/*_model'

    return (target_download_model_dir)


def main(argv=None):

    if argv is None:
        argv = sys.argv
    args = get_parsed_args()


    ######################################
    ##Check whether the files are provided
    if args.ipt_seq is None:
        print ('Cannot find input sequence, please provide the file!')
        return ##import to close the script if not find the te lib
    else:
        try:
            file = open(args.ipt_seq, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the file!')
            return


    if args.domain_file is not None:
        ##check the file
        try:
            file = open(args.domain_file, 'r')  ##check if the file is not the right file
        except IOError:
            print('There was an error opening the domain_file!')
            return

    if args.model_name is not None and args.model_dir is not None:
        print('Please only provide model_name \'-m\' or model_dir \'-m_dir\'')
        return
    else:
        if args.model_name is None and args.model_dir is None:
            print('Cannot find model dir, please provide the dir in \'-m_dir\' or the model name in \'-m\'!')
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


    ##updation 12.14
    ##check download model directory
    download_model_dir = ''
    if args.model_name is not None:
        model_name = args.model_name
        if model_name != 'P' and model_name != 'M' and model_name != 'F' and model_name != 'O':
            print("Please use one of 'P', 'M', 'F', 'O' to be model name")
            return
        else:
            ##download P model
            if model_name == 'P':
                ##create a dir in the working_dir to store the model dir
                google_drive_path = '1voj86STKcQH8lAhvY6yl5E65nzaM6o0B'
                download_model_dir = download_model(working_dir, model_name, google_drive_path)

            if model_name == 'M':
                ##create a dir in the working_dir to store the model dir
                google_drive_path = '1ExRwC3szJ4XMa3ikxM9Ccu31lY79rdw9'
                download_model_dir = download_model(working_dir, model_name, google_drive_path)

            if model_name == 'F':
                ##create a dir in the working_dir to store the model dir
                google_drive_path = '1uvnm99ypauIKtqCxoybdtT-mEMdoupip'
                download_model_dir = download_model(working_dir, model_name, google_drive_path)

            if model_name == 'O':
                ##create a dir in the working_dir to store the model dir
                google_drive_path = '1Q6HW1NhNs0a6Ykrw7jGEKKPWxawpWiuM'
                download_model_dir = download_model(working_dir, model_name, google_drive_path)

            if model_name == 'U':
                ##create a dir in the working_dir to store the model dir
                google_drive_path = '17JswzJkTakF4mpr1MObyTp0CNxXodNZK'
                download_model_dir = download_model(working_dir, model_name, google_drive_path)

    ##if users do not notify the model_name instead they want to use downloaded model they need provide model_dir
    else:
        if args.model_dir is None:
            print('Cannot find model dir, please provide the dir in \'-m_dir\' or the model name in \'-m\'!')
            return
        else:
            download_model_dir = args.model_dir


    #################################
    ##set the input model dir and seq
    model_dir = download_model_dir
    ipt_seq = args.ipt_seq
    sp_type = args.sp_type

    ##default
    if args.te_fam is not None:
        te_fam = args.te_fam
    else:
        te_fam = 'All'

    ##generate temp output in the working_dir
    temp_store_opt_dir = working_dir + '/store_temp_opt_dir'
    if not os.path.exists(temp_store_opt_dir):
        os.makedirs(temp_store_opt_dir)

    ##########################################
    ##transfer fasta data to CNN data inputset
    print('Step1: transfer fasta data to CNN input data')
    final_format_dic = generate_dataset.change_format_for_ncc(ipt_seq)
    final_format_line_dic = generate_dataset.generate_target_line(final_format_dic)

    with open(working_dir + '/opt_input_CNN_data.txt', 'w+') as opt:
        for eachid in final_format_line_dic:
            opt.write(final_format_line_dic[str(eachid)] + '\n')

    input_CNN_data_file = working_dir + '/opt_input_CNN_data.txt'

    ##updation 9.26
    ##if users call UNS model
    if args.yes is None:
        ##If users call domain argument
        print('Step2: classify TEs')
        if args.domain_file is not None:
            print('Step2: 1) domain information is exist')
            domain_file = args.domain_file
            te_domain_pattern_dic = pipeline_yes_m.store_domain_pattern_infor(domain_file)
            pipeline_yes_m.classify_pipeline(model_dir, input_CNN_data_file, temp_store_opt_dir, sp_type,te_domain_pattern_dic,te_fam)

        ##If users do not call domain argument
        else:
            print('Step2: 2) domain information is not exist')
            ##run the DeepTE_pipeline
            pipeline_no_m.classify_pipeline(model_dir, input_CNN_data_file, temp_store_opt_dir, sp_type,te_fam)

        ##write out final results
        print('Step3: generate final output')
        ##write out the name file
        combine_opt.extract_combine_infor(temp_store_opt_dir, output_dir)
        ##write out the fasta file
        combine_opt.generate_fasta(ipt_seq, output_dir + '/opt_DeepTE.txt', output_dir)

    else:
        print('Step2: classify unknown sequences')
        pipeline_uns.classify_pipeline(model_dir, input_CNN_data_file, temp_store_opt_dir)
        ##write out the name and fasta files
        combine_opt.generate_fasta_UNS(ipt_seq, temp_store_opt_dir, output_dir)







if __name__ == "__main__":

    start_time = time.time()
    print('start time is ' + str(start_time))
    main()

    print("--- %s seconds ---" % (time.time() - start_time))








