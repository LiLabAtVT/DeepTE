#!/usr/bin/env python


##detect the TR and RT for all the TEs
##detect all the six reading frames
##if there is a reading frame showing RT without TR, that will be classed into the classI
##if TR exit, no RT exits, that will be classed into the classII
##if no domain exit, leave it out
##if EN exit, with out other domains, that will be LINE


##This script is to detect the domain pattern for each tes with domains
##This script is to detect the completeness of the pattern for each te in order to calculate the score

##################
##import functions
##################

##import function of translation
from Bio.SeqUtils import six_frame_translations
import re
from Bio import SeqIO
import subprocess



##############################################################################
##Step 1: Translate each reading frame and identify the possible reading frame
##############################################################################
##update the translate for six reading frame for the no mite DNA
##define a reverse complement functions
def getRC(te_seq_rc):
    out = te_seq_rc.reverse_complement()
    return out

##define a function to translate all the six reading frame
def translate_all_te (te_seq):
    ##initial a dictionary to store the protein sequence for the six reading frame
    pro_seq_dic = {}
    ##Runs a six-frame translation for each gene in the string ts2seq
    f6 = six_frame_translations(te_seq)
    ##Splits each gene in the above mentioned translation with a new line
    f6tmp = str(f6).split('\n')
    num_m = re.match('(.+)5$', str(len(f6tmp)))
    num_cnt = num_m.group(1)
    ####################################################################
    ##the following is the first three reading frame of this te sequence
    ##get the location information of sequences
    l1 = 5
    l2 = 6
    l3 = 7
    ##Initializes a hand full of lists we will later fill with the genomic data to analyze
    list1 = []
    list2 = []
    list3 = []
    ##Sets range that the genomic data we are looking at can be based on the 6-frame translation
    for i in range(0, (int(num_cnt))):
        l1 = 5 + i * 10
        list1.append(l1)
        l2 = 6 + i * 10
        list2.append(l2)
        l3 = 7 + i * 10
        list3.append(l3)
    ##extract the sequence information from the location informationa
    ##the following is the first three list
    allist_f = [list1, list2, list3]
    ##create a string to be the index for the key of the dictionary
    fr_count = 0
    for eachlist in allist_f:
        fr_count = fr_count + 1
        ##extrac the sequence information
        ##For each gene, adds the word 'list' for each additional gene after the start + end codon
        framenm = 'f' + str(fr_count)
        ##initial a list to store all the protein sequences
        f6aa = []
        # Setting the frame protein for each translation
        for i in eachlist:
            pseq = f6tmp[i]
            f6aa.append(pseq)
        ##combination of the sequence
        comptlist = ' '.join(f6aa)
        ##split the '' between each aa and then replace it with ,
        comptstr = comptlist.replace(' ', '')
        ##store sequence information on the dictionary
        pro_seq_dic[framenm] = comptstr

    ###############################################################
    ##the following is last three reading frame of this te sequence

    ##Inputs the data to get the reversed compliment
    seq_reverse = getRC(te_seq)
    f6 = six_frame_translations(seq_reverse)
    f6tmp = str(f6).split('\n')
    num_m = re.match('(.+)5$', str(len(f6tmp)))
    num_cnt = num_m.group(1)
    l4 = 5
    l5 = 6
    l6 = 7
    list4 = []
    list5 = []
    list6 = []
    for i in range(0, (int(num_cnt))):
        l4 = 5 + i * 10
        list4.append(l4)
        l5 = 6 + i * 10
        list5.append(l5)
        l6 = 7 + i * 10
        list6.append(l6)
    allist_f = [list4, list5, list6]
    re_count = 3
    for eachlist in allist_f:
        re_count = re_count + 1
        framenm = 'f' + str(re_count)
        f6aa = []
        for i in eachlist:
            pseq = f6tmp[i]
            f6aa.append(pseq)
        comptlist = ' '.join(f6aa)
        comptstr = comptlist.replace(' ', '')
        pro_seq_dic[framenm] = comptstr

    return (pro_seq_dic)


########################################################
##step 2 identify domain information for the six rdframe
########################################################
##initial a function to identify the domain for the six reading frame
##define a function to identify te domain and generate the domain table
##import the pro_seq_dic generated from the translate function
def iden_domain_six_rdframe (pro_seq_dic,seqid,hmmscan_exe,working_dir):

    ##initial a dictionary to store the results
    te_domain_six_rdframe_dic = {}

    ##initial a dictionary to store the name of six reading frame
    six_frame_nm_dic = {}
    for eachframe in pro_seq_dic:
        te_frame_nm = seqid + '_' + eachframe
        te_opt_nm = eachframe + '.seq'
        te_pro_seq = open(working_dir + '/' + te_opt_nm, 'w')
        ##this will generate each reading frame file
        te_pro_seq.write('>' + te_frame_nm + '\n' + str(pro_seq_dic[eachframe]))
        te_pro_seq.close()  ##it is important to close the te_pro_seq file
        six_frame_nm_dic[eachframe] = 1

    ##use the hmm to detect domain information for each frame
    for eachnm in six_frame_nm_dic:
        tefm_seq = eachnm + '.seq'
        opt_tefm = eachnm + '.out'
        cmd = hmmscan_exe + ' ' + working_dir + '/minifam ' + working_dir + '/' + tefm_seq + ' > ' + working_dir + '/' + opt_tefm
        subprocess.call(cmd, shell=True)

        ##open the output file from hmmscan
        domain_te_file = open(working_dir + '/' + opt_tefm, 'r')
        for line in domain_te_file:
            line.strip()
            if re.match('.+\d+\s+PF\d+_.+', line):
                mt = re.match('(.+PF\d+_[A-Z]+).+', line)
                new_line = mt.group(1)
                te_line = seqid + '_' + eachnm + '\t' + new_line
                te_domain_six_rdframe_dic[te_line] = 1

    return (te_domain_six_rdframe_dic)


################################################
##step 3 collect domain information into a table
################################################
def collect_domain_from_six_rdframe (input_ltr_fas,hmmscan_exe,working_dir):

    ##import the input_ltr_fas from the arg in the begin of the script
    ##initial a dic to store all the information about the te and domain
    te_domain_pattern_dic = {}  ##key is the te name and value is the pattern

    ##initial a te_count to calculate the number of analyzed te
    te_count = 0

    for seq_record in SeqIO.parse(input_ltr_fas, "fasta"):

        seqid = seq_record.id
        te_count = te_count + 1
        print('the number of analysis is ' + str(te_count))

        ##function translate to get the pro seq for each te
        pro_seq_six_frame_dic = translate_all_te(seq_record.seq)
        ##function iden_domain_six_rdframe to identify domain for each reading frame
        te_domain_six_rdframe_dic = iden_domain_six_rdframe(pro_seq_six_frame_dic, seqid, hmmscan_exe,working_dir)

        temp_rd_nm_dic = {}  ##={'f1':TR,TR,TR,'f2':TR,RT}
        store_rd_nm_list = []
        ori_nm = ''
        if len(te_domain_six_rdframe_dic) != 0:  ##if te_domain_six_rdframe_dic has no content we do not consider

            dm_string = ''
            for eachline in te_domain_six_rdframe_dic:
                col = eachline.strip().split()
                ##get the original name
                mt = re.match('(.+)_(f\d)',col[0])
                ori_nm = mt.group(1)
                rd_nm = mt.group(2)
                ##get the domain name
                mt = re.match('.+_(.+)',col[9])
                dm_nm = mt.group(1)

                if rd_nm not in store_rd_nm_list:
                    store_rd_nm_list.append(rd_nm)
                    dm_string = dm_nm
                    temp_rd_nm_dic[rd_nm] = dm_string
                else:
                    dm_string = dm_string + ',' + dm_nm
                    temp_rd_nm_dic[rd_nm] = dm_string

            ##generate a final line
            final_line_list = [] ##contain f1:TR,TR
            for eachrd_nm in temp_rd_nm_dic:
                eachrd_string = eachrd_nm + ':' + temp_rd_nm_dic[eachrd_nm]
                final_line_list.append(eachrd_string)

            final_line = ';'.join(final_line_list)
            te_domain_pattern_dic[ori_nm] = final_line

    return (te_domain_pattern_dic)


































