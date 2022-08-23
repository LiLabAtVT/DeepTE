#!/usr/bin/env python3.7

####################################################################################
##this script is to generate seqs and labels that will be ready for the CNN training
def load_data(fname):
    seqs = []
    labels = []
    f = open(fname)
    for line in f:
        line_no_wspace = line.replace(" ","")
        line_no_nwline = line_no_wspace.replace("\n","")
        line_arr = line_no_nwline.split(",")
        label = line_arr[0]
        seq = line_arr[1]
        # sequence cleaning
        seq = seq.upper()    # b/c rep matrix built on uppercase
        seq = seq.replace("\t","")      # present in promoter

        seq = seq.replace("Y","C")  # undetermined nucleotides in splice
        seq = seq.replace("D","G")
        seq = seq.replace("S","C")
        seq = seq.replace("R","G")
        seq = seq.replace("V","A")
        seq = seq.replace("K", "G")
        seq = seq.replace("N", "T")
        seq = seq.replace("H", "A")
        seq = seq.replace("W", "A")
        seq = seq.replace("M", "C")
        seq = seq.replace("X", "G")
        seq = seq.replace("B", "C")
        #####
        labels.append(label)
        seqs.append(seq)
    f.close()
    return seqs, labels



