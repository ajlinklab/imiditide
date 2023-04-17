#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Sat Mar 20 14:22:39 2021

@author: trucdogunther

This is the fifth version of the RiPP-scanner program.
Like v4, this program finds ORFs nearby the query protein of interest and uses a reformatted 
blastp file from the "blastp-reformatter.py" program which converts protein hits with multiple 
scinames and taxids into individual protein hits. This version uses taxid instead of sciname to 
look up genome files. v4 uses pairwise alignment to find the single POI-associated ORF most 
likely to be correct. Pairwise alignment compares against an experimentally-determined 
RiPP peptide seqence.
Version v5 has additional functionality, recording the amino acid sequence of the predicted
POI (i.e., in this case methyltransferase protein sequence). Start and stop position of the
query-subject alignment in original blastp file is also recorded.

"""

#%% Import required modules and files

# input system module to obtain command-line argument from user
# arguments: upstream_window; downstream_window; min_ORF_len; max_ORF_len; working directory
# arguments: blastp_fiflename; control_RiPP_sequence
import sys

# allows user to interact with OS-dependent functionality
import os

# to check if path and file exists
import os.path

# import SeqIO subpackage from Bio package
from Bio import SeqIO

# import math module for rounding
import math

# import pairwise2 module and format_alignment method for pairwise sequence alignment
from Bio import pairwise2

# system arguments are all str
# RiPP-scanner.py upstream_window downstream_window min_ORF_len max_ORF_len working_directory blastp_file_name_ext control_RiPP_sequence

# user indicates length of region (in nucleotides) to scan upstream and downstream of POI for new ORFs
upstream_window = int(sys.argv[1])
downstream_window = int(sys.argv[2])

# user indicates the minimum and maximum length (in amino acids) for identified ORFs
min_ORF_len = int(sys.argv[3])
max_ORF_len = int(sys.argv[4])

# change the current working directory
os.chdir(sys.argv[5])

# user indicates path to blastp results file
blastp_filepath = sys.argv[5]

# indicate path to gbk files
gbk_filepath = sys.argv[5]

# user indicates name and extension (file.txt) of blastp results file
blastp_filename = sys.argv[6]

# user indicates peptide sequence of experimentally-determined RiPP ORF for pairwise alignment control
control_RiPP_sequence = sys.argv[7]

print("***************************************************************************************")
print("Note that start/end nucleotide position printed is offset by -1 to SnapGene coordinates")
print("***************************************************************************************")

#%% Class definitions

class blastp_record(object):
# define blastp_record class: each instance is one hit protein blastp result
# arguments: one blastp hit for one protein where columns are data fields in -outfmt 7
# all attributes are data-type str

    # constructor retrieves Blastp results and defines attributes
    def __init__(self, blastp_record):
        
        blastp_record_tabsplit = blastp_record.split("\t")
        
        self.qacc = blastp_record_tabsplit[0]
        self.saccver = blastp_record_tabsplit[1]
        self.stitle = blastp_record_tabsplit[2]
        self.evalue = blastp_record_tabsplit[3]
        self.bitscore = blastp_record_tabsplit[4]
        self.pident = blastp_record_tabsplit[5]
        self.ppos = blastp_record_tabsplit[6]
        self.qcov = blastp_record_tabsplit[7]
        self.qcovhsp = blastp_record_tabsplit[8]
        self.align_len = blastp_record_tabsplit[9]
        self.gap_opens = blastp_record_tabsplit[10]
        self.num_gaps = blastp_record_tabsplit[11]
        self.q_start = blastp_record_tabsplit[12]
        self.q_end = blastp_record_tabsplit[13]
        self.s_start = blastp_record_tabsplit[14]
        self.s_end = blastp_record_tabsplit[15]
        self.strand = blastp_record_tabsplit[16]
        self.taxid = blastp_record_tabsplit[17]
        self.sciname = blastp_record_tabsplit[18]
        self.comname = blastp_record_tabsplit[19]
        self.blastname = blastp_record_tabsplit[20]
   
    
class CDS_record(object):
# define CDS_record class: each instance is a SeqRecord feature match to one hit protein blastp result
# required argument: Bool flag indicating if CDS for hit protein is found in Genbank file
# optional arguments: one SeqRecord object and one SeqFeature object
    
    # constructor retrieves SeqRecord feature and defines attributes
    def __init__(self, CDS_flag, seq_record = None, seq_feature = None):
       
        
        # when a CDS is found in Genbank file with protein ID matching hit protein RefSeq ID
        if CDS_flag == True:
            
            # flag to indicate CDS match to blastp hit
            self.CDS_match = True
        
            # DNA sequence of contig/chromosome containing this POI        
            self.chr_seq = seq_record.seq
            
            # all information for this hit seq_feature
            self.seq_feature = seq_feature
            
            # corresponding DNA sequence coding this POI
            self.CDS_seq = seq_feature.extract(seq_record.seq)
            
            # sequence of query protein of interest
            self.POI_seq = seq_feature.qualifiers["translation"][0]
            
            # extract start and end position of the CDS encoding POI
            # use to determine default orientation of CDS (top or bottom strand) in original Genbank file
            # start and end positions are relative to 5' end of top strand in original Genbank file
            # but must account for whether top or negative strand
            
            # for top strand
            if seq_feature.strand == 1:
                self.CDS_start = seq_feature.location.start
                self.CDS_end = seq_feature.location.end
            # for bottom strand
            else:
                self.CDS_end = seq_feature.location.start
                self.CDS_start = seq_feature.location.end
            
            if self.CDS_start < self.CDS_end:
                self.CDS_strand = "top"
            else:
                self.CDS_strand = "bottom"
                
        # when a CDS is not found in Genbank file with protein ID matching hit protein RefSeq ID
        # possibly due to, e.g., misannotation in Genbank file vs. blastp database
        else:
            
            self.CDS_match = False
            self.chr_seq = "n/a"
            self.seq_feature = "no CDS matching blastp hit found in gbk file"
            self.CDS_seq = "n/a"
            self.CDS_strand = "n/a"
            self.POI_seq = "n/a"
                       
          
class ORF_record(object):
# define ORF_record class: each instance is an ORF with start codon to stop codon found in ROI
# arguments: ORF sequence as Bio.Seq.Seq obj, flag for plus/minus strand, start pos, end pos
# arguments: DNA sequence of ORF as Bio.Seq.Seq obj
    
    # constructor retrieves arguments and defines attribtes
    
    def __init__(self, ORF, strand_flag, ORF_start, ORF_end, ORF_DNA):
        
        self.ORF_seq = ORF
        self.ORF_strand = strand_flag
        self.ORF_start = ORF_start
        self.ORF_end = ORF_end
        self.ORF_len = len(ORF)
        self.ORF_DNA = ORF_DNA
    
#%% Function definitions
    
def gene_locator(POI_ID, gbk_filename, gbk_filepath):
# define a function that takes the RefSeq ID of a protein of interest and returns corresponding CDS information
# arguments: RefSeq ID of protein of interest (str), name and path of gbk file for organism (str)
# return: CDS_record object containing SeqRecord and SeqFeature information for POI

    # navigate to folder containing genbank file for subject organism
    os.chdir(gbk_filepath)
    
    # open Genbank file using organism name of current protein hit
    input_handle = open(gbk_filename, "r")
    
    # set flag to not exit loop unless POI record identified
    exit_flag = False
    
    # step through a SeqRecord iterator
    for seq_record in SeqIO.parse(input_handle, "genbank"):
        
        #print("*****************")
        #print("Dealing with Genbank record %s" % seq_record.name)
        
        # step through each feature within the current SeqRecord (i.e. chromosome)
        for seq_feature in seq_record.features:
            
            # search only coding sequences for matches to protein of interest
            if seq_feature.type == "CDS":
                
                # I added line 38 to handle Genbank files containing CDS with missing protein ID
                # ...which throws a KeyError because there is no "protein_id" key in the qualifiers dict within the seq_feature...
                # ...object. So CDSs with missing protein_id are excluded from the search
                
                if seq_feature.qualifiers.get("protein_id", "no protein_id") != "no protein_id":
                    
                    # extract and store RefSeq ID for current protein
                    protein_ID = seq_feature.qualifiers["protein_id"][0]
                    
                    # check if current protein ID matches RefSeq ID of protein of interest
                    if protein_ID == POI_ID:
                        
                        # extract and store relevant information of seq feature containing this POI
                        POI_record = CDS_record(True, seq_record, seq_feature)

                        # print the POI seq_feature information to check
                        # note that .location attribute values can be used directly to splice out POI DNA...
                        # ...sequence from POI_chromosome (i.e. same indexing scheme)
                        #print(POI_record.seq_feature)   
                        #print(POI_record.seq_feature.qualifiers["protein_id"][0])  
                        #print(POI_record.seq_feature.location)
                        #print(POI_record.seq_feature.strand)
                        #print(POI_record.CDS_seq)
                        #print("*****************")
                                        
                        # flag to exit out of seq_record for loop once hit protein reccord identified
                        # and exit out of seq_feature for loop
                        exit_flag = True
                        break
            
        if exit_flag == True:
            break

    # if no CDS protein ID matches RefSeq ID of protein of interest
    if exit_flag == False:
        POI_record = CDS_record(False)
         
    # close Genbank file             
    input_handle.close()
    
    # return dict storing genome information for POI
    return POI_record


def ORF_locator(CDS_obj, upstream_window, downstream_window, min_ORF_len, max_ORF_len):
# define a function to scan a genomic region of interest for ORFs on rev/for strands
# remember to import math module for mathematical manipulation
# arguments: CDS_record object corresponding to CDS of POI
# arguments: SeqFeatre.ExactPosition variables indicating CDS start and end position for POI
# arguments: ints indicating length of upstream and downstream window to POI
# arguments: ints indicating minimum and maximum length to filter ORFs of desired length
# return: list of ORF_record object with length, strand, location, sequence for each ORF
    
    # extract start and end position of the CDS encoding POI
    # position references 5' end of top strand
    CDS_start = CDS_obj.CDS_start
    CDS_end = CDS_obj.CDS_end
    
    # extract default orientation of CDS (top or bottom strand) in original Genbank file
    if CDS_obj.CDS_strand == "top":
        CDS_strand = "top"
    else:
        CDS_strand = "bottom"
        
    # use CDS start and end position and user-defined upstream/downstream windows
    # extract out entire genomic sequence encompassing upstream/POI/downstream ROIs
    # account for orientation of CDS_strand
    if CDS_strand == "top":
        genomic_ROI = CDS_obj.chr_seq[CDS_start-upstream_window:CDS_end+downstream_window]
    else:
        genomic_ROI = CDS_obj.chr_seq[CDS_end-downstream_window:CDS_start+upstream_window]
    
    # sequence of top strand is by default input genomic ROI (5' to 3') from Genbank file
    top_seq = genomic_ROI
    
    # sequence of bottom strand is by default reverse complement of genomic ROI (5' to 3')
    bottom_seq = genomic_ROI.reverse_complement()
    
    # store sequence of each strand in list to facilitate iteration
    strand_seq = [top_seq, bottom_seq]
    
    # define list to store identified ORFs of any length from both strands
    prelim_ORF_list = []
    
    # define list to store identified ORFs meeting certain criteria from both strands
    final_ORF_list = []
    
    # define list to store filtered ORFs for those closest to the POI
    final_ORF_list_2 = []
    
    # flag to start with top strand first
    strand_flag = "top"
    
    # iterate through each strand
    for strand in strand_seq:
                
        # iterate through each frame of current strand
        # each strand has 3 reading frames
        for frame in range(3):
    
            # calculate the number of codons (multiple of 3) for current frame     
            no_codons = math.floor((len(strand) - frame) / 3)
            
            # extract DNA sequence of current reading frame
            frame_seq = strand[frame : frame + no_codons*3]
            
            # translate entire extracted DNA sequence
            # use NCBI translation table = 11 for prokaryotes
            translated_seq = frame_seq.translate(table=11)
                        
            # obtain ORFs by splitting translated_seq with "*" delimiter denoting stop codon
            # ORFs are Bio.Seq.Seq objects
            # note that ORFs may encompass more than just a coding sequence
            all_ORFs = translated_seq.split("*")
            
            # delete last ORF in list because it is truncated by upstream/downstream window
            del all_ORFs[-1]
            
            #print(all_ORFs)
                              
            # calculate positions (inclusive boundaries) of nucleotides coding for ORF
            # add 2 to account for stop codon
            # with respect to 5' start of frame
            # store in dict that pairs ORF with its nucleotide sequence and nt start position
            # key is string; value is tuple of nucleotide sequence (Bio.seq.seq obj) and start index (int)
            ORF_nt_dict = {}
            
            # reset starting nucleotide index of ORF to find index of each ORF as interating
            ORF_nt_start = 0
            
            for ORF in all_ORFs:
    
                ORF_nt_end = ORF_nt_start + (len(ORF) * 3) + 2 
                ORF_nt_seq = frame_seq[ORF_nt_start : ORF_nt_end + 1]
                ORF_nt_dict[str(ORF)] = (ORF_nt_seq, ORF_nt_start)
                ORF_nt_start = ORF_nt_end + 1
                
            # iterate through identified ORFs and extract out ORFs from Met1 to stop codon   
            for ORF in all_ORFs:
                                
                # reset amino acid position to index 0 to scan from beginning of sequence for each ORF
                amino_acid_pos = 0
                
                # reset codon flag to False for each ORF to find Met1 and other alternative start codons
                codon_triplet_flag = False
                
                # find the first methionine in ORF
                for amino_acid in ORF:
                                
                    # if current amino acid is the canonical Met1 or less common alternative start residues
                    if amino_acid == "L" or amino_acid == "I" or amino_acid == "V" or amino_acid == "M":
                        
                        # trim ORF with potential start codon
                        trimmed_ORF = ORF[amino_acid_pos:]
                                                                        
                        # use dict pairing ORF with its DNA sequence to extract corrresponding trimmed DNA sequence
                        trimmed_nt = ORF_nt_dict[str(ORF)][0][amino_acid_pos * 3 : ]
                        
                        # calculate starting codon from extracted DNA sequence
                        codon_triplet = trimmed_nt[0:3]
                                                
                        # extract starting nucleotide index of trimmed ORF accounting for frame
                        # note that nt position here is relative to start of current frame for genomic ROI                    
                        ORF_nt_start_rel_to_ROI = ORF_nt_dict[str(ORF)][1] + (amino_acid_pos * 3)
                        
                        # flag that initiation codon found if Met1
                        if amino_acid == "M":
                            
                            codon_triplet_flag = True
                            
                        # check if alternative start sites have appropriate codons
                        elif amino_acid == "L":
                            
                                if codon_triplet == "TTG" or codon_triplet == "CTG":
                                
                                    # flag that initiation codon is found
                                    codon_triplet_flag = True
                                
                        elif amino_acid == "I":
                                
                                if codon_triplet == "ATT" or codon_triplet == "ATC" or codon_triplet == "ATA":
                                    
                                    codon_triplet_flag = True
                                    
                        else:
                            
                            if codon_triplet == "GTG":
                                
                                codon_triplet_flag = True
                    
                    # if proper initiation codon is identified above
                    if codon_triplet_flag == True:
                        
                        # back-calculate trimmed ORF (tORF) position relative to chromosome index 0 using CDS POI
                        # account for orientation of POI CDS to calculate 
                        # account for orientation of current frame
                        # remember that ORF_nt_start_rel_to_ROI is distance from 5' end of ROI strand start
                        # end position must account for stop codon (add 3)
                        
                        # for frames facing in the same orientation as the top strand
                        if strand_flag == "top":
                              
                            if CDS_strand == "top":
                                tORF_nt_start = CDS_start - upstream_window + ORF_nt_start_rel_to_ROI
                                
                            else:
                                tORF_nt_start = CDS_end - downstream_window + ORF_nt_start_rel_to_ROI
                            
                            # add 2 positions for stop codon in top strand (inclusive bounds)
                            tORF_nt_end = tORF_nt_start + (len(trimmed_ORF) * 3) + 2
                                                                                                            
                        # for frames facing in the same orientation as the bottom strand
                        # subtract 2 from start position so index scheme matches top strand
                        else:
                            
                            if CDS_strand == "top":
                                tORF_nt_start = CDS_end + downstream_window - ORF_nt_start_rel_to_ROI - 2
                            else:
                                tORF_nt_start = CDS_start + upstream_window - ORF_nt_start_rel_to_ROI - 2
                            
                            # subtract 2 positions for stop codon in bottom strand (inclusive bounds)
                            tORF_nt_end = tORF_nt_start - (len(trimmed_ORF) * 3) - 2
                            
                        # store detected ORF                                                                  
                        ORF_obj = ORF_record(trimmed_ORF, strand_flag, tORF_nt_start, tORF_nt_end, trimmed_nt)
                        prelim_ORF_list.append(ORF_obj)
                        
                        #print("Full ORF: " + str(ORF))
                        #print("Strand direction: " + ORF_obj.ORF_strand)
                        #print("Trimmed ORF start: " + str(ORF_obj.ORF_start))
                        #print("Trimmed ORF end: " + str(ORF_obj.ORF_end))
                        #print("Trimmed ORF: " + str(ORF_obj.ORF_seq))
                        #print("DNA seq of trimmed ORF: " + str(ORF_obj.ORF_DNA))
                        #print("-----------")
                              
                        break
              
                    amino_acid_pos = amino_acid_pos + 1
            
        strand_flag = "bottom"
        
    # filter out ORFs not satisfying desired conditions                 
    for ORF_obj in prelim_ORF_list:
                   
        # must meet user-defined length requirements
        if ORF_obj.ORF_len <= max_ORF_len and ORF_obj.ORF_len >= min_ORF_len:
                        
            # must have Asp in the amino acid sequence
            if "D" in ORF_obj.ORF_seq:
                              
                # must be in the same direction as POI gene
                if ORF_obj.ORF_strand == CDS_strand:
                                
                    # ORF start codon must not be within methyltransferase POI sequence itself
                    # position listed in Genbank relative to index 0 site of plus strand
                    # subtract coordinates of minus strand ORFs from length of chromosome
                    if ORF_obj.ORF_strand == "top":
                        
                        if ORF_obj.ORF_start < CDS_start or ORF_obj.ORF_end > CDS_end:
                            final_ORF_list.append(ORF_obj)
                            
                    else:
                        
                        if ORF_obj.ORF_start > CDS_start or ORF_obj.ORF_end < CDS_end:
                            final_ORF_list.append(ORF_obj)
                            
    # only keep ORFs that lie closest to the original protein of interest
    for ORF_obj in final_ORF_list:
        
        # # for ORFs found entirely within POI CDS
        # if ORF_obj.ORF_strand == "top":
            
        #     if ORF_obj.ORF_start >= CDS_start and ORF_obj.ORF_end <= CDS_end:
        #         final_ORF_list_2.append(ORF_obj)
        #         continue
                
        # else:
            
        #     if ORF_obj.ORF_start <= CDS_start and ORF_obj.ORF_end >= CDS_end:
        #         final_ORF_list_2.append(ORF_obj)
        #         continue
                
        # for ORFs whose start position are within the POI CDS but end position are outside POI CDS
        if ORF_obj.ORF_strand == "top":
            
            if ORF_obj.ORF_start > CDS_start and ORF_obj.ORF_start < CDS_end and ORF_obj.ORF_end > CDS_end:
                # only retain ORFs whose overlap with POI CDS are <= 25 nucleotides
                if (CDS_end - ORF_obj.ORF_start) <= 25:
                    final_ORF_list_2.append(ORF_obj)
                    continue
                else:
                    continue
                
        else:
            
            if ORF_obj.ORF_start < CDS_start and ORF_obj.ORF_start > CDS_end and ORF_obj.ORF_end < CDS_end:
                # only retain ORFs whose overlap with POI CDS are <= 25 nucleotides
                if (ORF_obj.ORF_start - CDS_end) <= 25:
                    final_ORF_list_2.append(ORF_obj)
                    continue
                else:
                    continue
                
        # for ORFs whose end position are within the POI CDS but start position are outside POI CDS
        if ORF_obj.ORF_strand == "top":
            if ORF_obj.ORF_end < CDS_end and ORF_obj.ORF_end > CDS_start and ORF_obj.ORF_start < CDS_start:
                # only retain ORFs whose overlap with POI CDS are <= 25 nucleotides
                if (ORF_obj.ORF_end - CDS_start) <= 25:
                    final_ORF_list_2.append(ORF_obj)
                    continue
                else:
                    continue
                
        else:    
            if ORF_obj.ORF_end > CDS_end and ORF_obj.ORF_end < CDS_start and ORF_obj.ORF_start > CDS_start:
                # only retain ORFs whose overlap with POI CDS are <= 25 nucleotides
                if (CDS_start - ORF_obj.ORF_end) <= 25:
                    final_ORF_list_2.append(ORF_obj)
                    continue
                else:
                    continue
                
        # for ORFs found in their entirety outside but close to POI CDS  
        if ORF_obj.ORF_strand == "top":
            
            # for ORFs found outside and upstream of POI CDS within 50 bp window of POI CDS start
            if ORF_obj.ORF_start < CDS_start and ORF_obj.ORF_end < CDS_start:
                if (CDS_start - ORF_obj.ORF_end) <= 50:
                    final_ORF_list_2.append(ORF_obj)
                    continue
                else:
                    continue
                
            # for ORFs found outside and downstream of POI CDS within 50 bp window of POI CDS end
            if ORF_obj.ORF_start < CDS_end and ORF_obj.ORF_end < CDS_end:
                if (ORF_obj.ORF_start - CDS_end) <= 50:
                    final_ORF_list_2.append(ORF_obj)
                else:
                    continue
        
        else:
            # for ORFs found outside and upstream of POI CDS within 50 bp window of POI CDS start
            if ORF_obj.ORF_start > CDS_start and ORF_obj.ORF_end > CDS_start:
                if (ORF_obj.ORF_end - CDS_start) <= 50:
                    final_ORF_list_2.append(ORF_obj)
                    continue
                else:
                    continue
                
            # for ORFs found outside and downstream of POI CDSs within 50 bp window of POI CDS end
            if ORF_obj.ORF_start < CDS_end and ORF_obj.ORF_end < CDS_end:
                if (CDS_end - ORF_obj.ORF_start) <= 50:
                    final_ORF_list_2.append(ORF_obj)   
                    continue
                else:
                    continue
                          
    # return list of ORF_record object
    return final_ORF_list_2             


def ORF_pairwise_aligner(dict_ORF_list, control_peptide):
# define a function to iterate through a list of predicted ORFs for each POI
# and extract the one ORF with the best pairwise alignment score to a control peptide sequence
# remember to import Bio.pairwise2 module for this function to work
# arguments: dict of list of ORF_record object with length, strand, location, sequence for each ORF
# input dict key is POI ID and value is list of ORF_record objects
# arguments continued: str of control amino acid sequence for pairwise alignment
# return: dict of single ORF_record object corresponding to top ORF for each POI

    # dict-type variable to store processed CDS objects filtered for the closest ORF to POI
    # POI RefSeq ID is key and values are the closest ORF (ORF_record obj) identified for corresponding POI
    ORF_dict = {}
    
    for POI_ID, ORF_list in dict_ORF_list.items():
            
        # initialize dict variable to store best alignment for each ORF in list identified for each POI
        ORF_alignments = {}
        
        for ORF in ORF_list:
            
            # determine all pairwise alignments between control amino acid sequence and current ORF amino acid sequence
            alignments = pairwise2.align.globalxx(control_peptide, ORF.ORF_seq)
            
            # iterate through alignments for current ORF and keep only the alignment with the highest score
            # store the ORF_record object corresponding to the top alignment so far
            # initialize best_alignment to first alignment of current ORF
            best_alignment = alignments[0]
            
            for alignment in alignments:
                
                if alignment.score > best_alignment.score:
                    best_alignment = alignment
                    
            # store tuple of best pairwise alignment and ORF_record obj for current ORF of current POI
            ORF_alignments[ORF.ORF_seq] = (best_alignment, ORF)
            
        # now iterate through each ORF of current POI and extract ORF with highest pairwise alignment score
        # store ORF_record object corresponding to the best ORF in dict with POI ID as key
        # initialize int variable to find ORF with best alignment score
        best_alignment_score = 0
        
        for ORF_aa_seq, alignment_ORFobj_tuple in ORF_alignments.items():
            
            if alignment_ORFobj_tuple[0].score > best_alignment_score:
                best_alignment_score = alignment_ORFobj_tuple[0].score
                ORF_dict[POI_ID] = alignment_ORFobj_tuple[1]
                                
    # return dict of ORF_record objects
    return ORF_dict

                
#%% To search for proteins of interest in the genome, first process blastp file to extract RefSeq IDs

# blastp files are formatted with -outfmt 7 from BLAST+ package
# header and descriptions from original blastp file were removed in reformatting process
    
# navigate to folder containing blastp results file
os.chdir(blastp_filepath)
    
# open blastp file to read and then close
blastp_file_object = open(blastp_filename)
blastp_file = blastp_file_object.read()
blastp_file_object.close()

# split blastp_file with newline character, and blastp results start at line 1 in file (index 0 in list)
blastp_list =  blastp_file.rstrip("\n").split("\n")

# dict-type variable to store processed blast objects
blastp_dict = {}
  
# iterate through list of blastp results for each protein hit
# convert blastp results of each protein hit to a blastp_record object with fields as attributes
# store the blastp_record objects as values in a dict with corresponding RefSeq ID (saccver) as key
for protein_hit in blastp_list:
    temp_blastp_obj = blastp_record(protein_hit)
    blastp_dict[temp_blastp_obj.saccver] = temp_blastp_obj

#%% Locate and extract DNA sequence and region encoding protein of interest from RefSeq WGS assembly

# iterate through each protein blastp hit to retrieve its DNA sequence and CDS location
# store returned CDS_record objects as values in a dict with corresponding RefSeq ID (saccver) as key
  
# dict-type variable to store processed CDS objects
CDS_dict = {}
     
for POI_ID, blastp_obj in blastp_dict.items():
        
    # indicate Genbank file name using organism taxid current protein hit
    gbk_filename = blastp_obj.taxid
    
    # search for protein of interest in the genome with its RefSeq ID
    # retrieve corresponding CDS information and store as CDS_record object
    # only search for CDS if there is a Genbank file for the organism

    if os.path.isfile(gbk_filepath + "/" + gbk_filename + ".gbff") == True:
        CDS_dict[POI_ID] = gene_locator(POI_ID, gbk_filename + ".gbff", gbk_filepath)
    elif os.path.isfile(gbk_filepath + "/" + gbk_filename + "-0.gbff") == True:
        CDS_dict[POI_ID] = gene_locator(POI_ID, gbk_filename + "-0.gbff", gbk_filepath)
    elif os.path.isfile(gbk_filepath + "/taxid" + gbk_filename + ".gbff") == True:
        CDS_dict[POI_ID] = gene_locator(POI_ID, "taxid" + gbk_filename + ".gbff", gbk_filepath)
    else:
        print(gbk_filename)

    
#%% Delete unused variables up to this point
    
del blastp_file
del blastp_file_object
del blastp_list
del blastp_obj
del protein_hit
del temp_blastp_obj

#%% Extract DNA sequences upstream and downstream of CDS coding for protein of interest

# iterate through each CDS entry to retrieve sequences upstream/downstream of POI using user-indicated windows
# store upstream/downstream sequences in a dict with corresponding RefSeq ID (saccver) as key
# dict will also store information of any identified upstream and downstream ORFs
  
# dict-type variable to store processed CDS objects
ORF_dict_full = {}
     
for POI_ID, CDS_obj in CDS_dict.items():
        
    # scan genomic region of interest for ORFs in both top and bottom strands
    # record in dict with POI RefSeq ID as key
    # values are lists; items in list are all the ORFs identified for corresponding POI
    # each ORF is stored as ORF_record obj with length, strand, location, sequence info, corresponding DNA seq
    # function also retrieves only ORFs having desired filters
    # only search for ORFs if CDS was found in original Genbank file with POI RefSeq ID
    if CDS_obj.CDS_match == True:
        ORF_dict_full[POI_ID] = ORF_locator(CDS_obj, upstream_window, downstream_window, min_ORF_len, max_ORF_len)
        
# iterate through list of identified ORFs for each POI and extract the single ORF closest to
# and/or overlapping with the POI to remove false positives and obtain the best ORF hit

# dict-type variable to store processed CDS objects filtered for the POI-associated ORF with highest
# pairwise alignment score to experimentally-determined RiPP peptide sequence
# POI RefSeq ID is key and values are the best ORF (ORF_record obj) identified for corresponding POI
ORF_dict = ORF_pairwise_aligner(ORF_dict_full, control_RiPP_sequence)
        

#%% Print out results to .txt file

os.chdir(blastp_filepath)

out_file_name = "RiPP-scanner-ORFs-test-2.txt"
out_file = open(out_file_name, "w")

out_file.write("Upstream window (in nt): " + str(upstream_window) + "\n")
out_file.write("Downstream window (in nt): " + str(downstream_window)  + "\n")
out_file.write("Minimum ORF length (in amino acids): " + str(min_ORF_len)  + "\n")
out_file.write("Maximum ORF length (in amino acids): " + str(max_ORF_len)  + "\n")

out_file.write("Query protein RefSeq ID\t")
out_file.write("Query organism\t")
out_file.write("Query protein sequence\t")
out_file.write("Query protein strand\t")
out_file.write("Query protein start position\t")
out_file.write("Query protein end position\t")
out_file.write("ORF sequence\t")
out_file.write("ORF strand\t")
out_file.write("ORF start position\t")
out_file.write("ORF end position\t")
out_file.write("ORF length\t")
out_file.write("ORF DNA sequence\t")
out_file.write("Query protein MT start pos\t")
out_file.write("Query protein MT end pos\n")

for POI_ID, ORF in ORF_dict.items():
    
    q_organism = blastp_dict[POI_ID].sciname
    q_strand = CDS_dict[POI_ID].CDS_strand
    q_start = CDS_dict[POI_ID].CDS_start
    q_end = CDS_dict[POI_ID].CDS_end
    q_POIseq = CDS_dict[POI_ID].POI_seq
        
    out_file.write(POI_ID + "\t")
    out_file.write(q_organism + "\t")
    out_file.write(q_POIseq + "\t")
    out_file.write(q_strand + "\t")
    out_file.write(str(q_start) + "\t")
    out_file.write(str(q_end) + "\t")
    out_file.write(str(ORF.ORF_seq) + "\t")
    out_file.write(str(ORF.ORF_strand) + "\t")
    out_file.write(str(ORF.ORF_start) + "\t")
    out_file.write(str(ORF.ORF_end) + "\t")    
    out_file.write(str(ORF.ORF_len) + "\t")  
    out_file.write(str(ORF.ORF_DNA) + "\t")  
    out_file.write(str(blastp_dict[POI_ID].s_start) + "\t")
    out_file.write(str(blastp_dict[POI_ID].s_end) + "\n")

out_file.close()
 
#%% Delete unused variables up to this point
    
del CDS_obj
del POI_ID
del q_organism
del q_strand
del q_start
del q_end
del q_POIseq