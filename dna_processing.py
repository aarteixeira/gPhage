"""
    This script is part of Teixeira et al 2020
    Please cite: Teixeira AAR, et al. A refined genome phage display methodology delineates the human antibody response in patients with Chagas disease. iScience 2021; 24:102540. 
    
    Created by Andre Teixeira and Ricardo Giordano

    This program is free software: you can redistribute it and/or modify
    it under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
    GNU General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program.  If not, see <https://www.gnu.org/licenses/>.
"""

"""
Process NGS sequences from T. cruzi lib built on pG8SAET
- Reads multiple FastQ files (paired-assembled) - see line 46 to 53
- Extract insert sequence
- Blast sequence against genome
- Tranlate DNA to peptide sequences
- Blast Peptide sequences
- calculate hydropathy
- output: excel files containing dna and peptide information

3rd party software needed:
- Blast
"""
import os
import pandas as pd
import numpy as np
import time
from collections import Counter, OrderedDict
from Bio import SeqIO
from Bio.Blast import NCBIXML

##############################################################################
### Options                                                                ###
##############################################################################
BLASTN_EXE = 'ncbi-blast-2.12.0+/bin/blastn'
BLASTN_DB = 'example/ncbi_genomes.fasta'
BLASTP_EXE = 'ncbi-blast-2.12.0+/bin/blastp'
BLASTP_DB = 'example/all_proteins.fasta'

# Read sequences from FASTq, trim vector seqs and count
option_read_from_fastq = True
dict_sample_fastq = OrderedDict()
dict_sample_fastq['Sample1'] = 'example/example.fastq'

flank_upstream   = 'ATGACCATGGCAGTAC'  # vector flanking sequences
flank_downstream = 'GTACCCGGTGCGCCGG'

# remove singletons?
option_remove_singletons = True

# BLAST DNA vs. genome/database?
option_blast_dna = True

# exclude reads according to DNA identity cutoff?
option_remove_DNA_below_cutoff = True
cutoff_dna_identity = 0.9

# translate sequences?
option_translate_sequences = True

# Group sequences by peptide?
option_group_peptide_table = True

# hydrophilicity
option_hydrophilicity = True

# BLAST peptides vs. proteome/database?
option_blast_peptides = True

# exclude peptides according to similarity cutoff?
option_remove_peptide_below_cutoff = True
cutoff_peptide_similarity = 0.6

# get clustering data and create tables
option_process_clustering = False

# blastp stats
option_blastp_stats = True
blastp_cutoff = 0.6

##############################################################################

def hamming(s1, s2):
    return sum(ch1 != ch2 for ch1, ch2 in zip(s1, s2))

def translate(strand, frame=1):
    if frame == 2:
        strand = strand[1:]
    elif frame == 3:
        strand = strand[2:]
    codons = {
        'GCA': 'A', 'GCC': 'A', 'GCG': 'A', 'GCT': 'A', 'TGC': 'C',
        'TGT': 'C', 'GAC': 'D', 'GAT': 'D', 'GAA': 'E', 'GAG': 'E',
        'TTC': 'F', 'TTT': 'F', 'GGC': 'G', 'GGT': 'G', 'GGG': 'G',
        'GGA': 'G', 'CAC': 'H', 'CAT': 'H', 'ATA': 'I', 'ATC': 'I',
        'ATT': 'I', 'AAA': 'K', 'AAG': 'K', 'CTA': 'L', 'CTC': 'L',
        'CTG': 'L', 'CTT': 'L', 'TTA': 'L', 'TTG': 'L', 'ATG': 'M',
        'AAC': 'N', 'AAT': 'N', 'CCA': 'P', 'CCC': 'P', 'CCG': 'P',
        'CCT': 'P', 'CAA': 'Q', 'CAG': 'Q', 'AGA': 'R', 'AGG': 'R',
        'CGA': 'R', 'CGC': 'R', 'CGG': 'R', 'CGT': 'R', 'AGC': 'S',
        'AGT': 'S', 'TCA': 'S', 'TCC': 'S', 'TCG': 'S', 'TCT': 'S',
        'TAA': '_', 'TAG': 'Q', 'TGA': '_', 'ACA': 'T',
        'ACC': 'T', 'ACG': 'T', 'ACT': 'T', 'GTA': 'V', 'GTC': 'V',
        'GTG': 'V', 'GTT': 'V', 'TGG': 'W', 'TAC': 'Y', 'TAT': 'Y'}
    lenght = len(strand) - len(strand) % 3
    return ''.join(codons[strand[i:i + 3]] for i in range(0, lenght, 3))

def trim_generator(seq_iter, flank5, flank3):
    #
    # trim sequences based on up and downstream flanking sequences
    # return list of trimmed seqs, # processed seqs, # trimmed seqs 
    #
    n_reads = 0
    seq_list = []
    seq_list_append = seq_list.append
    for seq in seq_iter:
        try:
            a = str(seq.seq).split(flank5)[1].split(flank3)
            a[1]  # test if splitted successfully
            seq = a[0]
            seq_list_append(seq)
        except IndexError:
            pass

        n_reads += 1

    return seq_list, n_reads, len(seq_list)


def parker_hydropathy(pep):

    pep = pep.upper()

    energy = {
        'A': 2.1,
        'R': 4.2,
        'N': 7,
        'D': 10,
        'C': 1.4,
        'Q': 6,
        'E': 7.8,
        'G': 5.7,
        'H': 2.1,
        'I': -8,
        'L': -9.2,
        'K': 5.7,
        'M': -4.2,
        'F': -9.2,
        'P': 2.1,
        'S': 6.5,
        'T': 5.2,
        'W': -10,
        'Y': -1.9,
        'V': -3.7}

    return sum(energy[aa] for aa in pep)

#
# Reading sequences from FASTq
#
if option_read_from_fastq:

    '''
    Dataframe for statistics from the processing. It will be initiated
    only if reading from FASTq, otherwise I'll try to read from a excel file
    '''
    stats = pd.DataFrame()
    dna_df = pd.DataFrame()

    # go over all samples defined at the dictionary
    for sample in dict_sample_fastq.keys():
 
        # series to stores statistics from the sample processing
        sample_stats = pd.Series(name=sample)

        # open FASTQ and store all sequences in the FASTq in a generator
        print('Reading sample {} from {}'.format(sample,
                                                 dict_sample_fastq[sample]))
        itera_seqIO = SeqIO.parse(dict_sample_fastq[sample], 'fastq')

        # trim vector sequences
        trimmed_dna_seqs, n_reads, n_trimmed_reads = trim_generator(
                                                             itera_seqIO,
                                                             flank_upstream,
                                                             flank_downstream)

        # Count inserts and insert into DNA dataframe 
        insert_dna_count = pd.Series(Counter(trimmed_dna_seqs), name=sample)
        dna_df = dna_df.join(insert_dna_count, how='outer')
        dna_df = dna_df.fillna(value=0)

        # Output stats and store the counts
        print('# reads: {}'.format(n_reads))
        print('# reads trimmed: {}'.format(n_trimmed_reads))
        print('# different DNA {}'.format(len(insert_dna_count)))
        sample_stats['# reads'] = n_reads
        sample_stats['# trimmed reads'] = n_trimmed_reads
        sample_stats['# different DNA'] = len(insert_dna_count)
        stats = stats.join(sample_stats, how='outer')

    stats = stats.T
    # Write excel file with DNA dataframe
    print('Writing EXCEL dataframes...')
    stats.to_excel('statistics.xlsx')
    dna_df.to_excel('dna_dataframe.xlsx')

#
# Remove singletons
#
if option_remove_singletons:

    print('Removing singletons...')

    # check if DNA dataframe exists, otherwise read excel
    try:
        dna_df
        stats
    except NameError:
        try:
            print('Loading dataframe...')
            dna_df = pd.read_excel('dna_dataframe.xlsx')
            stats = pd.read_excel('statistics.xlsx')
            print('Done!')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')
            exit()

    # find sequence that are not singletons and keep them
    samples = [sample for sample in dict_sample_fastq]
    not_singleton = dna_df.loc[:, samples].max(axis=1) > 1
    dna_df = dna_df.loc[not_singleton, :]

    print('Singletons removed')

    # count unique and total reads no singletons in all samples
    for sample in samples:
        no_singleton = dna_df.loc[:, samples].max(axis=1) > 1
        no_singleton = dna_df.loc[no_singleton, sample] > 0

        n_unique_nosingleton = dna_df.loc[no_singleton, sample].count()
        n_reads_nosingleton = dna_df.loc[no_singleton, sample].sum()

        stats.loc[sample, '# reads no singleton'] = n_reads_nosingleton
        stats.loc[sample, '# unique no singleton'] = n_unique_nosingleton
        # stats.set_value(sample, '# reads no singleton',
        #                 value=n_reads_nosingleton)
        # stats.set_value(sample, '# unique no singleton',
        #                 value=n_unique_nosingleton)

    print(stats)
    # write to excel
    dna_df.to_excel('dna_dataframe.xlsx')
    stats.to_excel('statistics.xlsx')

    #print(not_singleton)

#
# BlastN 
#
if option_blast_dna:
    print('Blast DNA against genome')

    # check if DNA dataframe exists, otherwise read excel
    try:
        dna_df
        stats
    except NameError:
        try:
            dna_df = pd.read_excel('dna_dataframe.xlsx')
            stats = pd.read_excel('statistics.xlsx')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')
            exit()

    # create fasta
    print('Writing fasta')
    with open('dna_temp.fasta', 'w') as f:
        for seq in dna_df.index:
            if str(seq) != 'nan':
                f.write('>{0}\n{0}\n'.format(seq))

    # blastN
    print('blastN!!!')
    os.system('{exe} -query dna_temp.fasta -task blastn-short -db {db} -out blastn.xml -outfmt 5  -max_hsps 1 -max_target_seqs 1 -num_threads 8 2>/dev/null'.format(exe=BLASTN_EXE, db=BLASTN_DB))

    # parse Blast output
    print('Parsing blastN output')
    blast_records = NCBIXML.parse(open('blastn.xml'))
    for record in blast_records:

        # if there are hits, get only the first HSP
        if len(record.alignments) > 0:
            
            #print(record.__dict__)
            seq = record.query
            name = record.alignments[0].hit_def
            accession = record.alignments[0].accession

            #print(record.alignments[0].hsps[0].__dict__)
            n_identities = record.alignments[0].hsps[0].identities
            query_len = record.query_length
            identity_ratio = n_identities / float(query_len)
            evalue = record.alignments[0].hsps[0].expect
            
            # update dna dataframe with blast info
            dna_df.loc[seq, 'blastN hit'] = name
            dna_df.loc[seq, 'blastN identity ratio'] = identity_ratio
            dna_df.loc[seq, 'blastN expect'] = evalue

    # generate histograms of identity
    dna_df.loc[:, 'blastN identity ratio'].fillna(value=0, inplace=True)
    dna_df.to_excel('dna_dataframe.xlsx')


# 
# gen DNA blast stats, histograms, apply DNA cutoff
#
if option_remove_DNA_below_cutoff:

    print('remove DNA below cutoff!')

    try:
        dna_df
        stats
    except NameError:
        try:
            print('Loading DNA dataframe...')
            dna_df = pd.read_excel('dna_dataframe.xlsx')
            stats = pd.read_excel('statistics.xlsx')
            print('Done!')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')
            exit()

    dna_below_cutoff = dna_df.loc[:, 'blastN identity ratio'] < cutoff_dna_identity
    dna_below_cutoff_df = dna_df.loc[dna_below_cutoff, :]
    dna_below_cutoff_df.to_excel('dna_below_cutoff.xlsx')

    dna_above_cutoff = dna_df.loc[:, 'blastN identity ratio'] >= cutoff_dna_identity
    dna_above_cutoff_df = dna_df.loc[dna_above_cutoff, :]
    dna_above_cutoff_df.to_excel('dna_above_cutoff.xlsx')

    # count unique and total reads DNA above cutoff in all samples
    for sample in dict_sample_fastq:

        dna_in_sample = dna_above_cutoff_df.loc[:, sample] > 0
        n_unique_above_dna_cutoff = dna_above_cutoff_df.loc[dna_in_sample, sample].count()
        n_reads_above_dna_cutoff = dna_above_cutoff_df.loc[dna_in_sample, sample].sum()

        stats.loc[sample, '# reads above DNA cutoff'] = n_reads_above_dna_cutoff
        stats.loc[sample, '# unique above DNA cutoff'] = n_unique_above_dna_cutoff

    stats.to_excel('statistics.xlsx')
    dna_df = dna_above_cutoff_df

#
# TRANSLATE EVERYTHING
#
if option_translate_sequences:

    print('Translating peptides')
    try:
        dna_df
        stats
    except NameError:
        try:
            print('Loading dataframe...')
            dna_df = pd.read_excel('dna_above_cutoff.xlsx')
            stats = pd.read_excel('statistics.xlsx')
            print('done!')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')
            exit()


    # check number of nucleotides ( n%3 = 1)
    correct_frame = [a for a in dna_df.index if len(a) % 3 == 1]

    # translate the ones in the correct frame
    for seq in correct_frame:
        peptide = translate(seq, frame=2)
        
        # check for a stop codon
        if '_' not in peptide:
            # dna_df.set_value(seq, 'peptide', value=peptide)
            dna_df.loc[seq, 'peptide'] = peptide

    dna_df.to_excel('dna_above_cutoff.xlsx')

    # create peptide df
    samples = [sample for sample in dict_sample_fastq]
    peptide_df = dna_df[samples].groupby(dna_df['peptide']).sum()
    #peptide_df.drop('peptide')

    # drop the ones without peptide
    #peptides = [peptide for peptide in peptide_df.index if peptide != '']
    #peptide_df = peptide_df.loc['']

    # write dataframe
    print('Writing peptide dataframe')
    peptide_df.to_excel('peptide_dataframe.xlsx')

    # calculate unique and reads with peptides
    print('calculating peptite stats')
    for sample in dict_sample_fastq:
        # count only peptides that show up in sample
        peptides_in_sample = peptide_df.loc[:, sample] > 0
        n_unique_peptide = peptide_df.loc[peptides_in_sample, sample].count()
        n_reads_peptide = peptide_df.loc[peptides_in_sample, sample].sum()

        stats.loc[sample, '# reads with peptide'] = n_reads_peptide
        stats.loc[sample, '# different peptides'] = n_unique_peptide
        # stats.set_value(sample, '# reads with peptide',
        #                 value=n_reads_peptide)
        # stats.set_value(sample, '# different peptides',
        #                 value=n_unique_peptide)

    stats.to_excel('statistics.xlsx')

#
# calculate hydrophilicity
#
if option_hydrophilicity:
    print('calculate hydrophilicity')

    # check if DNA dataframe exists, otherwise read excel
    try:
        peptide_df
        stats
    except NameError:
        try:
            print('Loading peptide dataframe...')
            peptide_df = pd.read_excel('peptide_dataframe.xlsx', index_col=0)
            stats = pd.read_excel('statistics.xlsx')
            print('Done!')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')   
            exit()

    #print(peptide_df.head())
    print('calculating hydrophilicity...')
    for peptide in peptide_df.index:
        peptide_df.at[peptide, 'hydrophilicity'] = parker_hydropathy(peptide)
        
    # generate peptides.txt
    with open("peptides.txt", "w") as fout:
        for p in peptide_df.index:
            fout.write(f"{p}\n")

    print(peptide_df.head())
    peptide_df.to_excel('peptide_dataframe.xlsx')


#
# Blast peptides
#
if option_blast_peptides:
    print('Blast peptides against proteome')

    # check if DNA dataframe exists, otherwise read excel
    try:
        peptide_df
        stats
    except NameError:
        try:
            print('Loading peptide dataframe...')
            peptide_df = pd.read_excel('peptide_dataframe.xlsx', index_col=0)
            stats = pd.read_excel('statistics.xlsx')
            print('Done!')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')   
            exit()

    print(peptide_df.head())
    #exit()

    # create fasta
    print('Writing fasta')
    with open('peptide_temp.fasta', 'w') as f:
        for seq in peptide_df.index:
            if str(seq) != 'nan':
                f.write('>{0}\n{0}\n'.format(seq))

    # blastN
    print('blastp!!!')
    #os.system('{blastp} -outfmt 5 -evalue 100 -max_target_seqs 1 -max_hsps 1 -num_threads 8 -db {blastdb} -query peptide_temp.fasta > blastp.xml'.format(blastp=blastp_exe, blastdb=database))
    os.system('{blastp} -word_size 6 -outfmt 5 -gapopen 13 -evalue 100 -max_target_seqs 1 -max_hsps 1 -num_threads 10 -db {blastdb} -query peptide_temp.fasta > blastp.xml 2>/dev/null'.format(blastp=BLASTP_EXE, blastdb=BLASTP_DB))


    # parse Blast output
    print('Parsing blastP output')
    blast_records = NCBIXML.parse(open('blastp.xml'))
    for record in blast_records:

        # if there are hits, get only the first HSP
        if len(record.alignments) > 0:
            
            #print(record.__dict__)
            seq = record.query
            name = record.alignments[0].hit_def
            accession = record.alignments[0].accession

            #print(record.alignments[0].hsps[0].__dict__)
            n_identities = record.alignments[0].hsps[0].identities
            n_positives = record.alignments[0].hsps[0].positives
            query_len = record.query_length
            identity_ratio = n_identities / float(query_len)
            positive_ratio = n_positives / float(query_len)           
            evalue = record.alignments[0].hsps[0].expect
            name = record.alignments[0].hit_def
            geneId = name.split(' ')[0].replace('>', '')
            query = record.alignments[0].hsps[0].query
            match = record.alignments[0].hsps[0].sbjct                  

            # update dna dataframe with blast info
            '''peptide_df.set_value(seq, 'blastP hit', value=name)
            peptide_df.set_value(seq, 'blastP identity ratio', value=identity_ratio)
            peptide_df.set_value(seq, 'blastP positive ratio', value=positive_ratio)
            peptide_df.set_value(seq, 'blastP query length', value=record.query_length)
            peptide_df.set_value(seq, 'blastP expect', value=evalue)
            peptide_df.set_value(seq, 'blastP match name', value=name)
            peptide_df.set_value(seq, 'blastP gene id', value=geneId)
            peptide_df.set_value(seq, 'blastP query', value=query)
            peptide_df.set_value(seq, 'blastP match', value=match)'''

            peptide_df.at[seq, 'blastP hit'] = name
            peptide_df.at[seq, 'blastP identity ratio'] = identity_ratio
            peptide_df.at[seq, 'blastP positive ratio'] = positive_ratio
            peptide_df.at[seq, 'blastP query length'] = record.query_length
            peptide_df.at[seq, 'blastP expect'] = evalue
            peptide_df.at[seq, 'blastP match name'] = name
            peptide_df.at[seq, 'blastP gene id'] = geneId
            peptide_df.at[seq, 'blastP query'] = query
            peptide_df.at[seq, 'blastP match'] = match

    peptide_df.to_excel('peptide_dataframe.xlsx')

#
# Filter peptides below cutoff
#
if option_remove_peptide_below_cutoff:
    print('Remove peptides below cutoff')

    # check if peptide dataframe exists, otherwise read excel
    try:
        peptide_df
        stats
    except NameError:
        try:
            print('Loading peptide dataframe...')
            peptide_df = pd.read_excel('peptide_dataframe.xlsx')
            stats = pd.read_excel('statistics.xlsx')
            print('Done!')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')
            exit()

    peptide_below_cutoff = peptide_df.loc[:, 'blastP positive ratio'] < cutoff_peptide_similarity
    peptide_below_cutoff_df = peptide_df.loc[peptide_below_cutoff, :]
    peptide_below_cutoff_df.to_excel('peptide_below_cutoff.xlsx')

    peptide_above_cutoff = peptide_df.loc[:, 'blastP positive ratio'] >= cutoff_peptide_similarity
    peptide_above_cutoff_df = peptide_df.loc[peptide_above_cutoff, :]
    peptide_above_cutoff_df.to_excel('peptide_above_cutoff.xlsx')
    

    # count unique and total reads DNA above cutoff in all samples
    for sample in dict_sample_fastq:

        peptide_in_sample = peptide_above_cutoff_df.loc[:, sample] > 0
        n_unique_above_peptide_cutoff = peptide_above_cutoff_df.loc[peptide_in_sample, sample].count()
        n_reads_above_peptide_cutoff = peptide_above_cutoff_df.loc[peptide_in_sample, sample].sum()

        stats.loc[sample, '# reads above peptide blastP cutoff'] = n_reads_above_peptide_cutoff
        stats.loc[sample, '# different peptide blastP cutoff'] = n_unique_above_peptide_cutoff
        # stats.set_value(sample, '# reads above peptide blastP cutoff',
        #                 value=n_reads_above_peptide_cutoff)
        # stats.set_value(sample, '# different peptide blastP cutoff',
        #                 value=n_unique_above_peptide_cutoff)

    print(stats)
    stats.to_excel('statistics.xlsx')



if option_process_clustering:
    ''' get clustering data and create a new dataframe
    
    '''

    # check if peptide dataframe exists, otherwise read excel
    try:
        peptide_above_cutoff
        stats
    except NameError:
        try:
            print('Loading peptide dataframe...')
            peptide_above_cutoff = pd.read_excel('peptide_above_cutoff.xlsx', index_col=0)
            stats = pd.read_excel('statistics.xlsx')
            print('Done!')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')
            exit()

    #print(peptide_above_cutoff.head())
    samples = [sample for sample in dict_sample_fastq]

    # convert all counts to relative frequency
    print('Calculating relative frequency of peptides...')
    dna_df = pd.read_excel('dna_dataframe.xlsx')
    sum_counts = dna_df.loc[:, samples].sum(axis=0)

    for sample in samples:
        peptide_above_cutoff.loc[:, sample] = \
        peptide_above_cutoff.loc[:, sample] / np.float64(sum_counts.loc[sample])
    
    # load peptide-cluster data
    print('Loading peptide-cluster matrix...')
    # df_pep_cluster = pd.read_csv('clustering/peptide_cluster.txt', sep='\t',
    #                              index_col=0)
    df_pep_cluster = pd.read_excel('peptide_cluster.txt',
                                   index_col=0)

    print(df_pep_cluster.head())

    # load cluster-consensus data
    print('Loading peptide-cluster matrix...')
    # df_cluster_consensus = pd.read_csv('clustering/cluster_consensus.txt', 
    #     sep='\t', index_col=0)
    df_cluster_consensus = pd.read_excel('cluster_consensus.txt', 
                                         index_col=0)
    print(df_cluster_consensus.head())

    # align peptide-cluster to peptides
    peptide_above_cutoff = pd.merge(peptide_above_cutoff, 
        df_pep_cluster, left_index=True, right_index=True)

    # create epitope df
    epitope_df = peptide_above_cutoff[samples].groupby(
        peptide_above_cutoff['cluster']).sum()

    # Calculate the  sum of frequency of a peptide in all samples
    freq_sum = peptide_above_cutoff.loc[:, samples].sum(axis=1)
    freq_sum.rename('sum', inplace=True)
    peptide_above_cutoff = peptide_above_cutoff.join(freq_sum, how='outer')

    # find the most frequent peptide in each cluster
    for cluster_n in epitope_df.index:
        peps_in_cluster = peptide_above_cutoff.loc[:, 'cluster'] == cluster_n
        idxmax = peptide_above_cutoff.loc[peps_in_cluster, 'sum'].idxmax(axis=0)
        epitope_df.at[cluster_n, 'Most frequent peptide'] = idxmax
        #epitope_df.set_value(cluster_n, 'Most frequent peptide', value=idxmax)
        #print(cluster_n, idxmax)

    # add consensus to epitope df
    epitope_df = pd.merge(epitope_df, df_cluster_consensus, left_index=True,
        right_index=True)

    print(epitope_df)
    # merge with peptide_df using the most frequent peptide of each cluster
    a = peptide_above_cutoff.loc[:, ['blastP hit',
        'blastP identity ratio', 'blastP positive ratio', 'blastP query length',
        'blastP expect', 'blastP match name', 'blastP gene id', 
        'blastP query']]
    epitope_df = pd.merge(epitope_df, a, left_on='Most frequent peptide', right_index=True)
    print(epitope_df)

    # write to EXCEL
    epitope_df.to_excel('epitope_df.xlsx')


if option_blastp_stats:
    print('blastP stats')

    # check if DNA dataframe exists, otherwise read excel
    try:
        peptide_df
        stats
    except NameError:
        try:
            print('Loading peptide dataframe...')
            peptide_df = pd.read_excel('peptide_dataframe.xlsx')
            stats = pd.read_excel('statistics.xlsx')
            print('Done!')
        except FileNotFoundError:
            print('ERROR!! Dataframe excel file not found!!!')   
            exit()


    # normalize peptide reads
    print('Normalizing peptide reads...', time.asctime())
    norm_peptide_df = pd.DataFrame()
    norm_peptide_df = peptide_df.copy()
    for sample in dict_sample_fastq.keys():
        norm_peptide_df.loc[:, sample] = peptide_df.loc[:, sample] / stats.loc[sample, "# reads above DNA cutoff"]
        # norm_peptide_df.loc[:, sample] = norm_peptide_df.loc[:, sample] / norm_peptide_df.loc[:, sample].sum()
    print(norm_peptide_df.head())
    
    norm_peptide_df.to_excel('norm_peptide_df.xlsx')