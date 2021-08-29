"""
    This script is part of Teixeira et al 2021
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
Clustering Script
Takes a file named peptides.txt as input
Input file has 1 peptide sequence per line
Cluster peptides according to similarity
Outputs: aligment files, consensus sequences

3rd party software needed:
- mafft
- hmmbuild
- hmmemmit
"""
import os
import subprocess
import shutil
import pandas as pd
from collections import OrderedDict
from Bio import SeqIO
from fuzzywuzzy import fuzz

PEPTIDES_TXT = "example/peptides.txt"

# path to 3rd party binaries
MAFFT = 'mafft'
HMMBUILD = 'hmmbuild'
HMMEMIT = 'hmmemit'

K = 4  # K-mer length
min_len = 6
score_cutoff = 80

class Cluster(object):

    def __init__(self, seq, name, K=4):
        self.name = name
        self.consensus = seq
        self.seqs = [seq]
        self.kmers = self.get_kmers(K=K)
        self.K = K

    def __len__(self):
        return len(self.seqs)

    def merge(self, cluster, update=False):
        ''' merge a new cluster into this one '''
        self.seqs += cluster.seqs

        if update:
            self.update_consensus()

        return

    def update_consensus(self):

        if len(self.seqs) == 1:
            self.consensus = self.seqs[0]
            return

        #
        # Update alignment if there is already one
        #
        try:
            # get seqs in the old alignment
            old_seqs = list(SeqIO.parse('align/{}.mafft.fasta'.format(self.name), 'fasta'))
            old_seqs = [str(a.seq).replace('-', '') for a in old_seqs]

            # get new_seqs that have to be added and create fasta
            new_seqs = set(self.seqs) - set(old_seqs)

            # check if new_seqs is empty and raise error
            if not new_seqs: raise(IOError)

            self.fasta_from_list(new_seqs, 'align/temp.fasta')

            # add new sequences to alignment
            os.system('{mafft} --quiet --reorder --thread 8 --add align/temp.fasta align/{0}.mafft.fasta > align/{0}.mafft.fasta'.format(self.name, mafft=MAFFT))
            os.remove('align/temp.fasta')
        #
        # Create a new alignment if it has never been done
        #
        except IOError:
            # create a fasta
            self.fasta_from_list(self.seqs, 'align/temp.fasta')

            # align
            os.system('{mafft} --quiet --reorder --thread -1 align/temp.fasta > align/{0}.mafft.fasta'.format(self.name, mafft=MAFFT))
            os.remove('align/temp.fasta')

        #
        # Create consensus using HMMER
        #

        # create HMM
        with open(os.devnull, 'wb') as devnull:
            subprocess.call(
                '{hmmbuild} --amino --cpu 8 align/{0}.hmm align/{0}.mafft.fasta'.format(self.name, hmmbuild=HMMBUILD), shell=True, stdout=devnull)
            subprocess.call(
                '{hmmemit} -o align/{0}.consensus -c align/{0}.hmm'.format(self.name, hmmemit=HMMEMIT), shell=True, stdout=devnull)

        # read consensus output
        consensus = ''
        with open('align/{0}.consensus'.format(self.name), 'r') as f:
            for line in f:
                if line[0] != '>':
                    line = line.replace('\n', '')
                    consensus = '{}{}'.format(consensus, line)

        consensus = consensus.replace('\n', '')
        # update self.consensus
        self.consensus = consensus
        self.kmers = self.get_kmers(K=self.K)

        return

    def fasta_from_list(self, seq_list, fasta_name):
        ''' create new fasta '''
        with open(fasta_name, 'w') as f:
            for seq in seq_list:
                f.write('>{0}\n{0}\n'.format(seq))

    def get_kmers(self, K=4):
        ''' get kmers in consensus '''
        n_kmers = len(self.consensus) - K + 1
        kmers = set(self.consensus[i:i + K] for i in range(n_kmers))
        return kmers


def get_kmer_dict(cluster_dict, K):
    ''' Get a cluster list and return a dict containing all clusters that have 
    a given kmer - kmer:[cluster_name]'''

    kmer_seq_dict = {}
    for name in cluster_dict.keys():
        for kmer in cluster_dict[name].kmers:

            try:
                kmer_seq_dict[kmer].append(name)
            except KeyError:
                kmer_seq_dict[kmer] = [name]

    return kmer_seq_dict



PEPTIDES_TXT = open(PEPTIDES_TXT)
WORKING_DIR = os.path.dirname(os.path.realpath(PEPTIDES_TXT.name))
os.chdir(WORKING_DIR)
try:
    os.mkdir('align')
except OSError:
    # os.system('rm -rf align')
    shutil.rmtree('align')
    os.mkdir('align')

#
# Load sequences
#
print('Reading dataframe')


seq_list = [seq.strip() for seq in PEPTIDES_TXT]
seq_list = [seq for seq in seq_list if len(seq) >= min_len]  # exclude small seqs
seq_list = sorted(seq_list, key=len, reverse=False)  # sort by length
seq_list = [Cluster(str(seq), i) for i, seq in enumerate(seq_list)]
print('# of peptides {}'.format(len(seq_list)))

# convert to dict so I can access cluster by name, used OrderedDict to keep sorting
cluster_dict = OrderedDict()
for a in seq_list:
    cluster_dict[a.name] = a

fuzz_partial_ratio = fuzz.partial_ratio

# Keep clustering until is not possible anymore
cycle = 1
while True:

    print('Cycle {} -  Clusters {}'.format(cycle, len(cluster_dict)))
    cycle += 1
    merged = set()
    n_merges = 0
    for iterator_i in list(cluster_dict.keys()):
        if iterator_i in merged:
            continue

        # get the Kmer dict
        kmer_dict = get_kmer_dict(cluster_dict, K=K)

        # cluster that I should compare with I
        compare_with_I = set(name for kmer in cluster_dict[iterator_i].kmers for name in kmer_dict[kmer])
        compare_with_I.remove(iterator_i)

        for iterator_j in compare_with_I:

            score = fuzz_partial_ratio(cluster_dict[iterator_i].consensus, 
                                       cluster_dict[iterator_j].consensus)

            # merge cluster if above cutoff
            if score >= score_cutoff:

                # merge
                cluster_dict[iterator_i].merge(cluster_dict[iterator_j])

                # delete old
                del cluster_dict[iterator_j]
                
                # add to merged list
                merged.add(iterator_j)

                # merge count
                n_merges += 1


        # update cluster I after merges
        cluster_dict[iterator_i].update_consensus()
        print('Cluster {} length {} consensus {}'.format(iterator_i, len(cluster_dict[iterator_i].seqs), cluster_dict[iterator_i].consensus))
        print('# clusters {}'.format(len(cluster_dict)))
        print('consensus: {}'.format(cluster_dict[iterator_i].consensus))


    # if no merges happened, break
    #print(n_merges)
    if n_merges == 0:
        break

#
# sort by cluster number of peptides
#
cluster_dict = OrderedDict(sorted(cluster_dict.items(), key=lambda x: len(x[1].seqs), reverse=True))

shutil.rmtree('align')
os.mkdir('align')
new_ordered_clusters = OrderedDict()
for i, cluster in enumerate(cluster_dict):
    new_ordered_clusters[i] = cluster_dict[cluster]
    new_ordered_clusters[i].name = i
    new_ordered_clusters[i].update_consensus()

cluster_dict = new_ordered_clusters


#
# consensus
#
with open('cluster_consensus.txt', 'w') as consensus_file:
    print('cluster\tn_seqs\tconsensus')
    consensus_file.write('cluster\tn_seqs\tconsensus\n')
    for cluster in cluster_dict.keys():
        print('{}\t{}\t{}'.format(cluster, len(cluster_dict[cluster].seqs), cluster_dict[cluster].consensus))
        consensus_file.write('{}\t{}\t{}\n'.format(cluster, len(cluster_dict[cluster].seqs), cluster_dict[cluster].consensus))

df = pd.read_csv('cluster_consensus.txt', sep="\t")
df.to_excel('cluster_consensus.xlsx')
os.remove('cluster_consensus.txt')

#
# output peptide - cluster
#
with open('peptide_cluster.txt', 'w') as peptide_cluster_file:
    print('peptide\tcluster')
    peptide_cluster_file.write('peptide\tcluster\n')
    for cluster in cluster_dict.keys():
        for peptide in cluster_dict[cluster].seqs:
            print('{}\t{}'.format(peptide.replace('\n', ''), cluster))
            peptide_cluster_file.write('{}\t{}\n'.format(peptide.replace('\n', ''), cluster))
df = pd.read_csv('peptide_cluster.txt', sep="\t")
df.to_excel('peptide_cluster.xlsx')
os.remove('peptide_cluster.txt')