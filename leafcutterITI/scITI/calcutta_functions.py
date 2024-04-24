import pandas as pd
import collections
import time
from pathlib import Path
import numpy as np
import scipy.sparse as sp

import gzip

from collections import OrderedDict, namedtuple

REVERSER=str.maketrans("AGCT","TCGA")

# this scipts contain function from calcultta that will be used in scITI
# functions in this file were written by Dr. David A. Knowles

def smart_open(filename, *argene_set, **kwargene_set):
    return gzip.open(filename, *argene_set, **kwargene_set) if filename.suffix==".gz" else open(filename, *argene_set, **kwargene_set)


def read_alevin_ec(fn):
    """Read gene_eqclass.txt.gz from `alevin quant --dump-eqclasses`"""

    ecs = collections.OrderedDict()

    with smart_open(fn) as f: 
        for i,l in enumerate(f):
            if i==0: # first line gives number of features (normally genes but can be transcripts)
                num_genes = int(l.decode().strip())
                continue
            if i==1: # second line gives number of ECs
                num_ec = int(l.decode().strip())
                continue
            l = l.decode().strip().split()
            l = [int(g) for g in l]
            ec_idx = l[-1] # last element of line is EC index
            gene_idx = l[:-1]
            ecs[ec_idx] = gene_idx
    return num_genes, num_ec, ecs




def to_coo(x, shape = None):
    """x is a list, where each item corresponds to a different row. The item for a row gives the column IDs for the non-zero entries. """
    nnz = sum([len(g) for g in x])

    indices = np.zeros((2,nnz), dtype=int) # cell then EC idx

    nz_idx = 0
    for row_idx,col_ids in enumerate(x): 
        nhere = len(col_ids)
        indices[0,nz_idx:nz_idx+nhere] = row_idx
        indices[1,nz_idx:nz_idx+nhere] = list(col_ids) # might be a set
        nz_idx += nhere
    
    return sp.coo_matrix((np.ones(nnz), indices), shape = shape)


def sparse_sum(x, dim):
    return np.squeeze(np.asarray(x.sum(dim)))


def get_fasta(fasta_file, first_field = False):
    with smart_open(fasta_file) as f:
        F = f.read().split(">")
    dic = OrderedDict()
    for x in F:
        x = x.split("\n")
        seq =  "".join("".join(x[1:]).split("\r"))
        seq=seq.upper()
        if len(x) <= 1: continue
        if first_field:
            dic[x[0].split()[0].strip()] = seq
        else:
            dic[x[0].strip()] = seq

    return dic

def get_transcript_lengths(fasta_file):
    cdna = get_fasta(fasta_file, first_field = True)
    return OrderedDict([(transcript,len(seq)) for transcript,seq in cdna.items()])



def get_feature_weights(features, transcript_lengths_dic,  fragment_size = 300):
    feature_lengths = np.array([transcript_lengths_dic[g] for g in features])
    eff_lens = np.array([ 
        ((g-fragment_size) if (g>fragment_size) else g) 
        for g in feature_lengths ]) # discontinous :(
    return feature_lengths, 1. / eff_lens


def EM(counts, ec_transcript_mat, w, iterations = 30):
    
    n_transcripts = len(w)
    alpha = np.full(n_transcripts,1.0/n_transcripts) # initialize to uniform
    
    alpha_w = ec_transcript_mat.copy()

    for i in range(iterations):
        alpha_w.data = (alpha * w)[ec_transcript_mat.col] # alpha_w[e,t] = ec_transcript_mat[e,t] alpha_t w_t
        ec_sums = sparse_sum(alpha_w,1) # ec_sums[e] = sum_{t \in e} alpha_t w_t
        z = sp.diags(counts / ec_sums) @ alpha_w 
        alpha_new = sparse_sum(z,0)
        alpha_new /= alpha_new.sum()
        #print(i,np.mean(np.abs(alpha - alpha_new)),end="\r")
        alpha = alpha_new
        
    return alpha









