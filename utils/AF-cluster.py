import numpy as np
import argparse
import pandas as pd
import os
from polyleven import levenshtein
from sklearn.cluster import DBSCAN
import warnings
from Bio import BiopythonWarning
with warnings.catch_warnings():
    warnings.simplefilter('ignore', BiopythonWarning)
    from Bio import SeqIO

def load_fasta(fil):
    seqs, IDs =[], []
    with open(fil) as handle:
        with warnings.catch_warnings():
            warnings.simplefilter("ignore", BiopythonWarning)
            for record in SeqIO.parse(handle, "fasta"):
                seq = ''.join([x for x in record.seq])
                IDs.append(record.id)
                seqs.append(seq)
    return IDs, seqs

def write_fasta(names, seqs, outfile='tmp.fasta'):
        with open(outfile,'w') as f:
                for nm, seq in list(zip(names, seqs)):
                        f.write(">%s\n%s\n" % (nm, seq))

def consensusVoting(seqs):
    ## Find the consensus sequence
    consensus = ""
    residues = "ACDEFGHIKLMNPQRSTVWY-"
    n_chars = len(seqs[0])
    for i in range(n_chars):
        baseArray = [x[i] for x in seqs]
        baseCount = np.array([baseArray.count(a) for a in list(residues)])
        vote = np.argmax(baseCount)
        consensus += residues[vote]
    return consensus

def encode_seqs(seqs, max_len=108, alphabet=None):
    if alphabet is None:
        alphabet = "ACDEFGHIKLMNPQRSTVWY-"
    arr = np.zeros([len(seqs), max_len, len(alphabet)])
    for j, seq in enumerate(seqs):
        for i,char in enumerate(seq):
            for k, res in enumerate(alphabet):
                if char==res:
                    arr[j,i,k]+=1
    return arr.reshape([len(seqs), max_len*len(alphabet)])

if __name__=='__main__':

    p = argparse.ArgumentParser(description=
    """
    Cluster sequences in a MSA using DBSCAN algorithm and write .a3m file for each cluster.
    Assumes first sequence in fasta is the query sequence.

    H Wayment-Steele, 2022
    """)
    p.add_argument("keyword", action="store", help="Keyword to call all generated MSAs.")
    p.add_argument("-i", action='store', help='fasta/a3m file of original alignment.')
    p.add_argument("-o", action="store", help='name of output directory to write MSAs to.')
    p.add_argument('--eps_val', action='store', type=float, help="Use single value for eps instead of scanning.")
    p.add_argument('--resample', action='store_true', help='If included, will resample the original MSA with replacement before writing.')
    p.add_argument("--gap_cutoff", action='store', type=float, default=0.25, help='Remove sequences with gaps representing more than this frac of seq.')
    p.add_argument('--min_eps', action='store',default=3, help='Min epsilon value to scan for DBSCAN (Default 3).')
    p.add_argument('--max_eps', action='store',default=20, help='Max epsilon value to scan for DBSCAN (Default 20).')
    p.add_argument('--eps_step', action='store',default=.5, help='step for epsilon scan for DBSCAN (Default 0.5).')
    p.add_argument('--min_samples', action='store',default=3, help='Default min_samples for DBSCAN (Default 3, recommended no lower than that).')
    p.add_argument('--max_clust', action='store',default=50, help='Janniks way of finding hyperparamters with a maximum of this many clusters (Default 50).')
    args = p.parse_args()
    os.makedirs(args.o, exist_ok=True)
    IDs, seqs = load_fasta(args.i)
    seqs = [''.join([x for x in s if x.isupper() or x=='-']) for s in seqs] # remove lowercase letters in alignment
    df = pd.DataFrame({'SequenceName': IDs, 'sequence': seqs})
    query_ = df.iloc[:1]
    df = df.iloc[1:]
    if args.resample:
        df = df.sample(frac=1)
    L = len(df.sequence.iloc[0])
    N = len(df)
    df['frac_gaps'] = [x.count('-')/L for x in df['sequence']]
    former_len=len(df)
    df = df.loc[df.frac_gaps<args.gap_cutoff]
    new_len=len(df)
    ohe_seqs = encode_seqs(df.sequence.tolist(), max_len=L)
    n_clusters=[]
    eps_test_vals=np.arange(args.min_eps, args.max_eps+args.eps_step, args.eps_step)
    if args.eps_val is None: # performing scan
        selectable_eps = []
        for eps in eps_test_vals:
            #testset = encode_seqs(df.sample(frac=0.25).sequence.tolist(), max_len=L)
            clustering = DBSCAN(eps=eps, min_samples=args.min_samples).fit(ohe_seqs)
            n_clust = len(set(clustering.labels_))
            n_not_clustered = len(clustering.labels_[np.where(clustering.labels_==-1)])
            if n_clust<=args.max_clust:
                selectable_eps.append(eps)
                n_clusters.append(n_clust)
        eps_to_select = selectable_eps[np.argmax(n_clusters)]
        print(f"Selected eps: {eps_to_select:.1f} resulting in {max(n_clusters)} clusters.")

    else:
        eps_to_select = args.eps_val
    # perform actual clustering
    clustering = DBSCAN(eps=eps_to_select, min_samples=args.min_samples).fit(ohe_seqs)
    df['dbscan_label'] = clustering.labels_
    clusters = [x for x in df.dbscan_label.unique() if x>=0]
    unclustered = len(df.loc[df.dbscan_label==-1])
    avg_dist_to_query = np.mean([1-levenshtein(x, query_['sequence'].iloc[0])/L for x in df.loc[df.dbscan_label==-1]['sequence'].tolist()])
    avg_dist_to_query = np.mean([1-levenshtein(x, query_['sequence'].iloc[0])/L for x in df.loc[df.dbscan_label!=-1]['sequence'].tolist()])
    cluster_metadata=[]
    for clust in clusters:
        tmp = df.loc[df.dbscan_label==clust]
        cs = consensusVoting(tmp.sequence.tolist())
        avg_dist_to_cs = np.mean([1-levenshtein(x,cs)/L for x in tmp.sequence.tolist()])
        avg_dist_to_query = np.mean([1-levenshtein(x,query_['sequence'].iloc[0])/L for x in tmp.sequence.tolist()])
        tmp = pd.concat([query_, tmp], axis=0)
        cluster_metadata.append({'cluster_ind': clust, 'consensusSeq': cs, 'avg_lev_dist': '%.3f' % avg_dist_to_cs, 
            'avg_dist_to_query': '%.3f' % avg_dist_to_query, 'size': len(tmp)})
        write_fasta(tmp.SequenceName.tolist(), tmp.sequence.tolist(), outfile=args.o+'/'+args.keyword+'_'+"%03d" % clust+'.a3m')
