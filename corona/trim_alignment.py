import Bio.AlignIO
import Bio.Align
from Bio.SeqRecord import SeqRecord
import collections
import math

def load_msa(fn):
    if fn.endswith('.phy'):
        format = 'phylip-relaxed'
    elif fn.endswith('.fa'):
        format = 'fasta'
    else:
        raise UnkownFormatError('unknown format for '+fn)
    with open(fn, 'rt') as fh:
        msa = next(Bio.AlignIO.parse(fn, format))
    return msa


def write_msa(fn, msa):
    with open(fn, 'wt') as fh:
        Bio.AlignIO.write(msa, fh, 'phylip-relaxed')


def count_nucs(data):
    c = collections.Counter(data)
    valid = sum(c[x] for x in ('ATCGN'))
    return valid

def trim(msa, min_residue):
    colcnt = [0] * msa.get_alignment_length()
    for rec in msa:
        for pos, nuc in enumerate(rec.seq):
            if nuc in 'ATCGN':
                colcnt[pos] += 1

    keep = [i for i in range(msa.get_alignment_length()) if colcnt[i] > min_residue]
    trimmed = []
    for rec in msa:
        seq = "".join([nuc for c, nuc in enumerate(rec.seq) if colcnt[c] >= min_residue])
        trimmed.append(SeqRecord(seq, id=rec.id))
    trimmed = Bio.Align.MultipleSeqAlignment(trimmed)
    return keep, trimmed

def filter_taxa(msa, min_residue):
    filtered = Bio.Align.MultipleSeqAlignment(filter(lambda taxon: count_nucs(taxon) > min_residue, msa))
    return filtered

class UnknownFormatError(Exception):
    pass


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description="sample part of the alignment that contains enough data, and through out species which have too little data")
    parser.add_argument('alignment', help="path to multiple sequence alignment")
    parser.add_argument('--min-per-col', type=int, help="Min nr of taxa that need to have a nuc at a column to be included. Defaults to ceil(nr_taxa*0.3)")
    parser.add_argument('--min-res-per-species', default=400, type=int, help="Minimum number of residues for a taxon in the trimmed alignment to not be kicked out. Defaults to 400")
    parser.add_argument('--out', help="Outfile of trimmed alignment")
    conf = parser.parse_args()

    msa = load_msa(conf.alignment)
    if conf.min_per_col is None:
        conf.min_per_col = math.ceil(0.3*len(msa))
    if conf.out is None:
        conf.out = conf.alignment+".trimmed"

    print("Loaded MSA ({}x{}). Filter cols with less than {} residue"
          .format(len(msa), msa.get_alignment_length(), conf.min_per_col))
    keep, trimmed_msa = trim(msa, conf.min_per_col)
    print("  after filtering columns: {}x{}".format(len(trimmed_msa), trimmed_msa.get_alignment_length()))
    filtered = filter_taxa(trimmed_msa, conf.min_res_per_species)
    print("  after filtering taxa: {}x{}".format(len(filtered), filtered.get_alignment_length()))
    write_msa(conf.out, filtered)
