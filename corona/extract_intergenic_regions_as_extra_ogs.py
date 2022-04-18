#!/usr/bin/env python

# script to extract non-covered parts of the sars-cov-2 wuhan sequence as extra groups.
# essentially take chunks of at least contigious 30bp in the assembly that are not covered 
# by any of the exported oma groups. (intergenic regions and ORF8 and ORF10 genes, which are
# not part of any exported oma group.

import Bio.SeqRecord
import Bio.SeqIO
import itertools
import gzip


with gzip.open('MN908947.embl.gz','rt') as fh:
    rec = next(Bio.SeqIO.parse(fh, 'embl'))
fff = [f for f in rec.features if f.type == "CDS" and f.qualifiers['gene'] not in ('ORF8', 'ORF10')]

cov = sum([f.location for f in fff])

inf, b, s = False, [], 0
for i in range(len(rec)):
    if inf and i in cov:
        b.append(slice(s, i))
        inf = False
    elif not inf and i not in cov:
        s = i
        inf = True
if inf: 
    b.append(slice(s, len(rec)))
bb = [x for x in b if x.stop - x.start > 30]
for nr, b in enumerate(bb):
    with open('extra_sars_ogs/OG{}.fa'.format(1000+nr), 'wt') as fout:
        trec = Bio.SeqRecord.SeqRecord(rec.seq[b], id="SARS20100{}".format(nr), description=str(b))
        Bio.SeqIO.write([trec], fout, 'fasta')
