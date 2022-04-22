# Corona Virus dataset analysis

The following steps were performed in order to obtain the corona virus tree with Read2tree.

The bash scripts in this document refer to some script in the read2tree_paper repository. The path to the repo should be set as an environment variable prior to run the code snippets in this document:
```shell
export paper_repo="/full/path/to/paper/repo"
```

## Data preparing step

### Fetching the data from corona-omabrowser
- Export marker genes from https://corona.omabrowser.org/oma/export_markers for all the _Orthocoronavirinae_.
- Download the cds sequences from the download section: https://corona.omabrowser.org/All/viruses.cdna.fa.gz
- Download the metadata for species from the download section: https://corona.omabrowser.org/All/oma-species.txt

Expand the tarball of the marker_genes you obtained.


### Extra markers for non-coding parts of SARS-COV-2 
To increase the ability to classify minute differences in the sars-cov-2 samples, we created extra oma groups 
with the non-coding parts of the sars-cov-2 reference genome (MN908947 - Wuhan-Hu-1). 

```shell
wget -O - "https://www.ebi.ac.uk/ena/browser/api/embl/MN908947.3?download=true" | gzip > MN908947.embl.gz 
python $paper_repo/corona/extract_intergenic_regions_as_extra_ogs.py
``` 
This will add 4 extra group with intergenic regions from SARS2 reference genome into `extra_sars_ogs/` subdirectory.
Link them to the marker genes and add the sequences to the viruses.cdna.fa.gz:

```shell
for f in extra_sars_ogs/*fa; do 
  ln $f marker_genes/
done
cat extra_sars_ogs/*fa | gzip >> viruses.cdna.fa.gz
```

### Obtain sars-cov-2 reads from sra
Next, we download as many sars-cov-2 samples as useful. We obtained sra accessions and metadata from nextstrain open:

- obtain Metadata (global): https://data.nextstrain.org/files/ncov/open/global/metadata.tsv.xz
- select different samples with sra accession and spanning all different Nextstrain_clades with the script:
  ```bash
  python $paper_repo/corona/subsample_nextstrain_covid_genomes_with_sra_accession.py --out subset.metadata.txt --nr-per-clade 200 metadata.tsv.xz 
  ```
  and extract the column with the SRA accession from that file. The file `used_accessions.txt` contains the file with the accessions we used in the analysis.
- download and convert sra-accession runs with
  ```bash
  python $paper_repo/corona/sradownloader --nogeo --noena --outdir reads used_accession.txt --threads 8
  ```
- clean the results: remove single paired-end read files (not ending in _1|2.fastq.gz), any broken file, etc from the `reads` subdirectory

## Run Read2tree

We run read2tree on a slurm cluster. The following steps explain how to map all reads onto the reference groups in parallel, build the concated alignment and infer the tree.

- Convert the reference species (exported groups) with read2tree
  ```shell
  mkdir corona_run
  cd corona_run
  Read2Tree --reference --standalone ../marker_genes/ --dna_reference ../viruses.cdna.fa.gz
  ```

- Run scripts/slurm_submit.py to generate jobfiles (and submit them to slurm). Adjust template of script according to needs.
  ```shell
  find ../reads -type f | sed -e 's/_.*//' | xargs python $paper_repo/submissions/slurm_submit.py ./ --only-jobfiles --species
  ```
  Note that the script assumes long-read read files to be in a subdirectory `ONT/`.

- Once all read sets are successfully mapped, the merge step will produce the concat_merge_dna.phy file. As a post-processing step, we recommend filtering uninformative columns and rows of this matrix:
  ```shell
  python $paper_repo/corona/trim_alignment.py --out trimmed_concat_merge_dna.phy concat_merge_dna.phy
  ```

- Relabel the sequences with more readable labels:
  ```shell
  python $paper_repo/corona/relabel_msa.py --out trimmed_concat_merge_dna_nice.phy --oma-map ../oma-species.txt --nextstrain ../subset.metadata.txt trimmed_concat_merge_dna.phy
  ```

- Reconstruct species phylogeny using FastTree:
  ```shell
  FastTree -nt trimmed_concat_merge_dna_nice.phy > corona.trimmed_dna.nwk
  ```
