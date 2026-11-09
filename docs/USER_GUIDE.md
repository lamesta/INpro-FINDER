# INpro FINDER 1.1.0: user guide

Contents:

1. What the tool does
2. Getting started
3. Input files
4. Models
5. How a gene is found
6. Result classes
7. Output files
8. Parameters
9. Exploratory mode
10. Validation
11. Limits and good practice
12. Command-line reference
13. Sharing the tool
14. Provenance and references

## 1. What the tool does

INpro FINDER searches genome sequences for genes of ice-nucleation proteins (INA, InaZ/InaK/InaV/InaX family). It does not rely on the genome annotation: every genome is translated in all six reading frames, and every open reading frame (ORF) is searched with profile hidden Markov models (profile HMMs) of the conserved INA terminal domains. Many INA genes are annotated as "hypothetical protein", split into pieces or not annotated at all; the tool finds them in the same way as annotated ones.

For every INA locus it reports the structural state (intact, interrupted, fragment), the predicted protein, the gene and its ±20 kb neighbourhood, and, if an annotation is supplied, how the locus is annotated.

The tool comes as two files that run the same code:

- `inafinder.html`: the browser tool. Open it, drop in genome files, press *Scan*. Nothing is uploaded; all computing happens in the browser tab.
- `inafinder.py`: the command-line tool, for many genomes or for scripts. It can also write the browser tool (`python3 inafinder.py viewer`).

## 2. Getting started

### Browser

1. Open `inafinder.html` in a current web browser (Chrome, Firefox, Edge or Safari). No installation and no internet connection are needed.
2. Drag the genome FASTA files (and, optionally, their annotation files) onto the file area, or click *choose files*.
3. Choose the model: *Bacterial INA*, *Fungal INA (Mortierellaceae)* or *Custom profiles*.
4. Press *Scan genomes*. Genomes are scanned in parallel in background threads; a progress bar is shown per genome.
5. Click a locus in the results table to see its architecture, profile matches and protein.
6. Download the tables and sequences, or *save report*: an HTML file that opens directly on the results and can be sent to colleagues.

A bacterial genome (about 5 Mb) takes about 10 s per processor core; a fungal genome of 40 Mb takes about 2 min.

### Command line

Requirements: Python 3.8 or later (standard library only) and Node.js 18 or later.

```
python3 inafinder.py scan genome1.fna genome2.fna.gz -o results/
python3 inafinder.py scan *.fna --model fungal -o results_fungi/
python3 inafinder.py scan R10.79.fasta --annotation R10.79.fasta=R10.79.gff3 -o results_r1079/
```

Without Node.js, use the browser tool: it gives the same results.

## 3. Input files

| file | formats | notes |
|---|---|---|
| genome | FASTA (`.fa .fna .fasta .fsa`, also `.gz`) | one file per genome; any number of contigs |
| annotation (optional) | GFF3 or GenBank (`.gff .gff3 .gbff .gbk .gb`, also `.gz`) | used only to compare the loci with the annotated genes and to draw neighbouring genes; sequence names must match the FASTA contig names |
| profiles (optional, custom model) | HMMER3 ASCII (`.hmm`, made with `hmmbuild`) | one or more profiles per file |

In the browser, file types are recognised from their content, and annotation files are matched to genomes by file name (for example `R10.79.fasta` and `R10.79.gff3`). The table of loaded files shows which annotation is used for which genome; check it before scanning.

If none of the sequence names in an annotation match the genome, the annotation is ignored and a warning is shown. If an annotation contains no CDS features, a warning is shown and every locus is reported as "not annotated".

## 4. Models

### Bacterial INA (default)

| profile | role | positions | built from |
|---|---|---|---|
| INA_N150 | N-terminal domain | 151 | first 150 aa of 42 representative INA |
| INA_C60 | C-terminal cap | 60 | last 60 aa of the same 42 representatives |

The 42 representatives were selected from 491 validated intact INA proteins of *Pseudomonas*, *Xanthomonas*, *Pantoea* and *Erwinia* (four-genera catalogue, September 2026) by clustering at 90 % identity and 80 % coverage (MMseqs2), aligned with MAFFT L-INS-i and built with `hmmbuild` (HMMER 3.4). The defaults reproduce the published four-genera analysis.

### Fungal INA (Mortierellaceae lineage I)

| profile | role | positions | built from |
|---|---|---|---|
| INAf_N | N-terminal region | 277 | 10 lineage I proteins of 4 species, from the start Met to the first GYGST repeat |
| INAf_Ccap | Cys-rich C-terminal cap | 45 | last 45 aa of 14 lineage I proteins |
| INA_C60 | bacterial-type cap | 60 | as above; for the one fungal INA with a bacterial-type cap |

The rules differ from the bacterial model (`models/fungal_model_rules.md`). The most important one: the INAf_N match must reach profile position 250, because the first half of the region (signal peptide and the part after it) also matches unrelated fungal proteins. The model was built and tested on Mortierellaceae only. Fungal "lineage II" proteins are not modelled.

### Custom profiles

Drop one or more `.hmm` files into the file area and choose *Custom profiles*. For each profile choose its role: N-terminal domain, C-terminal cap, or ignore. Default rules are set from the profile length (N: alignment from position ≤ 10 to ≥ 66 % of the profile; C: from position ≤ 15 to ≥ 80 % of the profile; i-Evalue ≤ 1e-10) and can be edited. On the command line, a custom model is a JSON file with the same fields as `models/fungal_model.json` (`python3 inafinder.py extract DIR` writes the model files as examples).

## 5. How a gene is found

1. **ORFs.** Each contig is translated in all six frames with the standard genetic code. Every stretch of at least 50 codons between two stop codons (or between a contig end and a stop codon) is an ORF. The ORFs are identical to those of `esl-translate -l 50` (HMMER/Easel).
2. **Profile search.** Every profile is searched against all ORFs in local alignment mode, following HMMER 3 (`hmmsearch`): MSV, Viterbi and Forward filters (P 0.02, 1e-3, 1e-5), sequences reported at E ≤ 1e-5, domains at conditional E ≤ 1e-3, E-values computed over the number of ORFs searched. The search is a re-implementation in JavaScript; scores agree closely with `hmmsearch` but are not identical (section 10).
3. **Canonical domains.** A domain match is canonical when it starts near the beginning of the profile, covers most of it, is significant (i-Evalue ≤ 1e-10) and lies at the correct end of the protein. Bacterial defaults: N domain from profile position ≤ 10 to ≥ 100, starting within the first 300 aa; C cap from position ≤ 15 to ≥ 50, ending within the last 200 aa.
4. **Candidates and loci.** An ORF is a candidate if it has a canonical N domain or C cap, or at least 5 weaker matches and a length of at least 300 aa. Candidate ORFs on the same strand within 8 kb form one locus. The representative ORF of a locus is the one with most canonical domains, then the longest protein.
5. **Start codon.** If an annotated, complete CDS ends at the same stop codon, its start is used. Otherwise the start is the first Met between 60 aa before and 30 aa after the first profile match. If there is none, the bacterial model takes the first Met of the ORF (as in the published analysis), and the fungal model keeps the whole ORF (`ORF_start_no_Met`), because fungal genes are often cut by contig ends or introns.
6. **Repeat test.** INA proteins have a central domain of 16-aa repeats. The central part of the protein (from residue 151 to 50 residues before the end) is cut into 16-aa blocks at each of the 16 possible offsets; a consensus is made per position, and the mean identity of the blocks to their consensus is the *repeat identity* (best offset; at least 5 blocks and 240 aa are needed). Bacterial INA score 0.49 to 0.79 (median 0.73 in 491 validated proteins); the threshold is 0.45.
7. **Annotation comparison** (if an annotation is loaded). Overlapping annotated CDS on the same strand are compared with the predicted gene (section 6).

## 6. Result classes

### Structural class

| class | meaning |
|---|---|
| intact_full | canonical N domain and C cap in one ORF, repeat identity ≥ 0.45, protein ≥ 800 aa |
| intact_short | as intact_full, 500 to 799 aa |
| N_and_C_low_repeat_or_short | N domain and C cap in one ORF, but weak repeats or shorter than 500 aa |
| interrupted_N_and_C_in_separate_ORFs | N domain and C cap in different ORFs of one locus, in the right order: a stop codon or frameshift inside the gene (a mutation or a sequencing error) |
| intact_noN_canonical_domain | C cap and repeats, ≥ 500 aa, no canonical N domain |
| intact_noC_canonical_domain | N domain and repeats, ≥ 500 aa, no canonical C cap |
| fragment_N_only | N domain only |
| fragment_C_only | C cap only |
| repeat_only | profile matches and repeats (≥ 300 aa), no canonical terminal domain; low confidence |
| weak_repeat_hit | profile matches only; low confidence |
| repeat_only_no_profile_hit | exploratory mode only: 16-aa repeats, no profile match; a hypothesis (section 9) |

### Annotation state

| state | meaning |
|---|---|
| annotated_intact_CDS | an annotated CDS ends at the same stop codon; its start is used |
| annotated_CDS_different_start_or_stop | an annotated CDS overlaps the gene (more than 30 % of either), with another start or stop |
| annotated_pseudo_or_partial | the overlapping annotated feature is marked pseudo or partial |
| no_CDS_annotated | no annotated CDS overlaps the gene: an unannotated INA gene |
| not_checked | no annotation was loaded |

### Other columns

- **contig_edge**: `5prime_open`, `3prime_open` or `both_ends` when the ORF reaches the end of a contig without a stop codon. Such genes may continue beyond the assembled sequence; many fragments in draft assemblies are of this kind.
- **period4_ratio**: identity at a 4-aa period divided by identity at 16 aa. INA give about 0.5 (0.50 to 0.59 in 491 bacterial INA); low-complexity repeats with a period of 4 aa or less give about 1.
- **start_source**: `annotation_same_stop`, `Met_rule`, `ORF_start_no_Met` or, for exploratory calls, `first_Met_in_ORF`.

## 7. Output files

| file | content |
|---|---|
| `loci.tsv` | one row per locus: genome, locus id, contig, strand, CDS start and end, length, class, repeat identity, number of 16-aa blocks, repeat consensus, N and C profile and i-Evalue, number of ORFs, locus start and end, start source, contig edge, period-4 ratio, annotation state, annotated locus tag, protein id and product, model |
| `proteins.faa` | predicted proteins |
| `genes.fna` | predicted genes (CDS, including the stop codon) |
| `ina.gff3` | gene and CDS features, with class and E-values as attributes |
| `loci_20kb.fna` | the gene with 20 kb on each side |
| GenBank `.gbk` (zip in the browser, `genbank/` on the command line) | one file per locus: the ±20 kb region with the INA gene and the neighbouring genes, ready for clinker |
| `run.json` | the run record: software version, model, all parameters, checksums (SHA-256) of every input file and profile, per-genome numbers and warnings |
| `report.html` | the browser tool with the results built in; opens on the results |

Coordinates are 1-based and inclusive, on the forward strand of the contig; CDS coordinates in the tables, GFF3 and GenBank files include the stop codon.

Neighbouring genes in the GenBank files come from the annotation. Without an annotation they are simple ORF calls (ATG, GTG or TTG start, at least 100 codons, non-overlapping), labelled as such in a `/note`. In the GenBank files, a CDS that starts without a known start codon or runs into a contig end is marked partial (`<` or `>`). Low-confidence classes get a different product name ("weak INA profile and repeat match", "16-aa repeat protein without INA profile match").

## 8. Parameters

All parameters can be changed in the browser (*Parameters*) and are saved in `run.json`. *reset to model defaults* restores them.

| parameter | default | meaning |
|---|---|---|
| per profile: alignment starts at profile position ≤ | N 10, C 15 | canonical rule |
| per profile: reaches position ≥ | N 100, C 50 (bacterial) | canonical rule |
| per profile: i-Evalue ≤ | 1e-10 | canonical rule |
| per profile: position in protein | N within first 300 aa; C within last 200 aa | canonical rule |
| minimum ORF length | 50 aa | ORFs searched |
| weak matches for a candidate | 5, with ORF ≥ 300 aa | candidate rule |
| 16-aa repeat identity ≥ | 0.45 | repeat rule |
| intact, full length ≥ | 800 aa | class rule |
| intact, short ≥ | 500 aa | class rule |
| repeats only ≥ | 300 aa | class rule |
| ORFs in one locus within | 8,000 bp | locus rule |
| start Met window | 60 aa before to 30 aa after the first match | start rule |
| no Met near first match | bacterial: first Met in the ORF; fungal: ORF start | start rule |
| exploratory: length ≥, repeat identity ≥, 4/16 ratio ≤ | 500 aa, 0.45, 0.75 | section 9 |
| search filters | MSV 0.02, Viterbi 1e-3, Forward 1e-5; E ≤ 1e-5; domain c-E ≤ 1e-3 | as `hmmsearch` |

Results obtained with changed parameters are not comparable with the published analysis; `run.json` records what was used.

## 9. Exploratory mode

The profiles describe the INA that are already known. To look for proteins with INA-like repeats but without known terminal domains, tick *exploratory* (browser) or use `--exploratory` (command line). The tool then also reports ORFs without any profile match when:

- the protein is at least 500 aa long;
- the 16-aa repeat identity is at least 0.45;
- the repeats are not explained by a period of 4 aa or less (period-4 ratio ≤ 0.75);
- the ORF does not overlap a profile-based locus; of overlapping ORFs in different frames only the longest is kept (repetitive DNA is periodic in every frame).

These calls are labelled `repeat_only_no_profile_hit` and are hypotheses, not INA genes. Measured on 1,062 bacterial genomes, the mode reports 2 proteins in total, both unrelated to INA; on 102 Mortierellaceae genomes, 11 proteins in 10 genomes, none of them a known INA. With a minimum of 300 aa the same genomes gave 3,752 and 559 calls, which is why the default is 500 aa. The mode does not find the fungal lineage II proteins: their repeats have a 4-aa period.

## 10. Validation

Full numbers: `TEST_REPORT.md`. In short:

- **Bacterial model, 1,062 genomes against the hmmsearch catalogue:** 774 of 775 loci found with the same class, 771 with identical coordinates, all annotation states identical, no locus in 300 genomes without INA. Domain scores: r = 0.998 with HMMER, median 2.0 bits lower. The one missed locus is a 73-aa fragment whose cap alignment starts one profile position beyond the rule.
- **Fungal model, 102 Mortierellaceae genomes:** all 21 confident lineage I loci found; 20 of 20 found again with profiles built without the species concerned; no canonical locus in 79 genomes without INA.
- **Code:** ORFs identical to `esl-translate`, repeat test identical to the Python pipeline, browser and command line give the same results.

## 11. Limits and good practice

- **Reference bias.** The profiles are built from known INA. Proteins with unrelated terminal domains are not found; the exploratory mode is a partial check only.
- **Not HMMER.** The search follows HMMER 3 closely but is not HMMER. For publication, confirm key loci with `hmmsearch` using the same profiles (the `.hmm` files are written by `inafinder.py extract`).
- **Introns.** Six-frame translation does not splice. A fungal INA gene with an intron appears as fragments or as an interrupted locus.
- **Draft assemblies.** Genes at contig ends are fragments of the assembly, not necessarily of the genome; check the *contig edge* column before interpreting a fragment biologically.
- **Interrupted genes.** A stop codon or frameshift inside a gene can be a real mutation or a sequencing error; confirm with reads or a second assembly.
- **Start codons** of unannotated genes are predictions.
- **Fungal model.** Built from four Mortierellaceae species; its behaviour on distant fungi is unknown.
- **Memory.** Large eukaryotic genomes produce several hundred thousand ORFs; scan fewer genomes at a time if the browser becomes slow.

## 12. Command-line reference

```
python3 inafinder.py viewer [-o inafinder.html]
python3 inafinder.py scan GENOME [GENOME ...] -o OUTDIR [options]
python3 inafinder.py extract DIR
```

| option of `scan` | meaning |
|---|---|
| `--model bacterial` / `fungal` / `my_model.json` | model (default bacterial) |
| `--annotation GENOME=ANNOTATION` | annotation for one genome (repeatable) |
| `--annotation-dir DIR` | annotation files matched to genomes by file-name stem |
| `--jobs N` | genomes scanned in parallel (default: number of CPUs minus 1) |
| `--exploratory` | add exploratory calls (section 9) |
| `--flank N` | flank in bp for the region files (default 20,000) |

`extract` writes the source code and the model files (profiles and JSON) carried inside `inafinder.py`.

## 13. Sharing the tool

Both files are self-contained; send `inafinder.html` to colleagues who will use the browser and `inafinder.py` to those who use the command line. Some e-mail systems block attachments ending in `.py`; rename the file to `inafinder.py.txt` (the receiver renames it back) or share it through a cloud folder. A saved report (`report.html`) can be sent on its own: it contains the results and the tool.

## 14. Provenance and references

- Bacterial profiles and rules: PSAIN four-genera INA catalogue (Aarhus University, September 2026): `ina_locus_pipeline/analyses/four_genera_complete_ncbi_20260930`.
- Fungal profiles and catalogue: `INA_EVOLUTION/01_fungal_ina` (October 2026).
- Eddy SR (2011) Accelerated profile HMM searches. *PLoS Comput Biol* 7: e1002195 (HMMER 3).
- Katoh K, Standley DM (2013) MAFFT multiple sequence alignment software version 7. *Mol Biol Evol* 30: 772–780.
- Steinegger M, Söding J (2017) MMseqs2 enables sensitive protein sequence searching for the analysis of massive data sets. *Nat Biotechnol* 35: 1026–1028.
- Gilchrist CLM, Chooi YH (2021) clinker & clustermap.js: automatic generation of gene cluster comparison figures. *Bioinformatics* 37: 2473–2475.
- Green RL, Warren GJ (1985) Physical and functional repetition in a bacterial ice nucleation gene. *Nature* 317: 645–648.
