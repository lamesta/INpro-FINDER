# INpro FINDER 1.1.0: test report

Date: 21 and 22 October 2026. Machine: macOS, 10 CPU cores, 16 GB RAM, Node.js 26, Python 3.12.
All numbers below come from the files in `validation/` and `tests/out/` of the package folder.

## 1. Automated test suite (`bash tests/run_all_tests.sh --full`)

| test | what it checks | result |
|---|---|---|
| orfs_vs_esl_translate | six-frame ORFs of R10.79 identical to `esl-translate -l 50` (Easel 0.49) | PASS |
| repeat_profile_parity | 16-aa repeat test identical to `repeat_profile` of the published Python pipeline (521 proteins: 491 bacterial INA, 30 fungal) | PASS |
| cli_r1079 | `inafinder.py scan` on R10.79 with its GFF3: one locus, intact_full, 1,342 aa, MGGFHP_01089, all output files written, GenBank parses with Biopython | PASS |
| single_file_viewer | `inafinder.py viewer` writes a file byte-identical to `inafinder.html` | PASS |
| browser_r1079 | browser tool in jsdom: load files, scan, all seven downloads, saved report reopens on the results | PASS |
| regression_pilot34 | 34 genomes rescanned; profile-based loci identical to the stored 1.1.0 results | PASS |

Result: 6 passed, 0 failed, 0 skipped.

## 2. Bacterial model against the hmmsearch catalogue

Reference: four-genera catalogue `structural_catalog/ina_structural_loci.tsv` (HMMER 3.4 `hmmsearch` on `esl-translate` ORFs, classified by `classify_hmm6f_v2.py`; 775 loci in 761 of 5,562 genomes).
Test set: 1,062 genomes (all 761 with a catalogue locus, 300 without, R10.79), with the PGAP GenBank annotation as in the catalogue.

| check | 1.0.0 | 1.1.0 |
|---|---|---|
| catalogue loci found with the same structural class | 774 / 775 | 774 / 775 |
| same CDS coordinates | 771 / 774 | 771 / 774 |
| same annotation state | 774 / 774 | 774 / 774 |
| loci in the 300 genomes without INA | 0 | 0 |
| R10.79 (not in the catalogue) | intact_full, 1,342 aa | intact_full, 1,342 aa |
| median run time per genome (one core, 8 in parallel) | 8.9 s | 10.0 s |

Differences between 1.0.0 and 1.1.0 over all 775 profile-based loci: 2, both a label only. In GCF_059716085.1 and GCF_059716145.1 a 63-aa fragment_C_only ORF contains no Met; 1.0.0 labelled its start `Met_rule`, 1.1.0 labels it `ORF_start_no_Met`. Coordinates and classes are unchanged.

Domain scores against HMMER (1.0.0; the search code is unchanged in 1.1.0): 460 matched domains in 27 genomes, Pearson r = 0.998, median difference (JavaScript minus HMMER) −2.0 bits; terminal (canonical) domains differ by a median of 0.7 bits. The one missed locus (GCF_030687105.1, a 73-aa fragment_C_only) and the other three differences are explained in `validation/validation_bacterial_summary.md`: domain boundaries come from the Viterbi alignment instead of HMMER's posterior decoding.

## 3. Fungal model against the Mortierellaceae catalogue

Reference: `INA_EVOLUTION/01_fungal_ina/fungal_ina_catalog.tsv` (46 entries; 30 confident INA in 23 assemblies; lineage I 21 loci in 16 assemblies, lineage II 9 loci).
Test set: all 102 Mortierellaceae assemblies (29 to 58 Mb), no annotation.

| check | result |
|---|---|
| lineage I confident loci found | 21 / 21 (intact_full 5, intact_short 1, intact_noN 2, fragment_N_only 4, fragment_C_only 7, repeat_only 2) |
| lineage I loci that run into a contig end | 11 / 21 |
| leave-one-species-out (4 species, 20 loci) | 20 / 20 with the same class |
| lineage II loci | 0 / 9 (not modelled) |
| non-confident catalogue entries found | 2 / 16 (locus0006 fragment_C_only; locus0010 fragment_N_only with a complete N-region, i-Evalue 2.5e-154) |
| canonical loci in the 79 assemblies without INA | 0 |
| low-confidence loci (repeat_only, weak_repeat_hit) without a catalogue entry | 11 in 11 assemblies (9 without INA, 2 with INA); all His/Lys-rich repeat proteins, not INA |
| median run time per genome | 117 s |

Details: `validation/validation_fungal_summary.md`, `models/fungal_model_rules.md`.

## 4. Exploratory mode (off by default)

Exploratory calls are ORFs without any profile match whose 16-aa repeat identity is ≥ 0.45 and whose 4-aa / 16-aa identity ratio is ≤ 0.75. Results of the same 1,062 bacterial and 102 fungal genomes:

| minimum length | bacterial calls (genomes) | fungal calls (genomes) |
|---|---|---|
| 300 aa | 3,752 (714 of 1,062) | 559 (102 of 102) |
| 500 aa (default) | 2 (2 of 1,062) | 11 (10 of 102) |

At 500 aa the 2 bacterial calls are annotated "hypothetical protein" (5,320 and 2,417 aa; repeat consensus TTAQVAALTTAQIAAL and TAQVAALTTAQVAALT) in genomes without INA. None of the 11 fungal calls overlaps a confident catalogue INA; 2 overlap non-confident catalogue entries (locus0005, locus0008).

Sensitivity of the same filters on known proteins (protein level, `validation/period_check.tsv`): 491 / 491 bacterial INA, 9 / 9 lineage I proteins of at least 500 aa, 0 / 4 lineage II proteins of at least 500 aa (their repeats have a 4-aa period).

The default of 500 aa was chosen from these numbers. Exploratory calls remain hypotheses.

## 5. Browser tool (jsdom, Node.js)

| test | result |
|---|---|
| R10.79 + GFF3, bacterial model, exploratory on | 1 locus: intact_full, 1,342 aa, annotated_intact_CDS MGGFHP_01089; no script errors; 7 downloads; saved report reopens with 1 locus; 13.5 s |
| GCA_001021685.1, fungal model, exploratory on (minimum then 300 aa) | fragment_N_only 418 aa (3′ open) and fragment_C_only 225 aa (5′ open), as in the catalogue (locus0003, locus0004); 7 exploratory calls; the 9 GenBank files parse with Biopython, partial ends marked with < and >; 187 s |
| custom profiles (INA_N150.hmm and INA_C60.hmm dropped as files), R10.79 and *X. vasicola* GCF_000772715.2 | R10.79 intact_full 1,342 aa; *X. vasicola* 0 loci; 19.7 s |
| architecture figure | renders with cairosvg at 700 and 1,180 px; legend wraps to the figure width |

The browser tool was tested in jsdom, not by hand in each browser. It uses standard web features only (Web Workers, Blob, FileReader).

## 6. Command line against the batch engine (5 genomes)

`inafinder.py scan` on GCA_030582555.1, GCA_037234985.1, GCF_030687105.1, GCF_000772715.2 and R10.79, each with its annotation: the 3 loci found (intact_full 1,525 aa; interrupted 195 aa, annotated pseudo/partial; intact_full 1,342 aa) have the same class, length and coordinates as the batch run; the two other genomes have no locus, as in the batch run.

One difference found and fixed after this test: the GenBank file of GCA_030582555.1 contains no CDS features. The command line ignored it (state `not_checked`); the batch run reported `no_CDS_annotated`, as the catalogue does. Since the fix, an annotation without CDS features is kept and a warning is written. The 5-genome test was not repeated after the fix; the test suite (section 1) was. Run times in this test (about 580 s per genome) reflect a busy computer and are not representative.

## 7. Known limits found by the tests

- Domain envelopes differ slightly from HMMER; one 73-aa fragment of 775 is missed.
- Fungal lineage II proteins are not found by the profiles or by the exploratory mode.
- Fungal genes interrupted by introns or split between contigs appear as fragments.
