# Fungal INA model (Mortierellaceae lineage I): profiles and rules

File: `fungal_model.json` (INA finder 1.1.0). Validation: `validation/validation_fungal_summary.md`.

## Profiles

| profile | role | positions | sequences | built from |
|---|---|---|---|---|
| INAf_N | N-terminal region | 277 | 10 (4 species) | lineage I proteins from the signal-peptide Met to the first GYGST repeat; *Mortierella alpina* 7, *M. antarctica*, *Podila clonocystis*, *Linnemannia amoeboidea* |
| INAf_Ccap | C-terminal cap | 45 | 14 (4 species) | last 45 aa of the 14 lineage I proteins with the Cys-rich cap |
| INA_C60 | C-terminal cap | 60 | 42 | the bacterial cap profile, for the one fungal INA with a bacterial-type cap |

Build: MAFFT L-INS-i alignment, `hmmbuild` (HMMER 3.4), default weighting. Source sequences: `fungal_build/fungal_N_region.faa`, `fungal_build/fungal_Ccap45.faa`; metadata `fungal_build/lineageI_meta.json`.

## Rules for a canonical match

| profile | alignment starts at profile position ≤ | reaches position ≥ | i-Evalue ≤ | position in the protein |
|---|---|---|---|---|
| INAf_N | 10 | 250 | 1e-10 | starts within the first 400 aa |
| INAf_Ccap | 5 | 40 | 1e-10 | ends within the last 10 aa |
| INA_C60 | 15 | 50 | 1e-15 | ends within the last 20 aa |

General rules are those of the bacterial model (16-aa repeat identity ≥ 0.45, intact full length ≥ 800 aa, intact short 500 to 799 aa, loci = ORFs on one strand within 8 kb), except `met_fallback = orf_start`.

## Why these values

1. **INAf_N must reach position 250.** Positions 1 to about 156 (signal peptide and the region after it) also match non-INA fungal proteins in INA-negative genomes, down to E = 1e-21 (real HMMER, `fungal_build/hs/`). A full INA N-region scores about 1e-160 over positions 1 to 277. Requiring coverage to position 250 keeps the INA-specific part in every canonical match.
2. **The last part of INAf_N (about positions 212 to 277) also matches single repeat units** (i-Evalue down to about 1e-20 in the JavaScript engine). These matches start far from position 10 and are never canonical; they only add to the count of weak matches.
3. **INAf_Ccap.** True caps score 1e-21 to 1e-34; the best non-cap match in 7 genomes was E ≈ 200.
4. **INA_C60 in fungi.** The bacterial cap profile matches fungal repeat units: 35 non-terminal matches with hmm_to ≥ 50 in 7 genomes, best i-Evalue 5.2e-9 (real HMMER). The only known fungal bacterial-type cap (locus0026, *M. claussenii*) scores 2.1e-26 over positions 1 to 52. The threshold 1e-15 sits between the two, and the cap must end within the last 20 aa.
5. **Start codon.** Fungal assemblies are fragmented and fungal genes can contain introns. When no Met lies near the first match, the whole ORF is kept and flagged `ORF_start_no_Met`, instead of starting at a Met far downstream.

## What the model cannot find

- **Lineage II** proteins: no conserved terminal domains. Six of the seven of at least 300 aa have 4-aa rather than 16-aa repeats (see `INA_EVOLUTION/ADDENDUM_lineageII_period_check_20261021.md`).
- **INA outside the Mortierellaceae.** All training sequences come from four species of one family. Leave-one-species-out tests transfer between these related species only.
- **Introns.** Six-frame translation does not splice. An INA gene with an intron in the N-region or the cap appears as fragments or as an interrupted locus.
