<p align="center"><img src="branding/INpro_FINDER_logo_readme.png" alt="INpro FINDER" width="420"></p>

# INpro FINDER 1.1.0

Annotation-independent search for ice-nucleation protein (INA) genes in genome sequences (formerly *INA finder*).

INpro FINDER translates every genome in all six reading frames and searches every open reading frame with profile hidden Markov models of the conserved INA terminal domains. It does not rely on the genome annotation, so it also finds INA genes that are annotated as "hypothetical protein", split into pieces or not annotated at all. For every locus it reports the structural state (intact, interrupted, fragment), the predicted protein, the gene with its ±20 kb neighbourhood and, if an annotation is supplied, how the locus is annotated.

**For users:** download `inafinder.html` from the latest release, open it in a browser and drop in genome files. Start with `docs/USER_GUIDE.md`.

| Folder | Contents |
|---|---|
| `dist/inafinder.html` | the browser tool: one self-contained file, works offline, nothing is uploaded |
| `dist/inafinder.py` | the command-line tool as one self-contained file (Python ≥ 3.8 and Node.js ≥ 18); it can also write the browser tool |
| `docs/` | user guide, test report, script that builds the PDF manual |
| `models/` | profile HMMs and model definitions: bacterial INA (`INA_N150`, `INA_C60`) and fungal INA of Mortierellaceae (`INAf_N`, `INAf_Ccap`) |
| `src/` | source of the search engine (`inafinder-core.js`), output formats (`inafinder-io.js`) and browser interface (`ui/`) |
| `cli/` | command-line source and the scripts that assemble `dist/` |
| `tests/` | test scripts; the reference genomes they use are not part of this repository, missing inputs are skipped |
| `branding/` | INpro FINDER logo (SVG/PNG) |

## Quick start

Browser: open `dist/inafinder.html`, drop genome FASTA files (and optionally their GFF3/GenBank annotation), choose *Bacterial INA*, *Fungal INA* or *Custom profiles*, press *Scan genomes*.

Command line:

```
python3 dist/inafinder.py scan genome1.fna genome2.fna.gz -o results/
python3 dist/inafinder.py scan *.fna --model fungal -o results_fungi/
python3 dist/inafinder.py viewer -o inafinder.html
```

Rebuild `dist/` after editing `src/` or `models/`: `python3 cli/build_html.py && python3 cli/make_single.py`.

The search follows HMMER 3 closely but is not HMMER; for publication, confirm key loci with `hmmsearch` using the same profiles.

## License

INpro FINDER is released under the [MIT License](LICENSE).

The bundled [Open Sauce Sans](https://github.com/marcologous/Open-Sauce-Fonts) typeface (embedded in `dist/inafinder.html`) keeps its own license (SIL Open Font License 1.1). The protein drawn in the logo is the AlphaFold model of InaZ ([AF-P06620-F1](https://alphafold.ebi.ac.uk/entry/P06620), CC BY 4.0).
