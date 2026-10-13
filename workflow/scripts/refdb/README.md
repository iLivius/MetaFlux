# Reference-database preparation scripts

This folder holds the **database-specific transform logic** — one script per
reference that MetaFlux cannot use as downloaded. Keeping them here (rather than
inline in a Snakemake rule) means there is a single, obvious place to look when you
ask *"how is each reference turned into a classifier-ready file?"*

The rules that call these live in [`workflow/rules/amplicon/10_refdb.smk`](../../rules/amplicon/10_refdb.smk),
which also carries the full **reference preparation map** (which DB is a plain
download, which needs unpacking, which needs one of the scripts below).

| Script | Turns | Into | Used by (rule) |
|---|---|---|---|
| `silva_drop_kingdom.py` | SILVA 144 DADA2 trainset (7 ranks) | 6-rank DADA2 trainset + phylum → kingdom table | `convert_silva_train` |
| `silva_to_sintax.py` | SILVA 16S DADA2 trainset (6 ranks; stops on more) | VSEARCH SINTAX FASTA | `convert_silva_sintax` |
| `frogs_rpob_to_dada2.py` | FROGS rpoB release FASTA | DADA2 assignTaxonomy trainset | `convert_rpob_to_dada2` |
| `pr2_utax_domain_first.py` | PR2 18S UTAX FASTA (domain under `k:`, supergroup under `d:`) | VSEARCH SINTAX FASTA with the domain under `d:` and the supergroup under `k:`, the order VSEARCH reads; headers only | `convert_pr2_sintax` |

## Conventions

- **In-workflow, from canonical sources.** We fetch each database from its upstream
  home and transform it here, rather than re-hosting pre-baked copies. The tree is
  reproducible from a clean clone and tracks upstream when a URL is bumped.
- **Sequence lines are never altered** — only FASTA *headers* are reshaped, so the
  sequences an ASV is compared against are exactly the upstream ones.
- **One script = one database's quirk.** A transform that is generic (a plain
  DADA2↔SINTAX rank remap) should become a single rank-model-aware converter driven
  by the marker pack's `tax_levels` / `rank_letters` / `rank_prefixes` /
  `prefix_style`, not a copy per marker. `silva_to_sintax.py` is deliberately
  SILVA/16S-specific today (see its docstring); generalise it when a second marker
  needs SINTAX.
- **Rank models stay in the marker packs**, not here. These scripts fix *file
  format*; the number of ranks and the prefix style are marker facts that
  `workflow/markers/*.yaml` owns.
