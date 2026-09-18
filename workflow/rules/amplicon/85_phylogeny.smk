# Optional de novo 16S phylogeny — builds a tree from the run's own ASVs.
#
# WHAT THIS IS FOR. Some diversity measures are phylogenetic: Faith's PD adds up
# branch length across the tree spanned by the species present, and UniFrac compares
# two samples by how much of the tree they share. Both need a tree whose tips are the
# ASVs in your abundance table. Users have always built that tree by hand in R after
# a MetaFlux run; this module moves the build into the workflow, so it is versioned,
# resource-managed, seeded, and recorded like every other step.
#
# WHERE MetaFlux STOPS. When the Newick file is written. MetaFlux computes no
# UniFrac, no Faith's PD, no ordination, no PERMANOVA — you do that in R, on your own
# filtered feature set. That is a deliberate boundary, not a stage that has yet to be
# built, and it is why the tree is exported UNROOTED (see phylo_export below).
#
# WHAT IT IS NOT. This is an approximate, same-region marker-fragment phylogeny. A
# ~400 bp fragment of one gene carries limited deep signal, so it resolves close
# relatives reasonably and says little about deep relationships. Nothing here is a
# species tree.
#
# THE INPUT DECISION THAT MATTERS MOST. This module reads 6.taxonomy/, the ASV set
# AFTER the contaminant filter, and never 5.dada2/seqs_lenfilt.fasta. The filter
# (keep/discard: chloroplast, mitochondria, wrong domain) runs inside stage 80, so the
# stage-60 output still holds exactly the off-target sequences that show up in a de
# novo tree as long branches — the documented failure mode where a handful of
# divergent ASVs invent a sample separation that isn't real. Feeding the filtered set
# removes that by construction, and guarantees the tree's tips and the abundance table
# you pair it with in R are the same feature set.
#
# 16S ONLY, and enforced in 00_common.smk rather than here: for the other markers a
# de novo tree is not merely untested but positively contraindicated, each for its own
# published reason (see docs/amplicon/phylogeny.md).
#
# THE RULES, in data-flow order:
#   phylo_input   : 6.taxonomy/ tables      -> asv_16s.fasta + input.json   [CHECKPOINT]
#   phylo_align   : asv_16s.fasta           -> asv_16s.aln.fasta   (MAFFT)
#   phylo_tree    : the alignment           -> a backend-native tree
#                   (iqtree | fasttree | raxml-ng — exactly one rule defined, chosen
#                    at parse time, so the DAG always sees a single rule named
#                    phylo_tree; same idiom as target_extract and assign_taxonomy)
#   phylo_export  : tree + alignment        -> asv_16s.unrooted.nwk and
#                                              phylogeny.params.json, after the
#                                              correctness gates pass
#   phylo_qc      : everything above        -> stats/phylogeny/phylogeny_qc.json
#
# The whole file is inside `if PHYLO_ENABLED:`. With the feature off, these rules are
# not merely unused — they do not exist, so the DAG is byte-identical to a run before
# this module was written and `--sdm conda` never builds a phylogeny environment.

if PHYLO_ENABLED:

    # Where this stage writes. 7.phylogeny/ continues the N.stage/ numbering (taxonomy
    # is stage 6); stats/phylogeny/ holds the diagnostic report, matching every other
    # stage's habit of keeping deliverables and diagnostics apart.
    PHYLO_DIR    = OUT / "7.phylogeny"
    PHYLO_QC_DIR = OUT / "stats" / "phylogeny"

    # Minimum ASV count for a tree to mean anything. Below four tips there is only one
    # unrooted topology, so there is nothing to infer. Runs under this threshold are
    # skipped and the skip is recorded — see the checkpoint note on phylo_qc.
    PHYLO_MIN_ASVS = 4

    # Per-backend artifacts live in their own subdirectory, so switching backends
    # leaves the previous run's files visible for comparison rather than half-
    # overwriting them. phylo_export warns when it finds a directory from a backend
    # that is not the active one.
    PHYLO_BACKEND_DIR = PHYLO_DIR / PHYLO_BACKEND


    # ── Stage 85a: the eligible ASV set ────────────────────────────────────────
    # Pulls ASV_ID -> sequence out of taxon_seq_table.txt and writes a plain FASTA
    # whose headers are bare ASV_N — no taxonomy appended, ever. Two reasons: FastTree
    # truncates a name at the first space unless quoted, and both FastTree and RAxML-NG
    # reject the ":,()" that a lineage string is full of. Users join taxonomy back on
    # in R, keyed on the ASV ID.
    #
    # asv_table.txt is read as well, as the authority on which ASVs survived the
    # contaminant filter. The two tables must describe the same ASVs; if they do not,
    # something upstream is inconsistent and the script stops rather than quietly
    # building a tree of the intersection.
    #
    # WHY THIS IS A CHECKPOINT (Snakemake's term). Normally Snakemake plans the whole
    # run before executing anything. A checkpoint is a rule after which it stops and
    # re-plans, using what that rule actually wrote. It is needed here because whether
    # a tree gets built at all depends on the ASV count, which nobody knows until
    # taxonomy has run — see phylo_qc for how the two cases are handled. This is the
    # only checkpoint in MetaFlux.
    checkpoint phylo_input:
        input:
            taxon_seq_table = OUT / "6.taxonomy" / "taxon_seq_table.txt",
            asv_table       = OUT / "6.taxonomy" / "asv_table.txt",
        output:
            seqs       = PHYLO_DIR / "asv_16s.fasta",
            input_json = PHYLO_DIR / "input.json",
        params:
            min_asvs = PHYLO_MIN_ASVS,
        log:
            LOGS / "phylo_input.log",
        conda:
            "../../envs/python_utils.yaml"
        threads: lambda wc: threads_for("phylo_input")
        resources:
            mem_mb = lambda wc: mem_mb_for("phylo_input"),
        script:
            "../../scripts/amplicon/85a_phylo_input.py"


    # ── Stage 85b: alignment ───────────────────────────────────────────────────
    # MAFFT, fixed — the aligner is not user-selectable. Two strategies, chosen by how
    # many ASVs there are, at MAFFT's own documented ~200-sequence ceiling for its
    # accuracy methods:
    #   < 200  L-INS-i  (--localpair --maxiterate 1000 --threadit 0)
    #   >= 200 FFT-NS-2 (--retree 2 --maxiterate 0)
    # The count is only known at run time, so the choice is made inside the script and
    # written into the run record — that recording is what makes a data-dependent
    # choice acceptable here.
    #
    # Not used, deliberately: --auto (its thresholds are unpublished, live only in
    # MAFFT's source, and shift the algorithm under you when a dataset crosses a size
    # boundary), --adjustdirection (it renames flipped sequences with an _R_ prefix,
    # which would break ASV-ID matching against your abundance table; MetaFlux fixes
    # orientation upstream at primer trimming), and --thread -1 (takes every physical
    # core regardless of what Snakemake allocated).
    #
    # Output: the alignment, retained as a deliverable so it can be re-used or
    # inspected, plus aln_run.json, the provenance phylo_export folds into
    # phylogeny.params.json.
    rule phylo_align:
        input:
            seqs = PHYLO_DIR / "asv_16s.fasta",
        output:
            aln     = PHYLO_DIR / "asv_16s.aln.fasta",
            run_json = PHYLO_DIR / "aln_run.json",
        params:
            strategy   = PHYLO_ALIGNER_STRATEGY,
            extra_args = PHYLO_EXTRA_ARGS.get("mafft", ""),
            linsi_threshold = 200,
        log:
            LOGS / "phylo_align.log",
        conda:
            "../../envs/mafft.yaml"
        threads: lambda wc: threads_for("phylo_align")
        resources:
            mem_mb = lambda wc: mem_mb_for("phylo_align"),
        script:
            "../../scripts/amplicon/85b_phylo_align.py"


    # ── Stage 85c: tree inference ──────────────────────────────────────────────
    # Three backends, one rule name. Exactly one branch below is evaluated, decided at
    # parse time from phylogeny.backend, so the DAG always contains a single rule
    # called phylo_tree — the same pattern 70_extract.smk uses for Metaxa2/ITSx and
    # 80_taxonomy.smk for rdp/sintax. The branches are exhaustive: 00_common.smk has
    # already rejected any other value, so no fourth case can arrive here.
    #
    # Common to all three: the tree, and a tree_run.json recording the resolved model,
    # the thread count actually used, the tool version, the seed, and the exact command
    # line. Each script also clears its tool's leftovers before running. Snakemake
    # already deletes a rule's DECLARED outputs before re-running it; the manual clean
    # is for the undeclared files these tools scatter around their prefix, which are
    # what trigger their "this analysis already finished" refusals.
    #
    # Every backend is seeded from amplicon.seed — the same seed as DADA2 and the
    # taxonomy step. There is deliberately no second seed key.

    if PHYLO_BACKEND == "iqtree":
        # The default backend: full maximum likelihood, and the only one of the three
        # offering real model selection.
        #
        # -m defaults to GTR+F+G4 rather than IQ-TREE's own -m MFP. MFP re-runs
        # ModelFinder on every execution — 22 base models crossed with frequency and
        # rate-heterogeneity variants — which makes runtime depend on the data and
        # differ run to run. For a single short marker gene that search lands in the
        # GTR family essentially always, and IQ-TREE's own documentation recommends
        # selecting once and then pinning. Users who want the search set model: MFP;
        # the model it picks is written into the .iqtree report and recorded.
        #
        # Threads are deliberately modest (~4 by default). IQ-TREE parallelises across
        # alignment columns, so its efficiency depends on alignment LENGTH; at the
        # ~250-430 bp of a 16S amplicon the benefit collapses quickly and more threads
        # can be slower. -T AUTO is never used: it benchmarks the machine it happens
        # to be running on, under whatever load it happens to be under, so the same
        # input would behave differently on different hardware.
        #
        # Duplicate sequences need no special handling: IQ-TREE removes them, then
        # re-inserts them at the end, so every input label appears in the .treefile.
        # --keep-ident is therefore not passed.
        #
        # Declared outputs are the three files every run produces. The .log is left
        # undeclared on purpose — Snakemake deletes a failed job's declared outputs,
        # and that log is exactly what you need to find out why it failed; the script
        # tees it into this rule's own log as well. Also undeclared: .uniqueseq.phy
        # (written only when duplicates exist), .bionj/.mldist, and the extra files
        # -B/MFP produce. All of them are covered by the pre-run clean.
        rule phylo_tree:
            input:
                aln = PHYLO_DIR / "asv_16s.aln.fasta",
            output:
                tree     = PHYLO_BACKEND_DIR / "asv_16s.treefile",
                report   = PHYLO_BACKEND_DIR / "asv_16s.iqtree",
                ckp      = PHYLO_BACKEND_DIR / "asv_16s.ckp.gz",
                run_json = PHYLO_BACKEND_DIR / "tree_run.json",
            params:
                prefix     = lambda wc: str(PHYLO_BACKEND_DIR / "asv_16s"),
                model      = PHYLO_MODEL_RESOLVED,
                model_requested = PHYLO_MODEL,
                support    = PHYLO_SUPPORT,
                seed       = config["amplicon"]["seed"],
                extra_args = PHYLO_EXTRA_ARGS.get("iqtree", ""),
            log:
                LOGS / "phylo_tree.log",
            conda:
                "../../envs/iqtree.yaml"
            threads: lambda wc: threads_for("phylo_tree")
            resources:
                mem_mb = lambda wc: mem_mb_for("phylo_tree"),
            script:
                "../../scripts/amplicon/85c_phylo_tree_iqtree.py"

    elif PHYLO_BACKEND == "fasttree":
        # The fast approximate-ML backend, and the one QIIME 2 and nf-core/ampliseq
        # use by default — though both invoke it with no model flags at all, which
        # means their shipped 16S tree is Jukes-Cantor. MetaFlux passes -gtr
        # explicitly, because JC is FastTree's default and is not defensible for 16S.
        #
        # Single-threaded, always. FastTree's parallel build (FastTreeMP) is
        # documented as non-deterministic both run-to-run and against the serial
        # binary, and its speedup does not extend to the maximum-likelihood phase past
        # about three cores anyway. The serial binary costs ~93 s at 10,000 ASVs, so
        # there is nothing to gain by risking it. `threads: 1` is written here rather
        # than read from the config for that reason: not because a larger number would
        # change the result (it cannot — the serial binary ignores it and the script
        # pins OMP_NUM_THREADS=1), but so Snakemake does not reserve cores this job
        # can never use, and so the single-threaded guarantee holds no matter what the
        # config says. resources.threads.phylo_tree still applies to the other two
        # backends.
        #
        # -gamma rescales branch lengths by a single global factor. Worth knowing what
        # that means downstream: a global factor cancels in UniFrac, which is a ratio,
        # but not in Faith's PD, which is an absolute sum.
        #
        # -noml is never used: it is documented to return slightly negative branch
        # lengths, which breaks both metrics.
        rule phylo_tree:
            input:
                aln = PHYLO_DIR / "asv_16s.aln.fasta",
            output:
                tree     = PHYLO_BACKEND_DIR / "asv_16s.nwk",
                run_json = PHYLO_BACKEND_DIR / "tree_run.json",
            params:
                model      = PHYLO_MODEL_RESOLVED,
                model_requested = PHYLO_MODEL,
                support    = PHYLO_SUPPORT,
                seed       = config["amplicon"]["seed"],
                extra_args = PHYLO_EXTRA_ARGS.get("fasttree", ""),
            log:
                LOGS / "phylo_tree.log",
            conda:
                "../../envs/fasttree.yaml"
            threads: 1
            resources:
                mem_mb = lambda wc: mem_mb_for("phylo_tree"),
            script:
                "../../scripts/amplicon/85d_phylo_tree_fasttree.py"

    else:  # PHYLO_BACKEND == "raxml-ng" — the only remaining case (validated in 00_common.smk)
        # Full maximum likelihood, run in RAxML-NG's own two-step form:
        #   --parse  validates the alignment, writes a fast-loading binary copy, and
        #            reports how many threads it recommends and how much memory it
        #            expects to need
        #   --search does the actual inference, using that thread count
        #
        # The two steps are not an optimisation — they are how the rule avoids simply
        # crashing. RAxML-NG wants a minimum number of alignment patterns per thread
        # and enforces it: moderately over-allocated it warns, badly over-allocated it
        # TERMINATES ("Too few patterns per thread"). Measured on the 211-ASV 16S test
        # alignment (328 distinct patterns, recommendation 1 thread): 8 threads ran
        # with a warning, 16 and 32 terminated. A normal Snakemake allocation of 16
        # would therefore kill the rule outright. The script uses
        # min(allocated, recommended) and records both numbers.
        #
        # --workers is likewise an explicit integer, not `auto`: auto reads the
        # machine, and the same input must not behave differently on different
        # hardware.
        #
        # --seed is mandatory here in a way it is not elsewhere: RAxML-NG's default
        # seed is the wall clock, so an unseeded run is non-reproducible by design.
        #
        # --extra seq-dup-remove is never exposed. It drops duplicate tips, which
        # would leave the tree describing a different feature set than the abundance
        # table — the one thing that must never happen. The defaults keep every tip.
        rule phylo_tree:
            input:
                aln = PHYLO_DIR / "asv_16s.aln.fasta",
            output:
                tree     = PHYLO_BACKEND_DIR / "asv_16s.raxml.bestTree",
                run_json = PHYLO_BACKEND_DIR / "tree_run.json",
            params:
                prefix       = lambda wc: str(PHYLO_BACKEND_DIR / "asv_16s"),
                parse_prefix = lambda wc: str(PHYLO_BACKEND_DIR / "asv_16s_parse"),
                model        = PHYLO_MODEL_RESOLVED,
                model_requested = PHYLO_MODEL,
                seed         = config["amplicon"]["seed"],
                extra_args   = PHYLO_EXTRA_ARGS.get("raxml_ng", ""),
            log:
                LOGS / "phylo_tree.log",
            conda:
                "../../envs/raxml-ng.yaml"
            threads: lambda wc: threads_for("phylo_tree")
            resources:
                mem_mb = lambda wc: mem_mb_for("phylo_tree"),
            script:
                "../../scripts/amplicon/85e_phylo_tree_raxml.py"


    # ── Stage 85d: validate and export ─────────────────────────────────────────
    # Turns whichever backend ran into MetaFlux's canonical, backend-independent
    # deliverables — but only after the correctness gates pass. Nothing canonical is
    # written unless every tip is accounted for, the Newick parses, and every branch
    # length is finite and non-negative. A tree that fails any of those would be worse
    # than no tree, because it looks usable.
    #
    # THE TREE IS EXPORTED UNROOTED, AND ONLY UNROOTED — a decision, not an omission.
    # Faith's PD and UniFrac both depend on where the root is, and users prune the tree
    # in R before computing them: decontam against negative controls, abundance
    # filtering of rare ASVs. Any root MetaFlux placed on the full ASV set would stop
    # being valid the moment the first tip is pruned (a midpoint in particular moves).
    # So rooting belongs in R, after pruning, immediately before the metric — and no
    # rooted copy is offered, because a file that is only valid until the user does the
    # thing they will certainly do is a trap, not a convenience.
    #
    # phylogeny.params.json is the run's provenance record: resolved aligner strategy,
    # resolved model, effective thread counts, tool versions, seed, input checksums,
    # and any extra_args. It is what makes a data-dependent choice reproducible later.
    rule phylo_export:
        input:
            aln         = PHYLO_DIR / "asv_16s.aln.fasta",
            tree        = rules.phylo_tree.output.tree,
            input_json  = PHYLO_DIR / "input.json",
            aln_json    = PHYLO_DIR / "aln_run.json",
            tree_json   = rules.phylo_tree.output.run_json,
        output:
            unrooted = PHYLO_DIR / "asv_16s.unrooted.nwk",
            params_json = PHYLO_DIR / "phylogeny.params.json",
        params:
            backend      = PHYLO_BACKEND,
            phylo_dir    = lambda wc: str(PHYLO_DIR),
            extra_args   = PHYLO_EXTRA_ARGS,
        log:
            LOGS / "phylo_export.log",
        conda:
            "../../envs/phylo_export.yaml"
        threads: lambda wc: threads_for("phylo_export")
        resources:
            mem_mb = lambda wc: mem_mb_for("phylo_export"),
        script:
            "../../scripts/amplicon/85f_phylo_export.py"


    # ── Stage 85e: QC report, and the skip decision ────────────────────────────
    # This function is where the "too few ASVs" case is actually handled. Snakemake
    # calls it AFTER the phylo_input checkpoint has run, so unlike everything else in
    # this file it can look at a real number.
    #
    #   4 or more eligible ASVs -> ask for the alignment, tree and params file. Asking
    #                              is what causes them to be built: this is the single
    #                              thread by which the whole chain enters the run.
    #   fewer than 4            -> ask for nothing but the checkpoint's own output.
    #                              No tree is attempted and 85g records why. Three
    #                              sequences have only one possible unrooted topology,
    #                              so there is nothing to infer; failing the run over
    #                              it would be wrong, and inventing a tree would be
    #                              worse.
    def _phylo_qc_inputs(wildcards):
        checkpoint_output = checkpoints.phylo_input.get(**wildcards).output
        with open(checkpoint_output.input_json) as fh:
            n_eligible = json.load(fh)["n_eligible"]

        wanted = {"input_json": str(checkpoint_output.input_json)}
        if n_eligible >= PHYLO_MIN_ASVS:
            wanted["aln"]         = str(PHYLO_DIR / "asv_16s.aln.fasta")
            wanted["unrooted"]    = str(PHYLO_DIR / "asv_16s.unrooted.nwk")
            wanted["params_json"] = str(PHYLO_DIR / "phylogeny.params.json")
        return wanted


    # The run's phylogeny diagnostic, and the only phylogeny file `rule all` asks for
    # (see the note in _amplicon_targets, 00_common.smk). Reports tip accounting,
    # alignment shape, and — the single most useful check for this module's known
    # failure mode — a long-branch report: any tip whose root-to-tip distance stands
    # far above the median is flagged. That is what a contaminant or an unrelated
    # off-target sequence looks like in a de novo tree, and it is the thing most likely
    # to distort an unweighted UniFrac result. Flagged only, never removed: dropping
    # tips would desynchronise the tree from your abundance table.
    #
    # Written on every enabled run, including a skipped one, so the reason a tree is
    # missing is always recorded somewhere.
    rule phylo_qc:
        input:
            unpack(_phylo_qc_inputs),
        output:
            qc_json = PHYLO_QC_DIR / "phylogeny_qc.json",
        params:
            min_asvs             = PHYLO_MIN_ASVS,
            backend              = PHYLO_BACKEND,
            # Two different threshold rules, because the two measures behave
            # differently — see long_branch_report in 85g_phylo_qc.py for the measured
            # justification. Root-to-tip distances share a root and cluster tightly,
            # so a multiple of the median works (1.5x, roughly the ratio in the
            # published UniFrac-artifact case). Pendant edges cannot use the median:
            # real ASV sets carry many near-identical sequences whose pendant edges
            # are ~0, so the median is ~0 and any multiple of it flags a third of the
            # tree (measured: 81 of 211 tips on the 16S test set). The boxplot fence
            # Q3 + 3*IQR flags 7 there, all of them genuinely divergent singleton
            # lineages. Both rules only decide what the log calls out — the longest
            # five tips by each measure are always listed in the report.
            long_branch_multiple          = 1.5,
            long_branch_pendant_fence_iqr = 3.0,
        log:
            LOGS / "phylo_qc.log",
        conda:
            "../../envs/phylo_export.yaml"
        threads: lambda wc: threads_for("phylo_qc")
        resources:
            mem_mb = lambda wc: mem_mb_for("phylo_qc"),
        script:
            "../../scripts/amplicon/85g_phylo_qc.py"
