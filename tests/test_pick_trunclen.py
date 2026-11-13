"""The overlap check of rule pick_trunclen (workflow/scripts/amplicon/50a_pick_trunclen.py)
and the 16S length-probe setting it depends on (workflow/markers/16S.yaml).

What is tested, and why it matters:
  - When no pair of cuts lets R1 and R2 overlap enough to merge, the run stops, under each
    of the three resolve_policy values (raise_trunc exit 3, relax_q exit 4, error exit 5).
    Each stop text gives the numbers behind the check (how far the reads reach after
    primer removal, R1 + R2; the expected amplicon length plus amplicon.min_overlap) and
    names the settings that feed it: amplicon.expected_length, and with auto the probe
    statistic and its exact key, amplicon.probe_length_stat.<key>. Up to v2.6.1 the text
    gave only a deficit and two headrooms, and users could not tell that a 582 bp probe
    statistic, not their reads, was at fault.
  - A length typed in amplicon.expected_length is named as such; a gyrB run names the
    key it really reads (16S) and says that its reference is already cut to the amplicon.
  - Under relax_q and error, a hint names raise_trunc when the reads would reach far
    enough; a lower probe statistic is offered only when one exists down to the median.
  - A run that passes the check writes the same trunclen.json as before.
  - The 16S marker pack leaves products with a second forward-primer site out of the
    probe statistics for both SILVA releases (probe_exclude_second_site at its top level,
    not overridden by a release block); 18S does not.

  - What the cut is based on: the quality values of the primer-trimmed reads are named
    (binned or continuous), with what q_threshold means on them and a warning when it
    equals one of the binned values or lies above the top one; each quality cut is said
    to come from quality or from the read-length limit.
  - A short off-target product (reads under half the expected amplicon length) no
    longer sets the read-length limit, unless it is most of the sample; the share left
    out is logged per sample. Without the change, the run below stops at the overlap.
  - method: simulate: the longest cut at which a read still passes filterAndTrim
    (truncQ, maxN, maxEE, maxLen), the table of pairs kept at every cut pair, the
    choice of the cut pair, read pairs sampled over the whole file, the same
    trunclen.json keys, the quality rule logged next to it, and ITS left to the
    quality rule.

The script is run whole with a stand-in snakemake object (run_script in helpers.py) on
Falco reports written to a temporary folder: one per sample and read direction, with a
quality value per position and every read of the same length, and on small FASTQ files
of primer-trimmed reads. Single functions are tested on the script loaded without
running it (load_picker).
"""
import json

from helpers import SCRIPTS, fake_snakemake, marker_pack, run_script, write_text

PICK_SCRIPT = SCRIPTS / "amplicon" / "50a_pick_trunclen.py"

# A probe result like the one SILVA 144 gave for 341F/805R before the 16S pack left out
# products with a second forward-primer site (q3 427, p95 582).
PROBE_341F_805R = {"amplicon_type": "16S", "probe_mode": "pcr", "n_amplicons": 607963,
                   "min": 32, "q1": 403, "median": 422, "q3": 427, "p95": 582,
                   "p99": 640, "max": 1818}


def falco_report(length, good_until, lengths=None, total=1000):
    """A minimal Falco fastqc_data.txt: lower quartile 30 up to `good_until`, 10 after it,
    and every read `length` bases long, or the read lengths in `lengths`
    ({length: number of reads}); `total` is its "Total Sequences" line."""
    rows = [">>Basic Statistics\tpass", "#Measure\tValue", f"Total Sequences\t{total}",
            ">>END_MODULE",
            ">>Per base sequence quality\tpass",
            "#Base\tMean\tMedian\tLower Quartile\tUpper Quartile\t10th Percentile\t90th Percentile"]
    for pos in range(1, length + 1):
        q1 = 30.0 if pos <= good_until else 10.0
        rows.append(f"{pos}\t{q1}\t{q1}\t{q1}\t{q1}\t{q1}\t{q1}")
    rows += [">>END_MODULE", ">>Sequence Length Distribution\tpass", "#Length\tCount"]
    for read_length, count in sorted((lengths or {length: 1000.0}).items()):
        rows.append(f"{read_length}\t{float(count)}")
    rows += [">>END_MODULE"]
    return "\n".join(rows) + "\n"


def fastq_text(reads):
    """FASTQ text of `reads`, a list of (bases, list of Phred values), Phred+33."""
    lines = []
    for index, (bases, phreds) in enumerate(reads):
        lines += [f"@read{index}", bases, "+", "".join(chr(q + 33) for q in phreds)]
    return "\n".join(lines) + "\n"


def uniform_reads(n, length, phred=30):
    """n reads of `length` A bases, all at the same quality."""
    return [("A" * length, [phred] * length) for _ in range(n)]


def load_picker():
    """The picker's functions, from the script loaded WITHOUT running main()."""
    code = compile(PICK_SCRIPT.read_text(), str(PICK_SCRIPT), "exec")
    namespace = {"__name__": "pick_trunclen"}
    exec(code, namespace)
    return namespace


def run_pick(tmp_path, r1=(250, 250), r2=(250, 250), policy="raise_trunc",
             expected="auto", stat="p95", key="16S", amp_type="16S", probe=None,
             samples=None, method="quality", q_threshold=20, extra_overlap=20):
    """Write each sample's two Falco reports and two FASTQ files (and the probe JSON for
    auto), run the picker, and return (exit code or None, log text, trunclen.json path).

    `samples` is a list of dicts, one per sample, with any of: "r1"/"r2" (length,
    good_until) for the Falco reports, "lengths_r1"/"lengths_r2" (the read-length
    distribution), "reads_r1"/"reads_r2" (the FASTQ reads, see fastq_text). Left out,
    one sample is made from `r1` and `r2`, with 50 reads of that length at Q30."""
    if samples is None:
        samples = [{}]
    falco_r1, falco_r2, fastq_r1, fastq_r2 = [], [], [], []
    for index, sample in enumerate(samples, start=1):
        name = f"S{index}"
        s_r1, s_r2 = sample.get("r1", r1), sample.get("r2", r2)
        reads_r1 = sample.get("reads_r1") or uniform_reads(50, s_r1[0])
        reads_r2 = sample.get("reads_r2") or uniform_reads(50, s_r2[0])
        falco_r1.append(write_text(tmp_path / "falco" / f"{name}_R1_stripped" / "fastqc_data.txt",
                                   falco_report(*s_r1, lengths=sample.get("lengths_r1"),
                                                total=len(reads_r1))))
        falco_r2.append(write_text(tmp_path / "falco" / f"{name}_R2_stripped" / "fastqc_data.txt",
                                   falco_report(*s_r2, lengths=sample.get("lengths_r2"),
                                                total=len(reads_r2))))
        fastq_r1.append(write_text(tmp_path / "03.stripped" / f"{name}_R1.fastq", fastq_text(reads_r1)))
        fastq_r2.append(write_text(tmp_path / "03.stripped" / f"{name}_R2.fastq", fastq_text(reads_r2)))
    probe_json = []
    if expected == "auto":
        probe_json = [write_text(tmp_path / "probe.json", json.dumps(probe or PROBE_341F_805R))]
    out_json = tmp_path / "stats" / "trunclen.json"
    log_path = tmp_path / "logs" / "pick_trunclen.log"
    sm = fake_snakemake(
        input={"falco_r1": falco_r1, "falco_r2": falco_r2, "probe_json": probe_json,
               "stripped_r1": fastq_r1, "stripped_r2": fastq_r2},
        output={"json": str(out_json)},
        params={"amplicon_type": amp_type, "expected_length": expected, "min_overlap": 12,
                "mode": "auto", "q_threshold": q_threshold, "q_floor": 15,
                "resolve_policy": policy, "manual_r1": None, "manual_r2": None,
                "probe_stat": stat, "probe_stat_key": key, "min_read_coverage_pct": 95,
                "method": method, "simulate_extra_overlap": extra_overlap,
                "trunc_q": 2, "max_ee": [2, 5], "min_len": 20, "max_len": None},
        log=[str(log_path)],
    )
    code = None
    try:
        run_script(PICK_SCRIPT, sm)
    except SystemExit as stop:
        code = stop.code
    return code, log_path.read_text(), out_json


def flat(text):
    """The text with every run of spaces and line ends made one space."""
    return " ".join(text.split())


def test_raise_trunc_stop_names_the_numbers_and_the_settings(tmp_path):
    code, log, out_json = run_pick(tmp_path)
    text = flat(log)
    assert code == 3
    assert not out_json.exists()
    assert "R1 and R2 cannot be cut so that the pairs still merge" in text
    assert "the reads reach 500 bp, 94 bp short of the 594 bp that merging needs" in text
    assert "reach at most R1 250 + R2 250 = 500 bp" in text
    assert "amplicon.trunc_len.min_read_coverage_pct" in text
    assert "582 bp (expected amplicon length) + 12 bp of overlap (amplicon.min_overlap) = 594 bp" in text
    assert "the p95 of the lengths of the 607,963 in-silico PCR products" in text
    assert "amplicon.expected_length: auto; amplicon.probe_length_stat.16S: p95" in text
    assert "set amplicon.expected_length to its length in bp" in text
    assert "under amplicon.probe_length_stat.16S (q3: 427 bp, median: 422 bp)" in text
    assert "no truncation setting makes the pairs merge" in text
    # the old wording is gone
    assert "headroom" not in text


def test_relax_q_stop_says_where_the_threshold_ended(tmp_path):
    code, log, _ = run_pick(tmp_path, policy="relax_q")
    text = flat(log)
    assert code == 4
    assert "lowered to Q15 (amplicon.trunc_len.q_floor: 15; resolve_policy: relax_q)" in text
    assert "the cuts reach R1 250 + R2 250 = 500 bp, 94 bp short of the 594 bp" in text
    assert "amplicon.probe_length_stat.16S: p95" in text


def test_error_stop_points_to_raise_trunc_when_the_reads_reach_far_enough(tmp_path):
    # The reads reach 300 + 300 = 600 bp, but quality drops at 200: the quality cuts
    # give 400 bp, 39 short of q3 427 + 12; raise_trunc would pass.
    code, log, _ = run_pick(tmp_path, r1=(300, 200), r2=(300, 200), policy="error", stat="q3")
    text = flat(log)
    assert code == 5
    assert "The quality cuts at Q20 (amplicon.trunc_len.q_threshold) reach R1 200 + R2 200 = 400 bp" in text
    assert "39 bp short of the 439 bp that merging needs, and resolve_policy: error stops here" in text
    assert "resolve_policy: raise_trunc (the default) extends the cuts" in text
    # q3 is the statistic in use, so only the median is offered as a lower one
    assert "(median: 422 bp)" in text
    assert "q3: 427 bp" not in text
    # the reads are long enough, so the text does not say that nothing helps
    assert "no truncation setting" not in text


def test_a_length_of_your_own_is_named_as_such(tmp_path):
    code, log, _ = run_pick(tmp_path, expected=600)
    text = flat(log)
    assert code == 3
    assert "600 bp, the number set in amplicon.expected_length" in text
    assert "amplicon.probe_length_stat is not used" in text
    assert "probe_length_stat.16S" not in text


def test_gyrb_names_the_key_it_reads_and_its_direct_reference(tmp_path):
    probe = {"amplicon_type": "gyrB", "probe_mode": "direct", "n_amplicons": 25000,
             "min": 200, "q1": 270, "median": 280, "q3": 285, "p95": 520, "p99": 600,
             "max": 900}
    code, log, _ = run_pick(tmp_path, amp_type="gyrB", key="16S", probe=probe)
    text = flat(log)
    assert code == 3
    assert "the 25,000 reference sequences, which are already cut to the amplicon" in text
    assert "amplicon.probe_length_stat.16S: p95" in text


def test_a_passing_run_writes_the_same_trunclen_json(tmp_path):
    # q3 427 + 12 = 439 bp; the reads give 250 + 250 = 500 bp at Q20: slack 61.
    code, _, out_json = run_pick(tmp_path, stat="q3")
    assert code in (None, 0)
    assert json.loads(out_json.read_text()) == {
        "r1": 250, "r2": 250, "mode": "auto", "amplicon_type": "16S",
        "q_threshold_requested": 20.0, "q_threshold_used": 20.0,
        "expected_length": 427, "expected_length_source": "auto:q3",
        "min_overlap": 12, "overlap_slack": 61, "resolved_via": "primary",
        "aggregation": "median_across_samples_per_bin",
    }


def test_16S_pack_leaves_second_site_products_out_for_both_releases():
    pack = marker_pack("16S")
    assert pack["probe_exclude_second_site"] is True
    for release, block in pack["releases"].items():
        # 00_common.smk copies a release block over the top level, so no block may
        # set the field to anything else
        assert "probe_exclude_second_site" not in (block or {}), release
    assert marker_pack("rpoB")["probe_exclude_second_site"] is True
    assert "probe_exclude_second_site" not in marker_pack("18S")


# ── What the cut is based on: quality values and the source of each cut ──────────

# Binned reads as a NextSeq 1000/2000 writes them with the 600-cycle kit: 9, 20, 34,
# and 2 at a no-call. Continuous reads: 2 and every value from 7 to 38, as a MiSeq.
BINNED_PHREDS = [34] * 200 + [20] * 30 + [9] * 18 + [2] * 2
CONTINUOUS_PHREDS = ([2] + list(range(7, 39))) * 7 + [30] * 19


def test_binned_values_are_named_and_the_edge_is_warned(tmp_path):
    reads = [("A" * 250, BINNED_PHREDS) for _ in range(20)]
    code, log, _ = run_pick(tmp_path, stat="q3", samples=[{"reads_r1": reads, "reads_r2": reads}])
    text = flat(log)
    assert code in (None, 0)
    assert "R1 quality values (first 10,000 reads of each of the 1 primer-trimmed files): 4 values (2, 9, 20, 34), binned" in text
    assert "R1 q_threshold 20 on these values: the cut falls where at least a quarter of the bases carry 9 or lower" in text
    assert "WARNING: R1 q_threshold 20 equals one of the quality values" in text
    # the Falco reports never fall below 20 within 250 bp: the cut is the length limit
    assert "R1 quality cut 250 bp is the read-length limit" in text


def test_continuous_values_give_no_edge_warning(tmp_path):
    reads = [("A" * 250, CONTINUOUS_PHREDS) for _ in range(20)]
    code, log, _ = run_pick(tmp_path, stat="q3", samples=[{"reads_r1": reads, "reads_r2": reads}])
    text = flat(log)
    assert "33 values (2-38), continuous" in text
    assert "the cut falls where at least a quarter of the bases are below Q20" in text
    assert "equals one of the quality values" not in text


def test_a_threshold_above_the_top_value_is_warned(tmp_path):
    reads = [("A" * 250, BINNED_PHREDS) for _ in range(20)]
    _, log, _ = run_pick(tmp_path, stat="q3", q_threshold=35,
                         samples=[{"reads_r1": reads, "reads_r2": reads}])
    assert "q_threshold 35 is above the highest quality value (34)" in flat(log)


def test_a_cut_set_by_quality_is_named_as_such(tmp_path):
    _, log, _ = run_pick(tmp_path, r1=(300, 200), r2=(300, 250), stat="q3")
    text = flat(log)
    assert "R1 quality cut 200 bp is set by quality: the median Q1 falls below 20 at position 201 (read-length limit 300 bp)" in text
    assert "R2 quality cut 250 bp is set by quality" in text


# ── A short off-target product does not set the read-length limit ────────────────

def test_a_short_offtarget_product_no_longer_stops_the_run(tmp_path):
    # Expected 431 bp + 12 = 443 bp. Sample S2 has 20 % of its reads in a 158 bp
    # product: over all reads its 95 % length is 158 bp, which capped the run at
    # 158 + 158 = 316 bp before. Half of 431 is 216, so the product is left out.
    product = {280: 800, 158: 200}
    samples = [{"r1": (280, 280), "r2": (280, 280)},
               {"r1": (280, 280), "r2": (280, 280), "lengths_r1": product, "lengths_r2": product}]
    code, log, out_json = run_pick(tmp_path, expected=431, samples=samples)
    text = flat(log)
    assert code in (None, 0)
    assert "not counting reads shorter than 216 bp" in text
    assert "S2_R1_stripped: 280 bp (20.00% of its reads are shorter than 216 bp and were left out)" in text
    assert "S1_R1_stripped: 280 bp (no read shorter than 216 bp)" in text
    assert json.loads(out_json.read_text())["r1"] == 280


def test_a_sample_made_mostly_of_short_reads_keeps_the_old_limit(tmp_path):
    # 60 % of S1's reads are under 216 bp: that is the sample, not an off-target
    # product, so its limit counts every read (158 bp) and the overlap check stops.
    product = {280: 400, 158: 600}
    samples = [{"r1": (280, 280), "r2": (280, 280), "lengths_r1": product, "lengths_r2": product}]
    code, log, _ = run_pick(tmp_path, expected=431, samples=samples)
    text = flat(log)
    assert code == 3
    assert "WARNING: S1_R1_stripped: 60.0% of its reads are shorter than 216 bp, more than half of the sample" in text


def test_no_short_reads_leave_the_limit_exactly_as_before():
    pick = load_picker()
    counts = [(271, 49512.0), (272, 135.0), (275, 37199.0), (281, 91042.0), (289, 1.0)]
    old = pick["read_coverage_cap"](counts, 95)
    assert pick["read_coverage_cap_with_floor"](counts, 95, 216) == (old, 0.0, False)
    assert pick["read_coverage_cap_with_floor"](counts, 95, None) == (old, 0.0, False)


# ── method: simulate ─────────────────────────────────────────────────────────────

def test_longest_passing_cut_follows_filterandtrim():
    cut = load_picker()["longest_passing_cut"]
    good = [40] * 100
    assert cut(b"A" * 100, bytes(q + 33 for q in good), 2, 2.0) == 100
    # truncQ: cut before the first base at Q2 or lower (base 51)
    low = good[:50] + [2] + good[51:]
    assert cut(b"A" * 100, bytes(q + 33 for q in low), 2, 2.0) == 50
    # maxN 0: an N at base 31 leaves 30 bases
    assert cut(b"A" * 30 + b"N" + b"A" * 69, bytes(q + 33 for q in good), 2, 2.0) == 30
    # maxEE: Q10 bases carry 0.1 expected errors each, so 20 of them stay within 2.05
    assert cut(b"A" * 100, bytes([10 + 33]) * 100, 2, 2.05) == 20
    # maxLen is applied to the whole read first
    assert cut(b"A" * 100, bytes(q + 33 for q in good), 2, 2.0, max_len=90) == 0


def test_the_cut_pair_keeps_the_most_pairs_that_still_merge():
    pick = load_picker()
    # (M_R1, M_R2) of four pairs; cuts must reach 300 bp together
    pairs = [(200, 150), (200, 120), (150, 200), (180, 180)]
    kept = pick["survival_table"](pairs, 200)
    assert pick["kept_at"](kept, 150, 120) == 4
    assert pick["kept_at"](kept, 180, 150) == 2
    assert pick["kept_at"](kept, 201, 1) == 0
    # 180/120 keeps 3 pairs (all but (150, 200)), the most any pair reaching 300 bp can
    assert pick["best_simulated_cuts"](kept, 300, 20) == (180, 120, 3)
    assert pick["best_simulated_cuts"](kept, 500, 20) is None
    # With a tolerance of one pair, any cut pair keeping 2 counts as good as 3, and the
    # longest of those is 180/150 (330 bp)
    assert pick["best_simulated_cuts"](kept, 300, 20, tolerance=1) == (180, 150, 2)


def test_pairs_are_sampled_over_the_whole_file(tmp_path):
    simulate_sample = load_picker()["simulate_sample"]
    # 1,000 pairs; only the second half of R2 has an N at base 11 (as a tile written late)
    good = [30] * 100
    r1 = [("A" * 100, good)] * 1000
    r2 = [("A" * 100, good)] * 500 + [("A" * 10 + "N" + "A" * 89, good)] * 500
    p1 = write_text(tmp_path / "r1.fastq", fastq_text(r1))
    p2 = write_text(tmp_path / "r2.fastq", fastq_text(r2))
    pairs = simulate_sample(p1, p2, 1000, 100, 2, 2.0, 5.0, None)
    assert len(pairs) == 100
    assert sum(1 for _m1, m2 in pairs if m2 == 10) == 50


def test_simulate_writes_the_same_keys_and_logs_the_quality_rule(tmp_path):
    # q3 427 + 12 = 439 bp, + 20 bp extra = 459 bp. R1 reads are good over 250 bp;
    # 30 % of R2 reads hit a Q2 base at 211, so they pass up to 210. The quality rule
    # (Falco reports clean to 250) cuts 250/250 and loses those 30 %; simulate keeps
    # them all with 250/210 (210 rather than 209: equal pairs kept, more overlap).
    good = [30] * 250
    bad = [30] * 210 + [2] + [30] * 39
    reads_r2 = [("A" * 250, good)] * 70 + [("A" * 250, bad)] * 30
    samples = [{"reads_r1": [("A" * 250, good)] * 100, "reads_r2": reads_r2}]
    code, log, out_json = run_pick(tmp_path, stat="q3", method="simulate", samples=samples)
    text = flat(log)
    assert code in (None, 0)
    result = json.loads(out_json.read_text())
    assert (result["r1"], result["r2"], result["resolved_via"]) == (250, 210, "simulate")
    quality_keys = {"r1", "r2", "mode", "amplicon_type", "q_threshold_requested",
                    "q_threshold_used", "expected_length", "expected_length_source",
                    "min_overlap", "overlap_slack", "resolved_via", "aggregation"}
    assert set(result) == quality_keys
    assert result["overlap_slack"] == 21
    assert "Quality rule for comparison: R1=250, R2=250 via primary, keeps 70.00% of the sampled pairs" in text
    assert "Simulated cut: R1=250, R2=210, keeps 100.00% of the sampled pairs" in text
    assert "S1: 100 pairs sampled, 100.0% kept at the simulated cut, 70.0% at the quality rule's cut" in text


def test_simulate_leaves_its_to_the_quality_rule(tmp_path):
    code, log, out_json = run_pick(tmp_path, amp_type="ITS", expected=None, method="simulate")
    assert code in (None, 0)
    assert "method=simulate does not apply to ITS" in log
    assert json.loads(out_json.read_text())["resolved_via"] == "primary"


def test_an_unknown_method_stops_with_a_readable_message(tmp_path):
    code, log, _ = run_pick(tmp_path, stat="q3", method="figaro")
    assert code == 2
    assert "amplicon.trunc_len.method must be quality or simulate, got 'figaro'" in log
