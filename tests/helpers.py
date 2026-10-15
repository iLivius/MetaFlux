"""Shared helpers of the MetaFlux test suite.

HOW TO RUN THE TESTS
    pytest tests            # from the repository root; see CONTRIBUTING.md
  The suite needs pytest and PyYAML; the OTU-table tests also need biom-format and
  numpy (the two packages of workflow/envs/biom.yaml) and are skipped without them.
  No conda, no Snakemake and no reference database is needed: every test runs on
  tiny inputs written on the spot, in a temporary folder pytest provides.

WHAT IS TESTED, AND HOW THE TESTS REACH THE CODE
  MetaFlux's Python lives in two kinds of file, and neither can simply be imported:
  the rule files (workflow/rules/**/*.smk) read the run config as soon as they are
  loaded, and the scripts (workflow/scripts/**/*.py) expect the `snakemake` object
  that Snakemake hands a running job. The three helpers below get around that
  without changing a line of the pipeline:

  load_functions()  pulls named functions and constants out of a file by parsing it
                    (Python's ast module reads the file as code without running it)
                    and runs only those definitions. Used for the start-up checks in
                    00_common.smk, which is plain Python, and for the parsers inside
                    the read-tracking scripts.
  run_script()      runs a whole script with a stand-in `snakemake` object that
                    carries the input, output, params and log a rule would pass.
                    The script then does exactly what it does in a run, on the tiny
                    files the test wrote. A script that stops the run calls
                    sys.exit(message); a test catches that as SystemExit and reads
                    the message.
  inline_program()  cuts a Python program out of a rule's shell command (the taxon
                    lookup of the read extraction lives there), so it can be run
                    with the Python interpreter as the rule runs it.

  marker_pack() reads a marker pack (workflow/markers/<name>.yaml) so that tests
  check tokens and references against the real rank models and value rules.
"""
import ast
import gzip
import json
import re
import struct
import sys
from collections import Counter
from pathlib import Path
from types import SimpleNamespace

import yaml

# Where things are, relative to this file: tests/ sits next to workflow/.
REPO_ROOT   = Path(__file__).resolve().parent.parent
WORKFLOW    = REPO_ROOT / "workflow"
COMMON_SMK  = WORKFLOW / "rules" / "shared" / "00_common.smk"
EXTRACT_SMK = WORKFLOW / "rules" / "shotgun" / "55_extract_taxa.smk"
SCRIPTS     = WORKFLOW / "scripts"
MARKERS     = WORKFLOW / "markers"
TEST_DIR    = REPO_ROOT / ".test"

# The standard-library names the extracted functions use. A function that needs
# another module fails with a clear NameError, which says what to add here.
STANDARD_NAMES = {
    "re": re, "struct": struct, "sys": sys, "json": json, "gzip": gzip,
    "Path": Path, "Counter": Counter,
}


def source_segments(path, names):
    """The source text of the top-level functions and constants called `names` in
    the Python file at `path`, as {name: text}, read with the ast module (the file
    is parsed, not run). Stops with a clear message when a name is not found."""
    source = Path(path).read_text()
    tree = ast.parse(source)
    found = {}
    for node in tree.body:
        if isinstance(node, ast.FunctionDef) and node.name in names:
            found[node.name] = ast.get_source_segment(source, node)
        elif isinstance(node, ast.Assign):
            for target in node.targets:
                if isinstance(target, ast.Name) and target.id in names:
                    found[target.id] = ast.get_source_segment(source, node)
    missing = [name for name in names if name not in found]
    if missing:
        raise AssertionError(f"{path}: no top-level definition of {missing}")
    return found


def load_functions(path, names):
    """Run the definitions of `names` from the file at `path` (functions and the
    constants they use) and return them as attributes of one object, e.g.
        checks = load_functions(COMMON_SMK, ["_value_passes"])
        checks._value_passes("Bacillati", {"ends_with": ["ati"]})
    Only these definitions run; the rest of the file (which would read a run config)
    is never executed."""
    segments = source_segments(path, names)
    namespace = dict(STANDARD_NAMES)
    # Definitions run in file order, so a constant is defined before the function
    # that uses it, exactly as when the file itself is loaded.
    ordered = []
    source = Path(path).read_text()
    for name, text in segments.items():
        ordered.append((source.index(text), name, text))
    ordered.sort()
    exec("\n\n".join(text for position, name, text in ordered), namespace)
    return SimpleNamespace(**{name: namespace[name] for name in names})


def fake_snakemake(input=None, output=None, params=None, log=None):
    """A stand-in for the `snakemake` object a running job gets: `sm.input.reports`,
    `sm.output.tsv`, `sm.params.samples`, `sm.log[0]` and so on. Each argument is a
    dict of the names the rule declares; `log` is a list of paths."""
    return SimpleNamespace(
        input=SimpleNamespace(**(input or {})),
        output=SimpleNamespace(**(output or {})),
        params=SimpleNamespace(**(params or {})),
        log=list(log or []),
    )


def run_script(script_path, snakemake_object):
    """Run a MetaFlux script top to bottom as Snakemake would, with `snakemake_object`
    standing in for the real one. Returns the script's namespace (its functions and
    variables), so a test can also call a function the script defines. A script that
    stops the run raises SystemExit, which the test catches."""
    script_path = Path(script_path)
    code = compile(script_path.read_text(), str(script_path), "exec")
    namespace = {"__name__": "__main__", "snakemake": snakemake_object}
    exec(code, namespace)
    return namespace


def inline_program(smk_path, start_marker, end_marker):
    """The text between `start_marker` and `end_marker` in a rule file: a Python
    program that a rule's shell command passes to `python -c`. The program starts
    right after the marker (the rest of that line) and ends where the end marker
    begins. Both markers must occur exactly once."""
    text = Path(smk_path).read_text()
    if text.count(start_marker) != 1 or text.count(end_marker) != 1:
        raise AssertionError(f"{smk_path}: the markers must occur exactly once")
    start = text.index(start_marker) + len(start_marker)
    end = text.index(end_marker, start)
    return text[start:end]


def marker_pack(name):
    """The shipped marker pack workflow/markers/<name>.yaml as a dict."""
    with open(MARKERS / f"{name}.yaml") as handle:
        return yaml.safe_load(handle)


def write_text(path, text):
    """Write `text` to `path`, creating the folders above it; returns the path as a
    string, which is what a rule hands a script."""
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)
    return str(path)
