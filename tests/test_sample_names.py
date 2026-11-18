"""Sample names that R or pandas would read back as numbers.

_names_read_as_numbers() (00_common.smk) is the start-up check behind two refusals:
in shotgun mode with sylph (sylph-tax reads sample names with pandas) and in amplicon
mode (MetaFlux's rdp taxonomy step reads them with data.table's fread, and so do the
downstream readers of MetaFlux's tables). A name is refused when such a reader would
turn it into a number, a missing value or true/false. It is loaded from 00_common.smk
with load_functions (see helpers.py), together with the helpers and the constant it
uses; which mode it runs in is decided by two lines of 00_common.smk that a dry run
checks, not these tests.
"""
from helpers import COMMON_SMK, load_functions

names = load_functions(COMMON_SMK, ["_PANDAS_NA", "_pandas_reads_as_number",
                                    "_pandas_rewrites", "_names_read_as_numbers"])


def test_run_accessions_and_names_with_letters_pass():
    safe = ["SRR10043840", "S001", "soil01", "Pea-Hel-Loc7-Soil1", "12_3", "2019_05", "1_000"]
    assert names._names_read_as_numbers(safe) == []


def test_digit_only_names_are_refused():
    # '01' and '1' would become one sample, '002' would become 2.
    assert names._names_read_as_numbers(["01", "1", "002", "S2"]) == ["01", "1", "002"]


def test_other_number_forms_are_refused():
    for name in ["+5", "-3", "5.", ".5", "1e5", "2E-3", "Infinity", "-inf", "nan"]:
        assert names._names_read_as_numbers([name]) == [name], name


def test_missing_value_and_true_false_words_are_refused():
    for name in ["NA", "NULL", "None", "null", "n/a", "true", "FALSE", "True"]:
        assert names._names_read_as_numbers([name]) == [name], name


def test_the_order_of_the_names_is_kept():
    assert names._names_read_as_numbers(["10", "A", "2", "B", "NA"]) == ["10", "2", "NA"]
