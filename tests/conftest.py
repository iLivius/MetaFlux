"""pytest reads this file before the tests. It only makes sure that the helper
modules next to the tests (helpers.py, make_test_kraken2_db.py) can be imported
whatever folder pytest is started from; the helpers themselves are in helpers.py."""
import sys
from pathlib import Path

TESTS_DIR = str(Path(__file__).resolve().parent)
if TESTS_DIR not in sys.path:
    sys.path.insert(0, TESTS_DIR)
