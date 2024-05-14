import pytest
from src.rule_utils import read_reports as rr


def test_read_hostile():
    report = rr.read_hostile_report("tests/test_data/hostile_report.report")
    assert report["reads_out"] == 59944