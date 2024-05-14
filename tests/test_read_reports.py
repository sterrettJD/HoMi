import pytest
from src.rule_utils import read_reports as rr


def test_read_hostile():
    report = rr.read_hostile_report("tests/test_data/hostile_report.report")
    assert report["version"] == "1.1.0"


@pytest.fixture
def report():
    return rr.read_hostile_report("tests/test_data/hostile_report.report")


def test_hostile_reads_in(report):
    assert rr.get_reads_into_hostile(report) == 177014


def test_hostile_reads_nonhost(report):
    assert rr.get_nonhost_reads(report) == 59944


def test_hostile_percent_host(report):
    assert rr.get_percent_host(report) == 0.66136
