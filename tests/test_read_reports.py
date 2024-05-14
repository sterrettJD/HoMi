import pytest
import pandas as pd
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


@pytest.fixture
def metadata():
    d = {"Sample": ["hostile_report"]}
    return pd.DataFrame(d)


def test_create_df_from_hostile_reports(metadata):
    actual = rr.create_df_from_hostile_reports(metadata, 
                                               hostile_directory="tests/test_data")
    expected = pd.DataFrame({"Reads passing QC": [177014],
                             "Nonhost reads": [59944],
                             "Percent host": [0.66136]},
                            index=["hostile_report"])
    equal = actual.values == expected.values
    assert equal.all()
    