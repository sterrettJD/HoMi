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
    d = {"Sample": ["hostile_report", "hostile_report_1"]}
    return pd.DataFrame(d)


def test_create_df_from_hostile_reports(metadata):
    actual = rr.create_df_from_hostile_reports(metadata, 
                                               hostile_directory="tests/test_data")
    expected = pd.DataFrame({"Reads passing QC": [177014, 59900],
                             "Nonhost reads": [59944, 59900],
                             "Percent host": [0.66136, 0]},
                            index=["hostile_report", "hostile_report_1"])
    equal = actual.values == expected.values
    assert equal.all()


def test_get_unmapped_nonhost():
    unmapped = rr.get_unmapped_nonhost_from_humann("tests/test_data/all_genefamilies.tsv")
    expected = {"hostile_report": 4698.0, "hostile_report_1": 29696.0}
    assert unmapped == expected
    

def test_add_unmapped_nonhost_to_hostile():
    hostile = pd.DataFrame({"Reads passing QC": [177014, 59900],
                             "Nonhost reads": [59944, 59900],
                             "Percent host": [0.66136, 0]},
                            index=["hostile_report", "hostile_report_1"])
    # note that this is in "reverse" order just to make sure it is merging correctly 
    # independent of dict order
    unmapped_nonhost =  {"hostile_report_1": 29696.0, "hostile_report": 4698.0}
    
    actual = rr.add_unmapped_nonhost_to_hostile_reports(unmapped_nonhost, hostile)
    
    expected = pd.DataFrame({"Reads passing QC": [177014, 59900],
                             "Nonhost reads": [59944, 59900],
                             "Percent host": [0.66136, 0],
                             "Unmapped nonhost RPK": [4698.0, 29696.0]},
                            index=["hostile_report", "hostile_report_1"])
    
    assert (actual.values == expected.values).all()


def test_fastq_gz_read_counter():
    n_reads = rr.count_reads_fastq_gz("tests/test_data/small_R1.fq.gz")
    assert n_reads == 3