import pytest
import pandas as pd
from homi_pipeline.rule_utils import read_reports as rr


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
    expected = {"hostile_report": (1/3)*(10**6), "hostile_report_1": (1/2)*(10**6)}
    assert unmapped == expected
    

def test_add_unmapped_nonhost_to_hostile():
    hostile = pd.DataFrame({"Reads passing QC": [177014, 59900],
                             "Nonhost reads": [59944, 59900],
                             "Percent host": [0.66136, 0]},
                            index=["hostile_report", "hostile_report_1"])
    # note that this is in "reverse" order just to make sure it is merging correctly 
    # independent of dict order
    unmapped_nonhost =  {"hostile_report_1": (1/2)*(10**7), "hostile_report": (1/3)*(10**7)}
    
    actual = rr.add_unmapped_nonhost_to_hostile_reports(unmapped_nonhost, hostile)
    
    expected = pd.DataFrame({"Reads passing QC": [177014, 59900],
                             "Nonhost reads": [59944, 59900],
                             "Percent host": [0.66136, 0],
                             "Unmapped nonhost TPM": [(1/3)*(10**7), (1/2)*(10**7)]},
                            index=["hostile_report", "hostile_report_1"])
    
    assert (actual.values == expected.values).all()


def test_fastq_gz_read_counter():
    n_reads = rr.count_reads_fastq_gz("tests/test_data/small.R1.fq.gz")
    assert n_reads == 3


def test_add_raw_counts_to_df():
    metadata_raw_counts = pd.DataFrame({"Sample": ["small"]})
    df = pd.DataFrame({"Reads passing QC": [177014],
                       "Nonhost reads": [59944],
                       "Percent host": [0.66136],
                       "Unmapped nonhost RPK": [4698.0]},
                       index=["small"])
    actual = rr.add_original_read_counts_to_df(metadata_raw_counts, "tests/test_data", df)
    expected = pd.DataFrame({"Reads passing QC": [177014],
                             "Nonhost reads": [59944],
                             "Percent host": [0.66136],
                             "Unmapped nonhost RPK": [4698.0],
                             "Raw reads": [6]},
                             index=["small"])
    
    assert (actual.values == expected.values).all()