import argparse
import subprocess
import pytest

import homi_pipeline.HoMi as homi


def make_args(unlock=False):
    return argparse.Namespace(config="config.yaml", profile=None, cores=1,
                              conda_prebuilt=False, unlock=unlock,
                              workdir=None, snakemake_extra=None)


def make_completed_process(returncode):
    return subprocess.CompletedProcess(args=["snakemake"], returncode=returncode)


def patch_main_dependencies(monkeypatch, args, returncodes):
    completed_processes = iter([make_completed_process(code) for code in returncodes])
    monkeypatch.setattr(homi, "get_args", lambda: args)
    monkeypatch.setattr(homi, "read_config", lambda path: {})
    monkeypatch.setattr(homi, "run_checker", lambda config: None)
    monkeypatch.setattr(homi, "get_snakefile_path", lambda: "snakefile")
    monkeypatch.setattr(homi.subprocess, "run",
                        lambda command: next(completed_processes))


@pytest.mark.parametrize("returncode", [0, 1])
def test_main_exits_with_snakemake_returncode(monkeypatch, returncode):
    patch_main_dependencies(monkeypatch, make_args(), [returncode])
    with pytest.raises(SystemExit) as exit_info:
        homi.main()
    assert exit_info.value.code == returncode


def test_main_exits_if_unlock_fails(monkeypatch):
    # the second returncode is for the main run, which must not start
    patch_main_dependencies(monkeypatch, make_args(unlock=True), [1, 0])
    with pytest.raises(SystemExit) as exit_info:
        homi.main()
    assert exit_info.value.code == 1
