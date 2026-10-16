"""Unit tests for iflas.io.files (filesystem + simple file-parsing helpers)."""

import os
from types import SimpleNamespace

from iflas.io.files import (
    getDictFromFile, makeLink, removeFiles,
    resolveDir, validateDir, validateFaAndFqFile, validateFile,
)


def test_get_dict_from_file(tmp_path):
    f = tmp_path / "kv.tsv"
    f.write_text("a\t1\nb\t2\n")
    assert getDictFromFile(str(f)) == {"a": "1", "b": "2"}
    # inline-separated value -> list
    f2 = tmp_path / "kv2.tsv"
    f2.write_text("a\tx,y\n")
    assert getDictFromFile(str(f2), inlineSep=",") == {"a": ["x", "y"]}
    # valueCol=None -> whole row as value
    assert getDictFromFile(str(f), valueCol=None) == {"a": ["a", "1"], "b": ["b", "2"]}



def test_validate_file_and_dir(tmp_path):
    f = tmp_path / "exists.txt"
    f.write_text("hi")
    assert validateFile(str(f)) is True
    assert validateFile(str(tmp_path / "nope.txt")) is False
    assert validateFile([f]) is False          # list -> False
    assert validateFile(str(tmp_path)) is False  # a dir is not a file
    assert validateDir(str(tmp_path)) is True
    assert validateDir(str(f)) is False
    assert validateDir(str(tmp_path / "nope")) is False


def test_dir_and_link_ops(tmp_path):
    sub = tmp_path / "a" / "b"
    resolveDir(str(sub), chdir=False)
    assert sub.is_dir()

    src = tmp_path / "src.txt"
    src.write_text("data")
    link = tmp_path / "link.txt"
    makeLink(str(src), str(link))
    assert os.path.islink(str(link))

    removeFiles(fileList=[str(src)])
    assert not src.exists()


class _SeqkitRunner:
    def __init__(self, stat_stdout: str, head_returncode: int = 0, stat_returncode: int = 0):
        self.stat_stdout = stat_stdout
        self.head_returncode = head_returncode
        self.stat_returncode = stat_returncode
        self.calls = []

    def run(self, cmd, **kwargs):
        self.calls.append((cmd, kwargs))
        if cmd == ["seqkit", "head", "-n", "100", "reads.fa"]:
            return SimpleNamespace(stdout=">r1\nACGT\n", stderr="", returncode=self.head_returncode)
        if cmd == ["seqkit", "stat"]:
            assert kwargs["input"] == ">r1\nACGT\n"
            return SimpleNamespace(stdout=self.stat_stdout, stderr="", returncode=self.stat_returncode)
        raise AssertionError(f"unexpected command: {cmd!r}")


def test_validate_fa_and_fq_file_routes_seqkit_through_runner():
    runner = _SeqkitRunner("file\tformat\nreads.fa\tFASTA\n")

    assert validateFaAndFqFile("reads.fa", runner=runner) == "fasta"
    assert runner.calls == [
        (["seqkit", "head", "-n", "100", "reads.fa"], {"check": False}),
        (["seqkit", "stat"], {"input": ">r1\nACGT\n", "check": False}),
    ]


def test_validate_fa_and_fq_file_preserves_fastq_and_false_results():
    assert validateFaAndFqFile("reads.fa", runner=_SeqkitRunner("file\tformat\nreads.fa\tFASTQ\n")) == "fastq"
    assert validateFaAndFqFile("reads.fa", runner=_SeqkitRunner("file\tformat\nreads.fa\tUNKNOWN\n")) is False
    assert validateFaAndFqFile("reads.fa", runner=_SeqkitRunner("", head_returncode=1)) is False
