#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: CPC2.py
Author: Yu-jian Kang
Created on: 2022-01-21
Last modified by CrazyHsu at 2022-01-21
'''

# !/usr/bin/env python

# This version is edited at 20190124 for compatibility on OS system

import sys, os, re, time, gzip
import numpy as np
from Bio.Seq import Seq
from Bio import SeqIO

from iflas.pu.isoelectric import isoelectric_point as _bio176_isoelectric_point
from optparse import OptionParser, OptionGroup


class FindCDS:
    '''
    Find the most like CDS in a given sequence
    The most like CDS is the longest ORF found in the sequence
    When having same length, the upstream ORF is printed
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar
    '''

    def __init__(self, seq):
        self.seq = seq
        self.result = (0, 0, 0, 0, 0)
        self.longest = 0
        self.basepair = {"A": "T", "T": "A", "U": "A", "C": "G", "G": "C", "N": "N", "X": "X"}

    def _reversecompliment(self):
        return "".join(self.basepair[base] for base in self.seq)[::-1]

    def get_codons(self, frame_number):
        '''
        Record every nucleotide triplet and its coordinate position for input sequence in one frame
        '''
        coordinate = frame_number
        while coordinate + 3 <= len(self.seq):
            yield (self.seq[coordinate:coordinate + 3], coordinate)
            coordinate += 3

    def find_longest_in_one(self, myframe, direction, start_codon, stop_codon):
        '''
        find the longest ORF in one reading frame: from a start codon to the next stop codon (integrity 1) or to the
        last codon (integrity -1); among ORFs of the same length the most upstream one
        '''
        triplet_got = self.get_codons(myframe)
        for codon, index in triplet_got:
            if codon not in start_codon or codon in stop_codon:
                continue
            orf_start, integrity = index, -1
            for codon, index in triplet_got:
                if codon in stop_codon:
                    integrity = 1
                    break
            orf_end = index + 3
            Length = (orf_end - orf_start)
            if Length > self.longest:
                self.longest = Length
                self.result = [direction, orf_start, orf_end, Length, integrity]
            if Length == self.longest and orf_start < self.result[1]:
                self.result = [direction, orf_start, orf_end, Length, integrity]

    def longest_orf(self, direction, start_codon={"ATG": None}, stop_codon={"TAG": None, "TAA": None, "TGA": None}):
        return_orf = ""
        for frame in range(3):
            self.find_longest_in_one(frame, "+", start_codon, stop_codon)
        return_orf = self.seq[self.result[1]:self.result[2]][:]
        start_coordinate = self.result[1]
        strand_direction = "+"
        orf_integrity = self.result[4]
        '''
        Also check reverse chain if -r is chosen
        '''
        if direction == "-":
            self.seq = self._reversecompliment()
            for frame in range(3):
                self.find_longest_in_one(frame, "-", start_codon, stop_codon)
            if self.result[0] == "-":
                return_orf = self.seq[self.result[1]:self.result[2]][:]
                start_coordinate = self.result[1]
                strand_direction = "-"
                orf_integrity = self.result[4]
        return return_orf, start_coordinate, strand_direction, orf_integrity


class Fickett:
    '''
    calculate Fickett TESTCODE for full sequence
    NAR 10(17) 5303-531
    modified from source code of CPAT 1.2.1 downloaded from https://sourceforge.net/projects/rna-cpat/files/?source=navbar
    '''

    def __init__(self):
        '''new compiled Fickett look-up table'''
        self.position_parameter = [1.9, 1.8, 1.7, 1.6, 1.5, 1.4, 1.3, 1.2, 1.1, 0.0]
        self.content_parameter = [0.33, 0.31, 0.29, 0.27, 0.25, 0.23, 0.21, 0.19, 0.17, 0]
        self.position_probability = {
            "A": [0.51, 0.55, 0.57, 0.52, 0.48, 0.58, 0.57, 0.54, 0.50, 0.36],
            "C": [0.29, 0.44, 0.55, 0.49, 0.52, 0.60, 0.60, 0.56, 0.51, 0.38],
            "G": [0.62, 0.67, 0.74, 0.65, 0.61, 0.62, 0.52, 0.41, 0.31, 0.17],
            "T": [0.51, 0.60, 0.69, 0.64, 0.62, 0.67, 0.58, 0.48, 0.39, 0.24],
        }
        self.position_weight = {"A": 0.062, "C": 0.093, "G": 0.205, "T": 0.154}
        self.content_probability = {
            "A": [0.40, 0.55, 0.58, 0.58, 0.52, 0.48, 0.45, 0.45, 0.38, 0.19],
            "C": [0.50, 0.63, 0.59, 0.50, 0.46, 0.45, 0.47, 0.56, 0.59, 0.33],
            "G": [0.21, 0.40, 0.47, 0.50, 0.52, 0.56, 0.57, 0.52, 0.44, 0.23],
            "T": [0.30, 0.49, 0.56, 0.53, 0.48, 0.48, 0.52, 0.57, 0.60, 0.51]
        }
        self.content_weight = {"A": 0.084, "C": 0.076, "G": 0.081, "T": 0.055}

    def look_up_position_probability(self, value, base):
        '''
        look up positional probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx, val in enumerate(self.position_parameter):
            if (float(value) >= val):
                return float(self.position_probability[base][idx]) * float(self.position_weight[base])

    def look_up_content_probability(self, value, base):
        '''
        look up content probability by base and value
        '''
        if float(value) < 0:
            return None
        for idx, val in enumerate(self.content_parameter):
            if (float(value) >= val):
                return float(self.content_probability[base][idx]) * float(self.content_weight[base])

    def fickett_value(self, dna):
        '''
        calculate Fickett value from full RNA transcript sequence
        '''
        if len(dna) < 2:
            return 0
        fickett_score = 0
        dna = dna
        total_base = len(dna)
        A_content = float(dna.count("A")) / total_base
        C_content = float(dna.count("C")) / total_base
        G_content = float(dna.count("G")) / total_base
        T_content = float(dna.count("T")) / total_base

        phase_0 = dna[::3]
        phase_1 = dna[1::3]
        phase_2 = dna[2::3]

        phase_0_A = phase_0.count("A")
        phase_1_A = phase_1.count("A")
        phase_2_A = phase_2.count("A")
        phase_0_C = phase_0.count("C")
        phase_1_C = phase_1.count("C")
        phase_2_C = phase_2.count("C")
        phase_0_G = phase_0.count("G")
        phase_1_G = phase_1.count("G")
        phase_2_G = phase_2.count("G")
        phase_0_T = phase_0.count("T")
        phase_1_T = phase_1.count("T")
        phase_2_T = phase_2.count("T")

        A_content = float(phase_0_A + phase_1_A + phase_2_A) / total_base
        C_content = float(phase_0_C + phase_1_C + phase_2_C) / total_base
        G_content = float(phase_0_G + phase_1_G + phase_2_G) / total_base
        T_content = float(phase_0_T + phase_1_T + phase_2_T) / total_base
        A_position = np.max([phase_0_A, phase_1_A, phase_2_A]) / (np.min([phase_0_A, phase_1_A, phase_2_A]) + 1.0)
        C_position = np.max([phase_0_C, phase_1_C, phase_2_C]) / (np.min([phase_0_C, phase_1_C, phase_2_C]) + 1.0)
        G_position = np.max([phase_0_G, phase_1_G, phase_2_G]) / (np.min([phase_0_G, phase_1_G, phase_2_G]) + 1.0)
        T_position = np.max([phase_0_T, phase_1_T, phase_2_T]) / (np.min([phase_0_T, phase_1_T, phase_2_T]) + 1.0)

        fickett_score += self.look_up_content_probability(A_content, "A")
        fickett_score += self.look_up_content_probability(C_content, "C")
        fickett_score += self.look_up_content_probability(G_content, "G")
        fickett_score += self.look_up_content_probability(T_content, "T")

        fickett_score += self.look_up_position_probability(A_position, "A")
        fickett_score += self.look_up_position_probability(C_position, "C")
        fickett_score += self.look_up_position_probability(G_position, "G")
        fickett_score += self.look_up_position_probability(T_position, "T")

        return fickett_score


# ===================
def gz_file(fq_file,mode,level=6):
    try:
        if fq_file.endswith("gz"):
            fq_fp = gzip.open(fq_file,mode+"b",level)
        else:
            fq_fp = open(fq_file,mode)  # Py2 file() builtin removed in Py3 (class-a crash fix)
    except:
        sys.stderr.write("Error: Fail to IO file: %s\n"%(fq_file))
        sys.exit(1)
    return fq_fp



def fasta_read(fn):
    """
    ID: gi|2765658|emb|Z78533.1|CIZ78533
    Name: gi|2765658|emb|Z78533.1|CIZ78533
    Description: gi|2765658|emb|Z78533.1|CIZ78533 C.irapeanum 5.8S rRNA gene and ITS1 and ITS2 DNA
    Number of features: 0
    Seq('CGTAACAAGGTTTCCGTAGGTGAACCTGCGGAAGGATCATTGATGAGACCGTGG...GGG', SingleLetterAlphabet())
    """
    with gz_file(fn, "r") as f:
        for seq in SeqIO.parse(f, "fasta"):
            yield seq

########################
def mRNA_translate(mRNA):
    return Seq(mRNA).translate()


def protein_param(peptide):
    # BUG-145: pI must match the published Python-2 run (Biopython 1.76), not the installed Biopython.
    return _bio176_isoelectric_point(peptide)


def _find_utils_dir():
    """Locate the iFLAS ``utils/`` dir (LIBSVM binaries + cpc2 model) from the packaged
    location. Honors ``IFLAS_UTILS_DIR``; otherwise walks up from this file looking for a
    ``utils/`` holding ``cpc2.model`` (repo-root in editable installs)."""
    env = os.environ.get("IFLAS_UTILS_DIR")
    if env and os.path.isdir(env):
        return env
    d = os.path.dirname(os.path.abspath(__file__))
    for _ in range(6):
        cand = os.path.join(d, "utils")
        if os.path.exists(os.path.join(cand, "cpc2.model")):
            return cand
        d = os.path.dirname(d)
    return os.path.join(os.path.dirname(os.path.abspath(__file__)), "utils")


def calculate_potential(fasta, strand, output_orf, outfile, runner=None):
    '''
    Calculate three features: putative peptide length,pI and Fickett
    And assess coding potential based on SVM model
    '''
    strinfoAmbiguous = re.compile("X|B|Z|J|U", re.I)
    ptU = re.compile("U", re.I)
    with open(outfile + ".feat", "w") as ftmp_feat, open(outfile + ".tmp.1", "w") as ftmp_svm, \
            open(outfile, "w") as ftmp_result:
        if output_orf == 1:
            my_header = ["#ID", "transcript_length", "peptide_length", "Fickett_score", "pI", "ORF_integrity", "ORF_Start",
                         "coding_probability", "label"]
        else:
            my_header = ["#ID", "transcript_length", "peptide_length", "Fickett_score", "pI", "ORF_integrity",
                         "coding_probability", "label"]
        ftmp_result.write("\t".join(map(str, my_header)) + "\n")
        fickett_obj = Fickett()
        for seq in fasta_read(fasta):
            seqid = seq.id
            seqRNA = ptU.sub("T", str(seq.seq).strip())
            '''seqRNA:transcript full sequence'''
            seqRNA = seqRNA.upper()
            seqCDS, start_pos, orf_strand, orf_fullness = FindCDS(seqRNA).longest_orf(strand)
            '''seqCDS:longest ORF'''
            seqprot = mRNA_translate(seqCDS)
            pep_len = len(seqprot)  # pep_len = len(seqprot.strip("*"))
            newseqprot = strinfoAmbiguous.sub("", str(seqprot))
            '''exclude ambiguous amio acid X, B, Z, J, Y in peptide sequence'''
            fickett_score = fickett_obj.fickett_value(seqRNA)
            if pep_len > 0:
                start_pos = start_pos + 1
                isoelectric_point = protein_param(newseqprot.strip("*"))
            else:
                start_pos = 0
                orf_fullness = -1
                isoelectric_point = 0.0
            if output_orf == 1:
                output_line = [seqid, len(seqRNA), pep_len, fickett_score, isoelectric_point, orf_fullness, start_pos]
            else:
                output_line = [seqid, len(seqRNA), pep_len, fickett_score, isoelectric_point, orf_fullness]
            ftmp_result.write("\t".join(map(str, output_line)) + "\n")
            ftmp_feat.write(
                "\t".join(map(str, [seqid, len(seqRNA), pep_len, fickett_score, isoelectric_point, orf_fullness])) + "\n")
            ftmp_svm.write("".join(map(str, ["999", " 1:", pep_len, " 2:", fickett_score, " 3:", isoelectric_point, " 4:",
                                             orf_fullness])) + "\n")

    # coding probability with LIBSVM
    data_dir = _find_utils_dir() + os.sep   # resolve repo-root utils/ from packaged location
    lib_dir = data_dir
    app_svm_scale = os.path.join(lib_dir, "svm-scale")
    app_svm_predict = os.path.join(lib_dir, "svm-predict")
    for app in (app_svm_scale, app_svm_predict):
        if not os.access(app, os.X_OK):
            from ..runtime import ToolError
            raise ToolError(app, 127, "[ERROR] No executable {} on CPC2 path!".format(os.path.basename(app)))

    if runner is None:
        from ..runtime import ToolRunner
        runner = ToolRunner()

    with open(outfile + ".tmp.2", "w") as stdout:
        runner.run([app_svm_scale, "-r", data_dir + "cpc2.range", outfile + ".tmp.1"], stdout=stdout)
    runner.run([app_svm_predict, "-b", "1", "-q", outfile + ".tmp.2", data_dir + "cpc2.model",
                outfile + ".tmp.out"])

    prob = {}           # libsvm output, by input line
    with open(outfile + '.tmp.out', 'r') as tmp_file:
        i = 0
        for line in tmp_file:
            i = i + 1
            array = line.split(' ')
            if array[0] == "1":
                label = "coding"
            else:
                label = "noncoding"
            prob[i] = str(array[1]) + '\t' + label
    # paste to features
    with open(outfile, 'r') as tmp_file, open(outfile + '.txt', 'w') as out_file:
        i = 0
        for line in tmp_file:
            i = i + 1
            if i == 1:
                out_file.write(line)
            else:
                line = line.rstrip('\n')
                out_file.write(line)
                out_file.write('\t' + prob[i] + '\n')
    for suffix in ('.tmp.1', '.tmp.2', '.tmp.out', '', '.feat'):
        path = outfile + suffix
        if os.path.exists(path):
            os.remove(path)
    return 0


def __main():
    start_time = time.time()
    usage = "usage: %prog [options] -i input.fasta -o output_file"
    description = "Contact: Kang Yujian <kangyj@mail.cbi.pku.edu.cn>"
    parser = OptionParser(usage, version="%prog 0.1", description=description)
    Common_group = OptionGroup(parser, "Common Options")
    Common_group.add_option("-i", dest="fasta", help="input sequence in fasta format [Required]", metavar="FILE",
                            type="string", default=None)
    Common_group.add_option("-o", dest="outfile", help="output file [Default: cpc2output.txt]", metavar="FILE",
                            type="string", default="cpc2output.txt")
    Common_group.add_option("-r", dest="reverse", help="also check the reverse strand [Default: FALSE]",
                            action="store_true")
    Common_group.add_option("--ORF", dest="ORF", help="output the start position of longest ORF [Default: FALSE]",
                            action="store_true")
    parser.add_option_group(Common_group)
    (options, args) = parser.parse_args()
    if options.fasta == None:
        parser.print_help()
        return -1
    else:
        if not os.path.isfile(options.fasta):
            sys.stderr.write("[ERROR] %s is not a file\n" % options.fasta)
            return -1
    if options.reverse:
        strand = "-"
    else:
        strand = "+"
    if options.ORF:
        output_orf = 1
    else:
        output_orf = 0
    if calculate_potential(options.fasta, strand, output_orf, options.outfile):
        return 1
    sys.stderr.write("[INFO] cost time: %ds\n" % (time.time() - start_time))
    return 0


if __name__ == "__main__":
    sys.exit(__main())
