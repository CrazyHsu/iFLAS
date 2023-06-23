#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: preprocess.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

from commonFuncs import *

def processRnaseq(dataObj=None, threads=None, dirSpec=None, max_reads_length_tirmmed=30):
    print getCurrentTime() + " Filtering short-reads with fastp for project {} sample {}...".format(dataObj.project_name, dataObj.sample_name)

    logDir = os.path.join(dirSpec.out_dir, dataObj.project_name, dataObj.sample_name, "log")
    resolveDir(logDir, chdir=False)
    if dataObj.ngs_reads_paired == "paired" and dataObj.ngs_left_reads and dataObj.ngs_right_reads:
        if ("fastp" in dataObj.ngs_left_reads and "fastp" in dataObj.ngs_right_reads) or \
                (checkFastpReadsExist(dataObj.ngs_left_reads) and checkFastpReadsExist(dataObj.ngs_right_reads)):
            print getCurrentTime() + " It seems like you have filter the reads with fastp, so we skip the process!"
            print getCurrentTime() + " Filtering short-reads with fastp done!"
            return
        leftReadsRepeats = [i.strip() for i in dataObj.ngs_left_reads.split(";")]
        rightReadsRepeats = [i.strip() for i in dataObj.ngs_right_reads.split(";")]
        if len(leftReadsRepeats) != len(rightReadsRepeats):
            raise Exception("The repeats of your NGS data not match between your left reads and right reads")
        else:
            newLeftReadsRepeats = []
            newRightReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                rightReads = rightReadsRepeats[i].split(",")
                if len(leftReads) != len(rightReads):
                    raise Exception("please input the right paired NGS reads for the analysis")
                else:
                    newLeftReads = []
                    newRightReads = []
                    for j in range(len(leftReads)):
                        leftReadsDir = os.path.dirname(leftReads[j])
                        rightReadsDir = os.path.dirname(rightReads[j])
                        leftReadsBase = os.path.basename(leftReads[j])
                        rightReadsBase = os.path.basename(rightReads[j])
                        newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                        newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                        cmd = "fastp -i {} -I {} -o {} -O {} -w {} -q 20 -l {} -j {}/fastp.json -h {}/fastp.html 2>/dev/null"
                        cmd = cmd.format(leftReads[j], rightReads[j], newLeft, newRight, threads,
                                         int(dataObj.ngs_reads_length) - max_reads_length_tirmmed, logDir, logDir)
                        subprocess.call(cmd, shell=True)
                        newLeftReads.append(newLeft)
                        newRightReads.append(newRight)
                    newLeftReadsRepeats.append(",".join(newLeftReads))
                    newRightReadsRepeats.append(",".join(newRightReads))
            dataObj.ngs_left_reads = ";".join(newLeftReadsRepeats)
            dataObj.ngs_right_reads = ";".join(newRightReadsRepeats)
    else:
        if dataObj.ngs_left_reads and dataObj.ngs_right_reads == None:
            if "fastp" in dataObj.ngs_left_reads or checkFastpReadsExist(dataObj.ngs_left_reads):
                print getCurrentTime() + " It seems like you have filter the reads with fastp, so we skip the process!"
                print getCurrentTime() + " Filtering short-reads with fastp done!"
                return
            leftReadsRepeats = [i.strip() for i in dataObj.ngs_left_reads.split(";")]
            newLeftReadsRepeats = []
            for i in range(len(leftReadsRepeats)):
                leftReads = leftReadsRepeats[i].split(",")
                newLeftReads = []
                for j in range(len(leftReads)):
                    leftReadsDir = os.path.dirname(leftReads[j])
                    leftReadsBase = os.path.basename(leftReads[j])
                    newLeft = os.path.join(leftReadsDir, "fastp.{}".format(leftReadsBase))
                    cmd = "fastp -i {} -o {} -w {} -q 20 -l {} -j {}/fastp.json -h {}/fastp.html 2>/dev/null"
                    cmd = cmd.format(leftReads[j], newLeft, threads,
                                     int(dataObj.ngs_reads_length) - max_reads_length_tirmmed, logDir, logDir)
                    subprocess.call(cmd, shell=True)
                    newLeftReads.append(newLeft)
                newLeftReadsRepeats.append(",".join(newLeftReads))
            dataObj.ngs_left_reads = ";".join(newLeftReadsRepeats)
        elif dataObj.ngs_right_reads and dataObj.ngs_left_reads == None:
            if "fastp" in dataObj.ngs_right_reads or checkFastpReadsExist(dataObj.ngs_right_reads):
                print getCurrentTime() + " It seems like you have filter the reads with fastp, so we skip the process!"
                print getCurrentTime() + " Filtering short-reads with fastp done!"
                return
            rightReadsRepeats = [i.strip() for i in dataObj.ngs_right_reads.split(";")]
            newRightReadsRepeats = []
            for i in range(len(rightReadsRepeats)):
                rightReads = rightReadsRepeats[i].split(",")
                newRightReads = []
                for j in range(len(rightReads)):
                    rightReadsDir = os.path.dirname(rightReads[j])
                    rightReadsBase = os.path.basename(rightReads[j])
                    newRight = os.path.join(rightReadsDir, "fastp.{}".format(rightReadsBase))
                    cmd = "fastp -i {} -o {} -w {} -q 20 -l {} -j {}/fastp.json -h {}/fastp.html 2>/dev/null"
                    cmd = cmd.format(rightReads[j], newRight, threads,
                                     int(dataObj.ngs_reads_length) - max_reads_length_tirmmed, logDir, logDir)
                    subprocess.call(cmd, shell=True)
                    newRightReads.append(newRight)
                newRightReadsRepeats.append(",".join(newRightReads))
            dataObj.ngs_right_reads = ";".join(newRightReadsRepeats)
        else:
            raise Exception("The NGS data seem not to be paired or single, please check it")
    print getCurrentTime() + " Filtering short-reads with fastp for project {} sample {} done!".format(dataObj.project_name, dataObj.sample_name)


def correctWithFmlrc2(dataObj, dirSpec=None, useFmlrc2=True, threads=None):
    if useFmlrc2:
        projectName, sampleName = dataObj.project_name, dataObj.sample_name
        print getCurrentTime() + " Correct long reads in project {} sample {} with short-reads...".format(projectName, sampleName)
        prevDir = os.getcwd()
        baseDir = os.path.join(dirSpec.out_dir, projectName, sampleName)
        resolveDir(os.path.join(baseDir, "preprocess", "fmlrc"))
        logDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "log")
        resolveDir(logDir, chdir=False)
        if dataObj.ngs_left_reads != None and dataObj.ngs_right_reads != None:
            leftReadsList = [i.strip() for i in re.split("[;|,]", dataObj.ngs_left_reads)]
            rightReadsList = [i.strip() for i in re.split("[;|,]", dataObj.ngs_right_reads)]
            cmd = "zcat {} > leftReads.fastq".format(" ".join(leftReadsList))
            subprocess.call(cmd, shell=True)
            cmd = "zcat {} > rightReads.fastq".format(" ".join(rightReadsList))
            subprocess.call(cmd, shell=True)
            cmd = "cat leftReads.fastq rightReads.fastq | awk 'NR % 4 == 2' | tr NT TN | "
        elif dataObj.ngs_left_reads != None or dataObj.ngs_right_reads == None:
            leftReadsList = [i.strip() for i in re.split("[;|,]", dataObj.ngs_left_reads)]
            cmd = "zcat {} > leftReads.fastq".format(" ".join(leftReadsList))
            subprocess.call(cmd, shell=True)
            cmd = "cat leftReads.fastq | awk 'NR % 4 == 2' | tr NT TN | "
        else:
            rightReadsList = [i.strip() for i in re.split("[;|,]", dataObj.ngs_right_reads)]
            cmd = "zcat {} > rightReads.fastq".format(" ".join(rightReadsList))
            subprocess.call(cmd, shell=True)
            cmd = "cat rightReads.fastq | awk 'NR % 4 == 2' | tr NT TN | "
        cmd = cmd + "ropebwt2 -LR 2>{}/{}.ropebwt2.log | tr NT TN | fmlrc2-convert comp_msbwt.npy 1>{}/{}.fmlrc_convert.log 2>&1".format(
            logDir, sampleName, logDir, sampleName)
        subprocess.call(cmd, shell=True)
        if dataObj.data_processed_location == None and dataObj.data_location:
            dataObj.data_processed_location = os.path.join(baseDir, "preprocess", dataObj.tgs_plat.lower(), "rawFlnc.fq")
        else:
            makeLink(dataObj.data_processed_location, "rawFlnc.fq")
        cmd = "seqkit seq --rna2dna {} -w 0 > raw.dna.fastq".format(dataObj.data_processed_location)
        subprocess.call(cmd, shell=True)
        cmd = "fmlrc2 -t {} -k 25 59 -C 10 comp_msbwt.npy raw.dna.fastq fmlrc_corrected.fasta 1>{}/{}.fmlrc.log 2>&1".format(
            threads, logDir, sampleName)
        subprocess.call(cmd, shell=True)
        dataObj.data_processed_location = os.path.join(os.getcwd(), "fmlrc_corrected.fasta")
        # removeFiles(os.getcwd(), ["leftReads.fastq", "rightReads.fastq"])
        os.chdir(prevDir)
        print getCurrentTime() + " Correct long reads in project {} sample {} with short-reads done!".format(projectName, sampleName)


def retrievePacbio(dataObj=None, ccsParams=None, dirSpec=None, threads=10):
    projectName, strategy, sampleName = dataObj.project_name, dataObj.strategy, dataObj.sample_name
    print getCurrentTime() + " Extract PacBio reads for project {} sample {}...".format(projectName, sampleName)
    workDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "preprocess", "pacbio")
    prevDir = os.getcwd()
    resolveDir(workDir)

    if strategy.lower() == "rsii":
        print getCurrentTime() + " The reads in project {} entry {} are from RSII".format(projectName, sampleName)
        if not os.path.isdir(dataObj.dataLocation):
            raise ValueError("The RSii data should be stored in a directory, please check it!")

        tgsDir = dataObj.dataLocation
        cmd = "find {} -name '*.bax.h5' > {}.fofn".format(tgsDir, sampleName)
        subprocess.call(cmd, shell=True)
        cmd = "bax2bam -f {}.fofn -o {}".format(dataObj.uniqName, sampleName)
        subprocess.call(cmd, shell=True)

        inputBam = "{}.subreads.bam".format(dataObj.uniqName)
        print getCurrentTime() + " Extract CCS reads from bax.h5 files and demultiplexing for project {} sample {}...".format(projectName, sampleName)
        ccsOut = subprocess.Popen(["ccs", "--version"], stdout=subprocess.PIPE).communicate()[0].strip("\n")
        if "ccs 3." not in ccsOut:
            raise ValueError("The RSii data only support pbccs version 3.*, please use the correct version")

        cmd = "ccs {} CCS.bam -j {} --noPolish --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
              "{}".format(inputBam, threads, ccsParams.min_pass, ccsParams.min_subread_length,
                          ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
        subprocess.call(cmd, shell=True)
        print getCurrentTime() + " Extract CCS reads from bax.h5 files for project {} entry {} done!".format(projectName, sampleName)
    elif strategy.lower() == "sequel":
        print getCurrentTime() + " The reads in project {} entry {} are from Sequel".format(projectName, sampleName)
        if os.path.isdir(dataObj.data_location):
            tgsDir = dataObj.data_location
            xmlFiles = os.popen("find -L {} -name '*.subreadset.xml'".format(tgsDir))
            xmlFileStr = " ".join([i.strip("\n") for i in xmlFiles])
            if not xmlFileStr:
                bamFiles = os.popen("find -L {} -name '*.bam'".format(tgsDir))
                bamFileStr = " ".join([i.strip("\n") for i in bamFiles])
                if not bamFileStr:
                    raise Exception("You should put the *.bam in the {} directory".format(tgsDir))
                fileStr = bamFileStr
            else:
                fileStr = xmlFileStr
        else:
            tgsFiles = [i.strip() for i in dataObj.data_location.split(",")]
            for i in tgsFiles:
                if not os.path.exists(i + ".pbi"):
                    cmd = "pbindex {}".format(i)
                    subprocess.call(cmd, shell=True)
            tgsFilesStr = " ".join(tgsFiles)
            fileStr = tgsFilesStr

        subreads = "{}.subreadset.xml".format(sampleName)
        cmd = "dataset create --type SubreadSet --force --name {} {} {}".format(projectName, subreads, fileStr)
        subprocess.call(cmd, shell=True)

        print getCurrentTime() + " Extract CCS reads from subreads files for project {} sample {}...".format(projectName, sampleName)
        ccsOut = subprocess.Popen(["ccs", "--version"], stdout=subprocess.PIPE).communicate()[0].strip("\n")
        if "ccs 3." in ccsOut:
            cmd = "ccs {} CCS.bam -j {} --minPasses {} --minLength {} --minPredictedAccuracy {} --minReadScore " \
                  "{} --force".format(subreads, threads, ccsParams.min_pass, ccsParams.min_subread_length,
                                      ccsParams.min_predicted_accuracy, ccsParams.min_read_score)
            subprocess.call(cmd, shell=True)
        elif "ccs 4." in ccsOut or "ccs 5." in ccsOut or "ccs 6." in ccsOut:
            cmd = "ccs {} CCS.bam -j {} --min-passes {} --min-length {} " \
                  "--min-rq {}".format(subreads, threads, ccsParams.min_pass, ccsParams.min_subread_length,
                                       ccsParams.min_predicted_accuracy)
            subprocess.call(cmd, shell=True)
        else:
            raise Exception("You should use version of ccs up to 3.4 or higher")
        print getCurrentTime() + " Extract CCS reads from subreads files for project {} sample {} done!".format(projectName, sampleName)
    else:
        raise Exception("You should input the right PacBio sequencing strategy!!")

    print getCurrentTime() + " Demultiplex CCS bam for project {} sample {}...".format(projectName, sampleName)
    cmd = "lima CCS.bam {} fl.bam --isoseq --dump-clips -j {}".format(dataObj.primer, threads)
    subprocess.call(cmd, shell=True)
    print getCurrentTime() + " Demultiplex CCS bam for project {} sample {} done!".format(projectName, sampleName)

    print getCurrentTime() + " Refine CCS bam for project {} sample {}...".format(projectName, sampleName)
    tmpBamName = "fl.5p--{}_3p.bam".format(sampleName, sampleName)
    cmd = "isoseq3 refine {} {} flnc.bam --require-polya -j {}".format(tmpBamName, dataObj.primer, threads)
    subprocess.call(cmd, shell=True)
    cmd = "samtools fastq flnc.bam > rawFlnc.fq"
    subprocess.call(cmd, shell=True)
    dataObj.data_processed_location = os.path.join(os.getcwd(), "rawFlnc.fq")
    os.chdir(prevDir)
    print getCurrentTime() + " Refine CCS bam for project {} sample {} done!".format(projectName, sampleName)


def retrieveNanopore(dataObj=None, dirSpec=None, threads=10, flowcellType="FLO-MIN106", kitType="SQK-RNA002"):
    projectName, sampleName = dataObj.project_name, dataObj.sample_name
    print getCurrentTime() + " Calling fastq files from fast5 files for project {} sample {}...".format(projectName, sampleName)
    workDir = os.path.join(dirSpec.out_dir, projectName, sampleName, "preprocess", "nanopore")
    prevDir = os.getcwd()
    resolveDir(workDir)
    fast5Dir = dataObj.data_location
    nanoporeFileType = checkFast5Files(fast5Dir)
    if nanoporeFileType == "fast5":
        basecallingOut = "guppy_basecalling"
        cmd = "guppy_basecaller.py --flowcell {} --kit {} -i {} --save_path {} --cpu_threads_per_caller {}" \
              " --number_callers {}".format(flowcellType, kitType, fast5Dir, basecallingOut, threads, 1)
        subprocess.call(cmd, shell=True)
        rawFlncFq = "rawFlnc.fq"
        cmd = "cat {}/*.fastq > {}".format(basecallingOut, rawFlncFq)
        subprocess.call(cmd, shell=True)
        dataObj.data_processed_location = os.path.join(os.getcwd(), rawFlncFq)
    else:
        rawFlncFq = "rawFlnc.fq"
        cmd = "cat {}/* > {}".format(dataObj.data_location, rawFlncFq)
        subprocess.call(cmd, shell=True)
        dataObj.data_processed_location = os.path.join(os.getcwd(), rawFlncFq)
    os.chdir(prevDir)
    print getCurrentTime() + " Calling fastq files from fast5 files for project {} sample {} done!".format(projectName, sampleName)


def preprocess(dataObj=None, ccsParams=None, dirSpec=None, threads=10):
    if dataObj.data_processed_location != None and validateFile(dataObj.data_processed_location):
        if validateFaAndFqFile(dataObj.data_processed_location) in ["fastq", "fasta"]:
            print getCurrentTime() + " It seems like you provided an valid processed fasta/fastq file, we will not call reads from raw data!"
            if dataObj.ngs_left_reads != None or dataObj.ngs_right_reads != None:
                processRnaseq(dataObj=dataObj, threads=threads, dirSpec=dirSpec, max_reads_length_tirmmed=1)
                renameNGSdata2fastp(dataObj=dataObj)
                correctWithFmlrc2(dataObj, useFmlrc2=dataObj.use_fmlrc2, dirSpec=dirSpec, threads=threads)
            else:
                raise Exception(getCurrentTime() + " You should set the ngs reads to correct the flnc reads!")
        else:
            raise Exception(getCurrentTime() + " It seems like you provided an unvalid processed fasta/fastq file, please check it!")
    else:
        if validateFile(dataObj.data_location) or validateDir(dataObj.data_location):
            if dataObj.tgs_plat.lower() == "pacbio" and (dataObj.ngs_left_reads != None or dataObj.ngs_right_reads != None):
                retrievePacbio(dataObj=dataObj, ccsParams=ccsParams, dirSpec=dirSpec, threads=threads)
                processRnaseq(dataObj=dataObj, threads=threads, dirSpec=dirSpec, max_reads_length_tirmmed=1)
                renameNGSdata2fastp(dataObj=dataObj)
                correctWithFmlrc2(dataObj, useFmlrc2=dataObj.use_fmlrc2, dirSpec=dirSpec, threads=threads)
            elif dataObj.tgs_plat.lower() == "nanopore" and (dataObj.ngs_left_reads != None or dataObj.ngs_right_reads != None):
                retrieveNanopore(dataObj=dataObj, dirSpec=dirSpec, threads=threads, flowcellType=dataObj.flowcell_type, kitType=dataObj.kit_type)
                processRnaseq(dataObj=dataObj, threads=threads, dirSpec=dirSpec, max_reads_length_tirmmed=1)
                renameNGSdata2fastp(dataObj=dataObj)
                correctWithFmlrc2(dataObj, useFmlrc2=dataObj.use_fmlrc2, dirSpec=dirSpec, threads=threads)
            elif dataObj.tgs_plat == None and (dataObj.ngs_left_reads != None or dataObj.ngs_right_reads != None):
                processRnaseq(dataObj=dataObj, threads=threads, max_reads_length_tirmmed=1)
                renameNGSdata2fastp(dataObj=dataObj)
            else:
                raise Exception(getCurrentTime() + "The TGS platform you input can't be identified, please check it!")
        else:
            if dataObj.ngs_left_reads != None or dataObj.ngs_right_reads != None:
                processRnaseq(dataObj=dataObj, threads=threads, dirSpec=dirSpec, max_reads_length_tirmmed=1)
                renameNGSdata2fastp(dataObj=dataObj)
            else:
                raise Exception(getCurrentTime() + "Please check the data your input is valid!")
