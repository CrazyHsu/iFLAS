#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
File name: iflas.py
Author: CrazyHsu @ crazyhsu9527@gmail.com
Created on: 2022-01-18
Last modified: 2022-01-18
'''

import sys, argparse
import pybedtools

# from commonFuncs import *
from commonObjs import *
from Config import *

def splitCommandRun(args, dataToProcess, refInfoParams, dirSpec, ccsParams, minimap2Params, collapseParams, hqIsoParams, optionTools):
    if args.command == 'preproc':
        pool = MyPool(processes=len(dataToProcess))
        from preprocess import preprocess
        for dataObj in dataToProcess:
            dataObj.single_run_threads = int(optionTools.threads / float(len(dataToProcess)))
            # preprocess(dataObj=dataObj, ccsParams=ccsParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
            pool.apply_async(preprocess, (dataObj, ccsParams, dirSpec, dataObj.single_run_threads))
        pool.close()
        pool.join()

    strain2data = {}
    for dataObj in dataToProcess:
        if dataObj.project_name not in strain2data:
            strain2data[dataObj.project_name] = {dataObj.ref_strain: {dataObj.strain: [dataObj]}}
        elif dataObj.ref_strain not in strain2data[dataObj.project_name]:
            strain2data[dataObj.project_name][dataObj.ref_strain] = {dataObj.strain: [dataObj]}
        elif dataObj.strain not in strain2data[dataObj.project_name][dataObj.ref_strain]:
            strain2data[dataObj.project_name][dataObj.ref_strain][dataObj.strain] = [dataObj]
        else:
            strain2data[dataObj.project_name][dataObj.ref_strain][dataObj.strain].append(dataObj)

    if optionTools.merge_data_from_same_strain or args.merge:
        sampleMergedToProcess = mergeSample(strain2data)
        processNum = len(list(nestedDictValues(sampleMergedToProcess, returned="vamapping.pylue")))
        pool = MyPool(processes=processNum)
        for proj in sampleMergedToProcess:
            for ref_strain in sampleMergedToProcess[proj]:
                for strain in sampleMergedToProcess[proj][ref_strain]:
                    dataObj = sampleMergedToProcess[proj][ref_strain][strain]
                    refParams = refInfoParams[ref_strain]
                    dataObj.single_run_threads = int(optionTools.threads / float(processNum))
                    if args.command == 'mapping':
                        from mapping import mapping
                        # mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
                        pool.apply_async(mapping, (dataObj, minimap2Params, refParams, dirSpec, dataObj.single_run_threads, args.correction, args.juncCombSup))
                    if args.command == 'collapse':
                        from collapse import collapse
                        # collapse(dataObj=dataObj, collapseParams=collapseParams, refParams=refParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
                        pool.apply_async(collapse, (dataObj, collapseParams, refParams, dirSpec, dataObj.single_run_threads))
                    if args.command == 'refine':
                        from refine import refineJunc
                        pool.apply_async(refineJunc, (dataObj, refParams, dirSpec, args.refine, args.adjust))
                    if args.command == 'pu_filter':
                        from iso_pu import iso_pu1
                        # iso_pu1(dataObj, dirSpec, refParams, hqIsoParams)
                        pool.apply_async(iso_pu1, (dataObj, dirSpec, refParams, hqIsoParams))
                    if args.command == 'find_as':
                        # from find_charaterize_as_functions import *
                        from identify_as import identify_as
                        # identify_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                        pool.apply_async(identify_as, (dataObj, refParams, dirSpec, hqIsoParams, args))
                    if args.command == 'visual_as':
                        from visual_as import visual_as
                        targetGenes = args.genes
                        # visual_as(dataObj=dataObj, targetGenes=targetGenes, refParams=refParams, dirSpec=dirSpec)
                        pool.apply_async(visual_as, (dataObj, targetGenes, refParams, dirSpec))
                    if args.command == 'asas':
                        from asas import asas, ase
                        asas(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                        # pool.apply_async(asas, (dataObj, refParams, dirSpec))
                        if args.ase:
                            ase(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, refFa=args.refFa,
                                altFa=args.altFa, useFreeBayes=args.useFreebayes)
                            # pool.apply_async(ase, (dataObj, refParams, dirSpec, args.refFa, args.altFa, args.useFreebayes))
                    if args.command == 'palen_as':
                        from palen_as import palen_as
                        # palen_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, filterByCount=args.pa_support, dataToProcess=dataToProcess, confidentPac=args.confidentPac)
                        pool.apply_async(palen_as, (dataObj, refParams, dirSpec, args.pa_support, dataToProcess, args.confidentPac))
        pool.close()
        pool.join()
        if args.command == 'visual_as':
            from visual_as import visual_as_merge
            targetGenes = args.genes
            # plot the gene structure for all samples
            for proj in sampleMergedToProcess:
                for ref_strain in sampleMergedToProcess[proj]:
                    dataToProcess = sampleMergedToProcess[proj][ref_strain].values()
                    visual_as_merge(dataToProcess=dataToProcess, targetGenes=targetGenes, refParams=refInfoParams[ref_strain], dirSpec=dirSpec)
        if args.command == 'diff_as':
            from diff_as import diff_as
            # diff_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            compCond = args.compCond
            diff_as(sampleMergedToProcess, compCondFile=compCond, dirSpec=dirSpec, sampleMerged=True, args=args, optionTools=optionTools)
        if args.command == 'go':
            from go import go
            go(args, optionTools=optionTools, dirSpec=dirSpec)
        if args.command == 'report':
            from report import report
            # report(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            for proj in sampleMergedToProcess:
                for ref_strain in sampleMergedToProcess[proj]:
                    dataToProcess = sampleMergedToProcess[proj][ref_strain].values()
                    report(dataToProcess=dataToProcess, refInfoParams=refInfoParams, dirSpec=dirSpec, optionTools=optionTools, args=args)
    else:
        pool = MyPool(processes=len(dataToProcess))
        for dataObj in dataToProcess:
            refParams = refInfoParams[dataObj.ref_strain]
            dataObj.single_run_threads = int(optionTools.threads / float(len(dataToProcess)))
            if args.command == 'mapping':
                from mapping import mapping
                # mapping(dataObj=dataObj, minimap2Params=minimap2Params, refParams=refParams, dirSpec=dirSpec, threads=dataObj.single_run_threads, juncCombSup=args.juncCombSup)
                pool.apply_async(mapping, (dataObj, minimap2Params, refParams, dirSpec, dataObj.single_run_threads, args.correction, args.juncCombSup))
            if args.command == 'collapse':
                from collapse import collapse
                # collapse(dataObj=dataObj, collapseParams=collapseParams, refParams=refParams, dirSpec=dirSpec, threads=dataObj.single_run_threads)
                pool.apply_async(collapse, (dataObj, collapseParams, refParams, dirSpec, dataObj.single_run_threads))
            if args.command == 'refine':
                from refine import refineJunc
                pool.apply_async(refineJunc, (dataObj, refParams, dirSpec, args.refine, args.adjust))
            if args.command == 'pu_filter':
                from iso_pu import iso_pu1
                # iso_pu1(dataObj, dirSpec, refParams, hqIsoParams)
                pool.apply_async(iso_pu1, (dataObj, dirSpec, refParams, hqIsoParams))
            if args.command == 'find_as':
                # from find_charaterize_as_functions import *
                from identify_as import identify_as
                # identify_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                pool.apply_async(identify_as, (dataObj, refParams, dirSpec, hqIsoParams, args))
            if args.command == 'visual_as':
                from visual_as import visual_as
                targetGenes = args.genes
                # visual_as(dataObj=dataObj, targetGenes=targetGenes, refParams=refParams, dirSpec=dirSpec)
                pool.apply_async(visual_as, (dataObj, targetGenes, refParams, dirSpec))
            if args.command == 'asas':
                from asas import asas, ase
                asas(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
                # pool.apply_async(asas, (dataObj, refParams, dirSpec))
                if args.ase:
                    ase(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, refFa=args.refFa, altFa=args.altFa,
                        useFreeBayes=args.useFreebayes)
                    # pool.apply_async(ase, (dataObj, refParams, dirSpec, args.refFa, args.altFa, args.useFreebayes))

            if args.command == 'palen_as':
                from palen_as import palen_as
                palen_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec, filterByCount=args.pa_support, dataToProcess=dataToProcess, confidentPac=args.confidentPac)
                # pool.apply_async(palen_as, (dataObj, refParams, dirSpec, args.pa_support, dataToProcess))

        pool.close()
        pool.join()
        if args.command == 'visual_as':
            from visual_as import visual_as_merge
            targetGenes = args.genes
            # plot the gene structure for all samples
            for ref_strain in refInfoParams:
                visual_as_merge(dataToProcess=dataToProcess, targetGenes=targetGenes, refParams=refInfoParams[ref_strain], dirSpec=dirSpec)
        if args.command == 'diff_as':
            from diff_as import diff_as
            # diff_as(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            compCond = args.compCond
            diff_as(dataToProcess, compCondFile=compCond, dirSpec=dirSpec, sampleMerged=False, args=args, optionTools=optionTools, refInfoParams=refInfoParams)
        if args.command == 'go':
            from go import go
            go(args, optionTools=optionTools, dirSpec=dirSpec)
        if args.command == 'report':
            from report import report
            # report(dataObj=dataObj, refParams=refParams, dirSpec=dirSpec)
            report(dataToProcess=dataToProcess, refInfoParams=refInfoParams, dirSpec=dirSpec, optionTools=optionTools, args=args)


def getDefaultOptions(subargs):
    default_args = {}
    for action in subargs._actions:
        if isinstance(action, argparse._StoreAction) or isinstance(action, argparse._StoreConstAction) or \
                isinstance(action, argparse._StoreTrueAction) or isinstance(action, argparse._StoreFalseAction):
            default_args[action.dest] = action.default
    return default_args


def iflas(args, changed_args):
    if not validateFile(args.default_cfg):
        print "The configuration file associated to the -cfg option is invalid. Please provide corrected one!"
        return
    # defaultCfg = Config(args.default_cfg)
    defaultCfg = Config1(args.default_cfg, args, changed_args)
    dataToProcess = defaultCfg.dataToProcess
    refInfoParams = defaultCfg.refInfoParams
    ccsParams = defaultCfg.ccsParams
    minimap2Params = defaultCfg.minimap2Params
    collapseParams = defaultCfg.collapseParams
    hqIsoParams = defaultCfg.hqIsoParams
    optionTools = defaultCfg.optionTools
    dirSpec = defaultCfg.dirParams
    print minimap2Params, hqIsoParams, args
    return
    for refStrain in refInfoParams:
        refParams = refInfoParams[refStrain]
        initRefSetting(refParams=refParams, dirSpec=dirSpec)
    initSysResourceSetting(optionTools=optionTools)

    pybedtools.set_tempdir(dirSpec.tmp_dir)
    if args.command == "all":
        oneCommandRun(args, dataToProcess, refInfoParams, dirSpec, ccsParams, minimap2Params, collapseParams, optionTools)
    else:
        splitCommandRun(args, dataToProcess, refInfoParams, dirSpec, ccsParams, minimap2Params, collapseParams, hqIsoParams, optionTools)
    pybedtools.cleanup(remove_all=True)

if __name__ == "__main__":
    USAGE = ' iFLAS: integrated Full Length Alternative Splicing analysis '

    parent_parser = argparse.ArgumentParser(add_help=False)
    parent_parser.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parent_parser.add_argument('-merge', dest="merge", action="store_true", default=False, help="Merge all samples from a same strain.")

    go_parser = argparse.ArgumentParser(add_help=False)
    go_parser.add_argument('-bg', dest="gene2goFile", type=str, default=None, help="The mapping file between gene and go term used for GO enrichment analysis.")
    go_parser.add_argument('-cutoff', dest="cutoff", type=float, default=0.05, help="The cutoff used to filter the output. Default: 0.05")
    go_parser.add_argument('-filterBy', dest="filterBy", type=str, choices=["pvalue", "p.adjust"], default="p.adjust", help="The value used to filter. Default: p.adjust.")
    go_parser.add_argument('-showCategory', dest="showCategory", type=int, default=20, help="The number of items to show off. Default: 20.")

    parser = argparse.ArgumentParser(usage='%(prog)s command [options]', description=USAGE)
    subparsers = parser.add_subparsers(title='command', metavar='', dest='command', prog=parser.prog)
    parser_preprocess = subparsers.add_parser('preproc', parents=[parent_parser], help='Pre-process the raw PacBio/NanoPore/NGS data. When TGS and NGS data both are provide, This step will use fmlrc2 to correct the TGS read with the information in NGS', usage='%(prog)s [options]')
    # parser_preprocess.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")

    parser_mapping = subparsers.add_parser('mapping', parents=[parent_parser], help='Mapping the TGS/NGS reads to the reference genome with minimap2', usage='%(prog)s [options]')
    # parser_mapping.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parser_mapping.add_argument('-c', dest="correction", action="store_true", default=False, help="Correct the flnc reads with fmlrc2.")
    parser_mapping.add_argument('-jcs', dest="juncCombSup", type=int, default=2, help="The number of junction combination supported by flnc reads. Default: 2.")

    parser_collapse = subparsers.add_parser('collapse', parents=[parent_parser], help='Collapse corrected reads into high-confidence isoforms', usage='%(prog)s [options]')
    # parser_collapse.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")

    parser_refine = subparsers.add_parser('refine', parents=[parent_parser], help='Refine the splice junction with the information in short reads', usage='%(prog)s [options]')
    # parser_refine.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parser_refine.add_argument('-adjust', dest="adjust", action="store_true", default=False, help="Adjust the strand orient by the information of junctions.")
    parser_refine.add_argument('-refine', dest="refine", action="store_true", default=False, help="Refine the junction position by the reads support.")

    parser_pu = subparsers.add_parser('pu_filter', parents=[parent_parser], help='Filter Low-Quality novel isoforms using a PU-learning based method', usage='%(prog)s [options]')
    parser_pu.add_argument('-filter_score', dest="filter_score", type=float, default=0.5, help="The PU-score that used to filter out low quality novel isoforms. Default: 0.5.")
    parser_pu.add_argument('-draw_auc', dest="draw_auc", action="store_true", default=False, help="To draw the AUC plot or not.")
    parser_pu.add_argument('-pos_fl_cov', dest="pos_fl_coverage", type=int, default=2, help="The minimal coverage that get the positive annotated isoforms. Default: 2.")
    parser_pu.add_argument('-pos_min_junc_rpkm', dest="pos_min_junc_rpkm", type=float, default=0.05, help="The minimal rpkm value of the junctions needed for a postive annotated isoform. Default: 0.05.")
    parser_pu.add_argument('-select_best_model', dest="select_best_model", action="store_true", default=False, help="Select the best model. If not, iFLAS will use GB.")
    parser_pu.add_argument('-auto_filter_score', dest="auto_filter_score", action="store_true", default=False, help="Auto determine pu_score when 'select_best_model' is selected.")

    parser_findAS = subparsers.add_parser('find_as', parents=[parent_parser], help='Identify alternative splicing(AS) type from high-confidence isoforms. Four common AS type are included: intron retention, exon skipping, alternative 3 end splicing and alternative 5 end splicing', usage='%(prog)s [options]')
    # parser_findAS.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parser_findAS.add_argument('-pa_rpkm', dest="paRPKM", type=float, default=0, help="Filter the pa cluster by RPKM(PAC). Default: 0.")
    parser_findAS.add_argument('-pa_sup', dest="pa_support", type=int, default=5, help="Filter the pa cluster by RPKM(PAC). Default: 5.")
    parser_findAS.add_argument('-conf_pa', dest="confidentPa", default=None, help="The confident PA file used for filtering the results.")

    parser_visualAS = subparsers.add_parser('visual_as', parents=[parent_parser], help='Visualize the specific gene structure with details including isoform mapping, short reads coverage and AS types identified', usage='%(prog)s [options]')
    # parser_visualAS.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parser_visualAS.add_argument('-g', dest="genes", type=str, help="The gene list separated by comma or a single file contain genes one per line used for visualization.")

    # parser_rankAS = subparsers.add_parser('rank_iso', parents=[parent_parser], help='Score the isoform by the produce of each inclusion/exclusion ratio in that isoform, and rank all the isoforms from high to low', usage='%(prog)s [options]')
    # # parser_rankAS.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    # parser_rankAS.add_argument('-coding', dest="coding", action="store_true", default=False, help="Filter the isoforms by coding or not.")
    # parser_rankAS.add_argument('-min_tpm', dest="min_tpm", type=float, default=0, help="Filter the isoforms by minimal TPM value.")
    # parser_rankAS.add_argument('-reads_freq', dest="reads_freq", type=float, default=0, help="Filter isoforms by the frequency of reads.")
    # parser_rankAS.add_argument('-read_support', dest="read_support", type=int, default=0, help="Filter isoforms by the count of reads.")

    parser_alleleAS = subparsers.add_parser('asas', parents=[parent_parser, go_parser], help='Identify allele-specific AS', usage='%(prog)s [options]')
    # parser_alleleAS.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parser_alleleAS.add_argument('-ase', dest="ase", action="store_true", default=False, help="Whether to Carry out ASE analysis.")
    parser_alleleAS.add_argument('-ref_fa', dest="refFa", default=None, help="The reference fasta file used to be quantified in ASE.")
    parser_alleleAS.add_argument('-alt_fa', dest="altFa", default=None, help="The alternative fasta file used to be quantified in ASE.")
    parser_alleleAS.add_argument('-fbs', dest="useFreebayes", action="store_true", default=False, help="Call the heterozygosity SNPs with freebayes in ASE.")

    parser_palenAS = subparsers.add_parser('palen_as', parents=[parent_parser, go_parser], help='Identify functional poly(A) tail length related to AS', usage='%(prog)s [options]')
    # parser_palenAS.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parser_palenAS.add_argument('-pa_sup', dest="pa_support", type=int, default=10, help="The pa cluster coverage supported by flnc reads. Default: 10.")
    parser_palenAS.add_argument('-conf_pac', dest="confidentPac", default=None, help="The confident PAC file used for filtering the results")

    parser_diffAS = subparsers.add_parser('diff_as', parents=[parent_parser, go_parser], help='Carry out differential AS ananlysis among conditions', usage='%(prog)s [options]')
    # parser_diffAS.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parser_diffAS.add_argument('-d', dest="compCond", type=str, help="The condition file used to detect differential AS between samples.")
    parser_diffAS.add_argument('-go', dest="go", action="store_true", default=False, help="Perform GO enrichment analysis for DSGs between samples.")
    parser_diffAS.add_argument('-pu_filter', dest="pu_filter", action="store_true", default=False, help="Get differentially spliced gene only from ranked isoforms.")
    # parser_diffAS.add_argument('-bg', dest="gene2goFile", type=str, default=None, help="The mapping file between gene and go term used for GO enrichment analysis or defined in config file at 'optionTools' section.")
    # parser_diffAS.add_argument('-cutoff', dest="cutoff", type=float, default=0.05, help="The cutoff used to filter the output. Default: 0.05")
    # parser_diffAS.add_argument('-filterBy', dest="filterBy", type=str, choices=["pvalue", "p.adjust"], default="p.adjust", help="The value used to filter. Default: p.adjust.")
    # parser_diffAS.add_argument('-showCategory', dest="showCategory", type=int, default=20, help="The number of items to show off. Default: 20.")

    parser_goAS = subparsers.add_parser('go', parents=[parent_parser, go_parser], help='Perform GO enrichment analysis and plot results for the specified gene set or multiple gene sets', usage='%(prog)s [options]')
    # parser_goAS.add_argument('-cfg', dest="default_cfg", type=str, help="The config file used for init setting.")
    parser_goAS.add_argument('-tg', dest="targetGeneFile", type=str, help="The target gene file or file list separated by comma used for GO enrichment analysis.")
    # parser_goAS.add_argument('-bg', dest="gene2goFile", type=str, default=None, help="The mapping file between gene and go term used for GO enrichment analysis.")
    # parser_goAS.add_argument('-cutoff', dest="cutoff", type=float, default=0.05, help="The cutoff used to filter the output. Default: 0.05")
    # parser_goAS.add_argument('-filterBy', dest="filterBy", type=str, choices=["pvalue", "p.adjust"], default="p.adjust", help="The value used to filter. Default: p.adjust.")
    # parser_goAS.add_argument('-showCategory', dest="showCategory", type=int, default=20, help="The number of items to show off. Default: 20.")
    parser_goAS.add_argument('-s', dest="sampleName", type=str, help="The sample name used plot the track, multi-sample should be separated by commma used for GO enrichment analysis.")
    parser_goAS.add_argument('-o', dest="outName", type=str, default="goEnrichment", help="The prefix of the GO enrichment output file.")

    parser_report = subparsers.add_parser('report', parents=[parent_parser], help='Automatic detect the plots generated in each step, and merge them into a report file', usage='%(prog)s [options]')
    # parser_report.add_argument('-cfg', dest="default_cfg", help="The config file used for init setting.")
    parser_report.add_argument('-all', dest="all", action="store_true", default=False, help="Generate all the plots.")
    parser_report.add_argument('-basic', dest="basic", action="store_true", default=False, help="Generate basic information plots.")
    parser_report.add_argument('-asp', dest="asp", action="store_true", default=False, help="Generate AS pattern plots.")
    parser_report.add_argument('-geneStruc', dest="geneStruc", action="store_true", default=False, help="Generate gene structure with AS events.")
    # parser_report.add_argument('-novelIso', dest="novelIso", action="store_true", default=False, help="Generate identity curve for novel isoforms.")
    parser_report.add_argument('-asas', dest="asas", action="store_true", default=False, help="Generate allele-specific AS events.")
    parser_report.add_argument('-palen', dest="palen", action="store_true", default=False, help="Generate AS events related differential poly(A) tail length.")
    parser_report.add_argument('-diff', dest="diff", action="store_true", default=False, help="Generate the statistics for differential spliced events.")
    parser_report.add_argument('-html', dest="html", action="store_true", default=False, help="Generate the html report for the results have been generated.")
    # parser_report.add_argument('-bg', dest="gene2goFile", type=str, default=None, help="The mapping file between gene and go term used for GO enrichment analysis.")
    # parser_report.add_argument('-cutoff', dest="cutoff", type=float, default=0.05, help="The cutoff used to filter the output. Default: 0.05")
    # parser_report.add_argument('-filterBy', dest="filterBy", type=str, choices=["pvalue", "p.adjust"], default="p.adjust", help="The value used to filter. Default: p.adjust.")
    # parser_report.add_argument('-showCategory', dest="showCategory", type=int, default=20, help="The number of items to show off. Default: 20.")

    parser_all = subparsers.add_parser('all', parents=[parent_parser, go_parser], help="Dynamically perform all analysis with the setting!")
    parser_all.add_argument('-c', dest="correction", action="store_true", default=False, help="Correct the flnc reads with fmlrc2.")
    parser_all.add_argument('-g', dest="genes", type=str, help="The gene list separated by comma or a single file contain genes one per line used for visualization.")
    parser_all.add_argument('-d', dest="compCond", type=str, help="The condition file used to detect differential AS between samples.")
    parser_all.add_argument('-tg', dest="targetGeneFile", type=str, help="The target gene file or file list separated by comma used for GO enrichment analysis.")
    # parser_all.add_argument('-bg', dest="gene2goFile", type=str, default=None, help="The mapping file between gene and go term used for GO enrichment analysis.")
    parser_all.add_argument('-s', dest="sampleName", type=str, help="The sample name used plot the track, multi-sample should be separated by commma used for GO enrichment analysis.")
    parser_all.add_argument('-o', dest="outName", type=str, default="goEnrichment", help="The prefix of the GO enrichment output file.")

    if len(sys.argv) <= 2:
        sys.argv.append('-h')
    args = parser.parse_args()
    default_args = {}
    for subargs in [parent_parser, go_parser, parser_preprocess, parser_mapping, parser_collapse, parser_refine,
                    parser_pu, parser_findAS, parser_visualAS, parser_alleleAS, parser_diffAS, parser_goAS, parser_report]:
        default_args.update(getDefaultOptions(subargs))
    # print default_args
    changed_args = {k: v for k, v in vars(args).items() if k in default_args and v != default_args[k]}
    # print args, changed_args
    iflas(args, changed_args)
