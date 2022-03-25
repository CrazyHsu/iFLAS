python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py preproc -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py mapping -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg -c -jcs 2
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py collapse -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py refine -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg -refine
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py find_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py visual_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg -g Zm00001d050245
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py rank_iso -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py diff_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg -d /data/CrazyHsu_data/Projects/isoseq/iFLAS_toolkit_20210510/iFLAS_toolkit/test_data/wangbo_compCond.txt
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py go -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg -tg test_data/targetGene1.txt,test_data/targetGene2.txt,test_data/targetGene3.txt -bg /data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.50/gene2goId.txt -s sample1,sample2,sample3 -filterby pvalue -cutoff 0.1
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py report -basic -asp -geneStruc -diff -html -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wangbo.cfg
