python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py preproc -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py mapping -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg -c -jcs 2
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py collapse -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py refine -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg -refine
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py find_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py visual_as -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg -g Zm00001d050245
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py rank_iso -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py go -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg -tg targetGene1.txt,targetGene2.txt,targetGene3.txt -bg /data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.50/gene2goId.txt -s sample1,sample2,sample3
python /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/iflas.py report -cfg /data/CrazyHsu_data/Projects/isoseq/iFLAS_simple/test_data/test_pacbio_wanglab.cfg
