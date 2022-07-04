run_dir="/data/CrazyHsu_data/Projects/isoseq/iFLAS_20220704"
python ${run_dir}/iflas.py preproc -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg
python ${run_dir}/iflas.py mapping -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg -c -jcs 2
python ${run_dir}/iflas.py collapse -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg
python ${run_dir}/iflas.py refine -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg -refine
python ${run_dir}/iflas.py find_as -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg
python ${run_dir}/iflas.py visual_as -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg -g Zm00001d050245
python ${run_dir}/iflas.py rank_iso -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg
python ${run_dir}/iflas.py diff_as -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg -d ${run_dir}_toolkit_20210510/iFLAS_toolkit/test_data/wangbo_compCond.txt
python ${run_dir}/iflas.py go -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg -tg test_data/targetGene1.txt,test_data/targetGene2.txt,test_data/targetGene3.txt -bg /data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.50/gene2goId.txt -s sample1,sample2,sample3 -filterby pvalue -cutoff 0.1
python ${run_dir}/iflas.py report -basic -asp -geneStruc -diff -html -cfg ${run_dir}/test_data/test_pacbio_wangbo.cfg
