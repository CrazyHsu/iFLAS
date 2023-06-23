#run_dir="/data/CrazyHsu_data/Projects/isoseq/iFLAS_20220704"
python iflas.py preproc -cfg test_data/test_pacbio_wangbo.cfg
python iflas.py mapping -cfg test_data/test_pacbio_wangbo.cfg -c -jcs 2
python iflas.py collapse -cfg test_data/test_pacbio_wangbo.cfg
python iflas.py refine -cfg test_data/test_pacbio_wangbo.cfg -refine
python iflas.py pu_filter -cfg test_data/test_pacbio_wangbo.cfg -draw_auc -select_best_model -auto_filter_score
python iflas.py find_as -cfg test_data/test_pacbio_wangbo.cfg
python iflas.py visual_as -cfg test_data/test_pacbio_wangbo.cfg -g Zm00001d050245
#python iflas.py rank_iso -cfg test_data/test_pacbio_wangbo.cfg
python iflas.py diff_as -cfg test_data/test_pacbio_wangbo.cfg -d test_data/wangbo_compCond.txt
python iflas.py go -cfg test_data/test_pacbio_wangbo.cfg -tg test_data/targetGene1.txt,test_data/targetGene2.txt,test_data/targetGene3.txt -bg /data/CrazyHsu_data/Genomes/maize/Zea_mays.AGPv4.50/gene2goId.txt -s sample1,sample2,sample3 -filterby pvalue -cutoff 0.1
python iflas.py report -basic -asp -geneStruc -diff -html -cfg test_data/test_pacbio_wangbo.cfg