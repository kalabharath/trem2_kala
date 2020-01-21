export top_dir=/home/kalabharath/PycharmProjects/trem2_kala/analys
export analys_dir=$top_dir/analys/testsystem
export mod_dir=$analys_dir/GSMs_cl-1
export name="testsys_-1"

cp $top_dir/analys/density.txt $mod_dir
cp $mod_dir/sample_A/sample_A.txt $mod_dir/Scores_A.txt
cp $mod_dir/sample_B/sample_B.txt $mod_dir/Scores_B.txt

ls -lta $analys_dir/GSMs_cl-1/sample_A | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_A_cluster-1_random.dat
ls -lta $analys_dir/GSMs_cl-1/sample_B | awk '{print $9}' | grep 'rmf3' > $analys_dir/selected_models_B_cluster-1_random.dat

/usr/bin/python /home/kalabharath/PycharmProjects/trem2_kala/analys/imp-sampcon/pyext/src/Master_Sampling_Exhaustiveness_Analysis.py --sysname $name --path $mod_dir --mode cpu_omp --cores 4 --align --density density.txt --gridsize 1.0 --gnuplot --scoreA /Scores_A.txt --scoreB /Scores_B.txt --cluster_threshold 18.00 > addBarrier_clustering_1.log
