analysis_dir="../EPICC_analysis_publication_test"


cp "$analysis_dir/analysis/ATAC/MEGABULKS/datasets/edger_summary/final_reccurence_data.rds" ./data/
cp "$analysis_dir/created_datasets/genotyping_estimates_per_sample.rds" ./data/

cp "$analysis_dir/created_datasets/website/figures/heatmap_scaas_in_driver_genes_martincorena_and_intogen_all.rds" ./data/pltA.rds
cp "$analysis_dir/created_datasets/website/figures/heatmap_recurrent_epigenetic_changes_few.rds" ./data/pltB.rds

cp "$analysis_dir/created_datasets/lowpass_cn_calls.rds" ./data/other/
cp "$analysis_dir/created_datasets/sequenza_cna_data.rds" ./data/other/
