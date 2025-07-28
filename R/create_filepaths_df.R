create_filepaths_df <- function(base_path, project_nr, TITLE) {
    filepaths <- data.frame(
        name = c("assembly", "coverage", "busco", "taxonomy", "assembly_info", "fasta_folder", "bed_folder", "bam_file", "telomere_left", "telomere_right"),
        path = c(
            file.path(base_path, paste0("flye_output/", project_nr, TITLE, "_flye/assembly.fasta")),
            file.path(base_path, paste0("coverm_output/", project_nr, TITLE, "_mean_cov/mean_coverage.tsv")),
            file.path(base_path, paste0("busco_output/", project_nr, TITLE, "_busco/run_ascomycota_odb12/full_table.tsv")),
            file.path(base_path, paste0("diamond_output/", project_nr, TITLE, "_taxonomy_prots/", TITLE, "_taxonomy_prots")),
            file.path(base_path, paste0("flye_output/", project_nr, TITLE, "_flye/assembly_info.txt")),
            file.path(base_path, paste0("binning_output/", project_nr, TITLE, "_binned/contigs")),
            file.path(base_path, paste0("binning_output/", project_nr, TITLE, "_binned/beds")),
            file.path(base_path, paste0("coverm_output/", project_nr, TITLE, "_mean_cov/bam_cache/assembly.fasta.",project_nr , TITLE, ".hifireads.fastq.gz.bam")),
            file.path(base_path, paste0("telfinder_output/",project_nr, TITLE, "_telomeres/NCR_left_score.txt")),
            file.path(base_path, paste0("telfinder_output/",project_nr, TITLE, "_telomeres/NCR_right_score.txt"))
        ),
        stringsAsFactors = FALSE
    )
    return(filepaths)
}
