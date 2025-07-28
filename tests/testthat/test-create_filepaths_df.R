test_that("create_filepaths_df returns correct data frame", {
    base_path <- "/base_output"
    project_nr <- "pr_01_"
    TITLE <- "001"

    df <- create_filepaths_df(base_path, project_nr, TITLE)

    # Check output is a data frame
    expect_s3_class(df, "data.frame")

    # Check columns
    expect_true(all(c("name", "path") %in% colnames(df)))

    # Check the length is as expected (10 items)
    expect_equal(nrow(df), 10)

    # Check some paths contain the inputs
    expect_true(any(grepl(project_nr, df$path)))
    expect_true(any(grepl(TITLE, df$path)))
    expect_true(any(grepl(base_path, df$path)))

    # Check specific path is constructed correctly
    expected_assembly_path <- file.path(base_path, paste0("flye_output/", project_nr, TITLE, "_flye/assembly.fasta"))
    expect_equal(df$path[df$name == "assembly"], expected_assembly_path)
})

