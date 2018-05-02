

dirs <- c("inst/extdata/moderns", "inst/extdata/ancients")
out <- prepare_archaic(dirs, delete = TRUE)
model <- model_archaic(out, K = 2, output_dir = "inst/extdata/archaic_results")
viz <- plot_archaic(model, output_dir = "inst/extdata/archaic_results",
                    background = "modern")


