# Source this file to load all callable R helpers from r_functions.

source_frames <- Filter(
  function(frame) !is.null(frame$ofile),
  sys.frames()
)
loader_path <- if (length(source_frames) > 0) {
  tail(source_frames, 1)[[1]]$ofile
} else {
  "r_functions/load_functions.R"
}
function_dir <- dirname(normalizePath(loader_path, mustWork = TRUE))

source(file.path(function_dir, "semi_manual_annotation.R"))
