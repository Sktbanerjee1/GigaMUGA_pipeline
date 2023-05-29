# clear env
rm(list = ls())


source("src/build_packages.R")

# arguments

args_list <- list(
  optparse::make_option(
    c("-i", "--input"),
    default = "data/",
    help = "Path to the folder with the `FinalReport` file"
  ),
  optparse::make_option(
    c("-o","--out_dir"),
    default = "results/",
    help = "folder to save the output"
  )
)


arg_parser <- optparse::OptionParser(
  option_list = args_list
)

args <- optparse::parse_args(
  arg_parser
)

if (is.null(args$input)) {
  optparse::print_help(arg_parser)
  stop("--input not supplied!")
}
# arguments end
#-------------

# -------
# check output directory
out_dir <- normalizePath(
  paste0(getwd(),"/",args$out_dir)
)

logger::log_info(
  "searching out_dir:\n",
  out_dir
)

if(!dir.exists(out_dir)){
  dir.create(
    out_dir
  )
  logger::log_warn("Not Found!\ncreated", out_dir)
} else{
  logger::log_info("Found..")
}
# check output directory end
#------------------

# find FinalReport File (geneseek)
# parse FinalReport data into a dataframe
logger::log_info("Searching for FinalReport file")
file_list <- list.files(args$input)
report_file <- file_list[grep(".FinalReport", file_list)]
report_files_path <- paste0(args$input, report_file)
# check number of FinalReport files found
if(length(report_files_path)==0){
  logger::log_error("No .FinalReport file was found")
  gigamuga_data <- NULL
  stop()
} else if (length(report_files_path)>1){
  logger::log_warn("More than one final FinalReport file found")
  print(report_files_path)
  gigamuga_data <- plyr::ldply(report_files_path, function(fp){
    # read FinalReport.txt skipping metadata lines
    # metadata lines = 9 in test data
    readr::read_delim(fp, delim = "\t", col_names = TRUE, skip = 9)
  }, .progress = "text")
} else {
  logger::log_info(paste("Found file", report_files_path))
  gigamuga_data <- readr::read_delim(report_files_path, delim = "\t", col_names = TRUE, skip=9)
}
# parse FinalReport file end
#--------------------------------

# write parsed FinalReport file
data.table::fwrite(gigamuga_data, paste0(args$out_dir,"parsed_FinalReport.csv"))
