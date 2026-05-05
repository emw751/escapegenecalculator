## Helps you setup a config file.
## The goal of this script is to eliminate the dependency of escapegenecalculator on
## using filenames to store metadata.

wd = dirname(dirname(this.path::here()))  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'list_files', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(

    make_option(c("-i", "--input-dir"), default="data/berletch-spleen",
                metavar="data/berletch-spleen", type="character",
                help="set the base directory"),

    make_option(c("-o", "--output-dir"), default="output-2",
                metavar="output-2", type="character",
                help="data will be output in this folder inside the input-dir"),

    make_option(c("-m", "--mat-mouse-strain"), default="Mus_musculus",
                metavar="Mus_musculus", type="character",
                help="set the maternal mouse strain"),

    make_option(c("-p", "--pat-mouse-strain"), default="Mus_musculus_casteij",
                metavar="Mus_castaneus", type="character",
                help="set the paternal mouse strain"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting = opt[['troubleshooting']]

in_dir = file.path(wd, opt[['input-dir']])
out_dir = file.path(in_dir, opt[['output-dir']])
mat_mouse_strain = opt[['mat-mouse-strain']]
pat_mouse_strain = opt[['pat-mouse-strain']]

# required defaults
mat_reads_dir = file.path(in_dir, "input", "mat_reads")
pat_reads_dir = file.path(in_dir, "input", "pat_reads")
rpkms_dir = file.path(in_dir, "input", "rpkms")
default_mouse_gender = "female"


# Start Log
start_time = Sys.time()
log <- log_open(paste0("generate_config-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))
if (!troubleshooting) {
    log_print(paste(Sys.time(), 'in_dir:', in_dir))
    log_print(paste(Sys.time(), 'out_dir:', out_dir))
} else {
    log_print(paste(Sys.time(), 'troubleshooting:', troubleshooting))
}


# ----------------------------------------------------------------------
# Generate Config File

mat_reads_files = list_files(mat_reads_dir, recursive=FALSE, full_name=FALSE)
pat_reads_files = list_files(pat_reads_dir, recursive=FALSE, full_name=FALSE)
rpkms_files = list_files(rpkms_dir, recursive=FALSE, full_name=FALSE)

log_print(paste(Sys.time(), 'mat_reads files found:', length(mat_reads_files)))
log_print(paste(Sys.time(), 'pat_reads files found:', length(pat_reads_files)))
log_print(paste(Sys.time(), 'rpkms files found:', length(rpkms_files)))


# Note: cbind can lead to some strange edge effects
# When no files are found, the column is missing entirely
# When one file is found, the entire column is filled with that one file
# Will think about how to solve this later
df = as.data.frame(do.call(cbind,
    list(mat_reads_filenames=mat_reads_files,
         pat_reads_filenames=pat_reads_files,
         rpkms_filenames=rpkms_files))
)

if (length(rpkms_files)==0) {
    df['rpkms_filenames'] = NA
}

df['mouse_id'] = paste('mouse_', seq(1, nrow(df), by=1), sep='')  # generate mouse_id if needed
df['mouse_gender'] = default_mouse_gender  # see above
df['mat_mouse_strain'] = mat_mouse_strain
df['pat_mouse_strain'] = pat_mouse_strain


# save
if (!troubleshooting) {

    log_print(paste(Sys.time(), "Writing data..."))

    if (!dir.exists(file.path(out_dir))) {
        dir.create(file.path(out_dir), recursive=TRUE)
    }

    write.table(
        df,
        file = file.path(out_dir, 'config.csv'),
        row.names = FALSE,
        sep = ','
    )

    log_print(paste(Sys.time(), "Remember to check the config file before using!"))
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
