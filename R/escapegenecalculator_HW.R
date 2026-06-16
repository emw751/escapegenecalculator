## This pipeline is based on Berletch and Disteche's 2015 paper on XCI escape in mouse tissues
## It takes SNP-specific mapped read counts as input,
## normalizes, computes confidence intervals, then outputs a list of escape genes.
## 
## Usage:
## Rscript R/escapegenecalculator.R -i data/berletch-spleen

wd = dirname(this.path::here())  # wd = '~/github/R/escapegenecalculator'
library('optparse')
library('logr')
import::from(file.path(wd, 'R', 'tools', 'df_tools.R'),
    'coalesce1', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'file_io.R'),
    'read_csv_or_tsv', .character_only=TRUE)
import::from(file.path(wd, 'R', 'tools', 'list_tools.R'),
    'items_in_a_not_b', .character_only=TRUE)
import::from(file.path(wd, 'R', 'functions', 'preprocessing.R'),
    'get_single_reads', .character_only=TRUE)
import::from(file.path(wd, 'R', 'configs', 'columns.R'),
    'escape_gene_cols', 'x_read_output_cols', .character_only=TRUE)


# ----------------------------------------------------------------------
# Pre-script settings

# args
option_list = list(

    make_option(c("-i", "--input-dir"), default="data/berletch-spleen",
                metavar="data/berletch-spleen", type="character",
                help="set the directory"),

    make_option(c("-o", "--output-subdir"), default="output-2",
                metavar="output-2", type="character",
                help="useful for running multiple scripts on the same dataset"),

    make_option(c("-z", "--zscore-threshold"), default=0.99,
                metavar="0.99", type="double",
                help="was 0.975 in Zack's version, but Berletch's paper requires 0.99"),

    make_option(c("-t", "--troubleshooting"), default=FALSE, action="store_true",
                metavar="FALSE", type="logical",
                help="enable if troubleshooting to prevent overwriting your files")
)
opt_parser = OptionParser(option_list=option_list)
opt = parse_args(opt_parser)
troubleshooting <- opt[['troubleshooting']]

# for readability downstream
base_dir = file.path(wd, opt[['input-dir']])
in_dir = file.path(base_dir, 'input')
out_dir = file.path(base_dir, opt[['output-subdir']])
zscore = qnorm(opt[['zscore-threshold']])  # 2.326 if zscore_threshold=0.99

# optional commands
gene_id_col='gene_id'  # could use 'locus' for Disteche's data

# in case num_total_reads is unavailable, use this to estimate the num_total_reads
# In Berletch's 2015 paper, they got 9258991 snp-specific reads from 88842032 total reads
estimated_pct_snp_reads=0.104  # 9258991/88842032


# Start Log
start_time = Sys.time()
log <- log_open(paste0("escapegenecalculator-",
                       strftime(start_time, format="%Y%m%d_%H%M%S"), '.log'))
log_print(paste('Script started at:', start_time))
if (!troubleshooting) {
    log_print(paste(Sys.time(), 'base_dir:', base_dir))
    log_print(paste(Sys.time(), 'in_dir:', in_dir))
    log_print(paste(Sys.time(), 'out_dir:', out_dir))
    log_print(paste(Sys.time(), 'zscore:', zscore))
} else {
    log_print(paste(Sys.time(), 'troubleshooting:', troubleshooting))
}


# ----------------------------------------------------------------------
# Configuration Files


config_file = file.path(out_dir, 'config.csv')
if (file.exists(config_file)) {
    config <- read.csv(config_file, header=TRUE, sep=',', check.names=FALSE)
} else {
    stop('config.csv is required to proceed')
}

metadata_file = file.path(in_dir, 'run_metadata.csv')
if (file.exists(metadata_file)) {
    run_metadata <- read.csv(metadata_file, header=TRUE, sep=',', check.names=FALSE)
    get_total_reads_from_metadata = TRUE
} else {
    get_total_reads_from_metadata = FALSE
}


# ----------------------------------------------------------------------
# Calculate Escape Genes

log_print(paste(Sys.time(), "Number of files found:", nrow(config)))

mouse_ids = config[['mouse_id']]
for (mouse_id in mouse_ids) {

    loop_start_time =  Sys.time()
    log_print(paste(loop_start_time, 'Loop started!'))

    row = (config['mouse_id']==mouse_id)
    if (sum(row) > 1) {
        stop("Mouse IDs must be unique!")
    }

    mat_reads_filename = config[row, "mat_reads_filenames"]
    pat_reads_filename = config[row, "pat_reads_filenames"]
    rpkms_filename = config[row, "rpkms_filenames"]

    mouse_gender = config[row, "mouse_gender"]
    mat_mouse_strain = config[row, "mat_mouse_strain"]
    pat_mouse_strain = config[row, "pat_mouse_strain"]


    # ----------------------------------------------------------------------
    # Read in raw data

    # get reads and merge exon lengths
    mat_reads = get_single_reads(
        mat_reads_filename,
        chromosomal_parentage='mat',
        mouse_strain=mat_mouse_strain
    )
    
    pat_reads = get_single_reads(
        pat_reads_filename,
        chromosomal_parentage='pat',
        mouse_strain=mat_mouse_strain
    )
    
    rpkms_file = file.path(in_dir, "rpkms", rpkms_filename)
    if (file.exists(rpkms_file) && !dir.exists(rpkms_file)) {
        rpkms = read_csv_or_tsv(rpkms_file)
    } else {
        # this will cause the rpkms to be estimated downstream
        rpkms = NULL
    }


    # ----------------------------------------------------------------------
    # Preprocessing

    # Calculate bias based on autosomal reads only. Should be pretty close to 1 if this is done right
    # In Zack's pipeline, the X chromosome is included in this calculation, whereas
    # Berletch specifies that this needs to be computed on the autosomal reads only
    num_snp_reads = sum(mat_reads[, 'count']) + sum(pat_reads[, 'count'])
    num_autosomal_reads_mat = sum(
        mat_reads[(mat_reads['chromosome']!='X') & (mat_reads['chromosome']!='Y'), 'count'])
    num_autosomal_reads_pat = sum(
        pat_reads[(pat_reads['chromosome']!='X') & (pat_reads['chromosome']!='Y'), 'count'])
    bias_xi_div_xa <- num_autosomal_reads_pat/num_autosomal_reads_mat

    log_print(paste(Sys.time(), 'SNP Reads Count:', num_snp_reads))
    log_print(paste(Sys.time(), 'Bias:', bias_xi_div_xa))

    # merge on gene name
    # only keep rows where gene names exist
    all_reads = merge(
        mat_reads[(mat_reads['gene_name']!=''), ],
        pat_reads[(pat_reads['gene_name']!=''), ],
        by='gene_name', suffixes=c('_mat', '_pat'),
        all.x=FALSE, all.y=FALSE,  # do not include null values
        na_matches = 'never'
    )  # we don't care about genes that don't have gene names, apparently
    
    # rename columns
    colnames(all_reads) <- lapply(
      colnames(all_reads),
      function(x) {
        step1 <- gsub('-', '_', x)  # convert mixed_case-col_names to fully snake_case
        step2 <- gsub('^chr_', 'chromosome_', step1)
        step3 <- gsub('^count_', 'num_reads_', step2)
        step4 <- gsub('^reads_', 'num_reads_', step3)
        return (step4)}
    )

    all_reads['mouse_id'] <- mouse_id
    all_reads['bias_xi_div_xa'] <- bias_xi_div_xa

    # reindex columns
    gene_id_cols = c(paste(gene_id_col, '_mat',  sep=''), paste(gene_id_col, '_pat',  sep=''))
    index_cols = c(
        'mouse_id', 'gene_name',
        gene_id_cols[0], 'chromosome_mat',
        gene_id_cols[1], 'chromosome_pat'
    )
    value_cols = items_in_a_not_b(colnames(all_reads), index_cols)
    all_reads <- all_reads[, c(index_cols, value_cols)]


    # ----------------------------------------------------------------------
    # Compute metrics

    # Duplicates
    # This is disabled, because the potential for this to break the pipeline is too high.
    # all_reads <- all_reads[!duplicated(all_reads[, c('gene_id_mat', 'chromosome_mat')]),]

    # filter Y chromosome, all the reads here should be 0 anyway
    all_reads <- all_reads[all_reads['chromosome_mat']!='Y',]

    # Estimate num_total_reads
    if (get_total_reads_from_metadata) {
        num_total_reads = run_metadata[(run_metadata['mouse_id']==mouse_id), 'num_total_reads']
    } else {
        num_total_reads = num_snp_reads / estimated_pct_snp_reads
    }
    all_reads['num_total_reads'] = num_total_reads
    all_reads['ratio_xi_over_xa'] = all_reads['num_reads_pat'] / all_reads['num_reads_mat']

    # Compute SRPM (allele-specific SNP-containing exonic reads per 10 million uniquely mapped reads)
    # This requires a-priori knowledge of how many reads
    all_reads[c('srpm_mat', 'srpm_pat')] = all_reads[c('num_reads_mat', 'num_reads_pat')] / num_total_reads * 1e7


    # RPKM (reads per kilobase of exon per million reads mapped)
    # Either merge RPKMs (good) or estimate them from the data (bad)
    if (!is.null(rpkms)) {

        # The RPKM is computed a priori and stored in anb RPKMs file
        all_reads = merge(
            all_reads,
            rpkms[c('gene_name', 'rpkm')],
            by='gene_name', suffixes=c('', '_'),
            all.x=TRUE, all.y=FALSE
            # na_matches = 'never'  # do not include null values
        )

    } else {

        # Estimate the RPKMs: num_total_reads = num_snp_reads / estimated_pct_snp_reads

        all_reads['rpm_mat'] = all_reads['num_reads_mat'] / num_total_reads * 1e6
        all_reads['rpm_pat'] = all_reads['num_reads_pat'] / num_total_reads * 1e6

        exon_lengths = coalesce1(all_reads["exon_length_mat"], all_reads["exon_length_pat"])
        all_reads['rpkm_mat'] = all_reads['rpm_mat'] / exon_lengths * 1000
        all_reads['rpkm_pat'] = all_reads['rpm_pat'] / exon_lengths * 1000

        # We're trying to estimate total reads, so rowSums is more appropriate here than rowMeans
        all_reads['rpkm'] = rowSums(all_reads[c('rpkm_mat', 'rpkm_pat')])

    }


    # ----------------------------------------------------------------------
    # Compute confidence intervals from binomial model

    # log_print('Computing confidence intervals...')

    # subset for downstream analysis
    x_reads <- all_reads[all_reads['chromosome_mat']=='X',]

    x_reads['total_reads'] = x_reads['num_reads_mat'] + x_reads['num_reads_pat']  # mat=Xa=n_i1, pat=Xi=n_i0
    x_reads['pct_xi'] = x_reads['num_reads_pat'] / x_reads['total_reads']
    x_reads[is.na(x_reads[, 'pct_xi']), 'pct_xi'] <- 0  # fillna with 0

    # human readable
    total_reads <- x_reads['total_reads']  # n_i
    pct_xi <- x_reads['pct_xi']

    x_reads['corrected_pct_xi'] <- pct_xi/(pct_xi+bias_xi_div_xa*(1-pct_xi))  # correction reduces to pct_xi if bias=1, eg. mouse_4
    corrected_pct_xi <- x_reads['corrected_pct_xi']

    x_reads['lower_confidence_interval'] <- corrected_pct_xi - zscore * (
        sqrt(corrected_pct_xi*(1-corrected_pct_xi)/total_reads))
    x_reads[is.na(x_reads[, 'lower_confidence_interval']), 'lower_confidence_interval'] <- 0  # fillna with 0

    x_reads['upper_confidence_interval'] <- corrected_pct_xi + zscore * (
        sqrt(corrected_pct_xi *(1-corrected_pct_xi )/total_reads))
    x_reads[is.na(x_reads[, 'upper_confidence_interval']), 'upper_confidence_interval'] <- 0  # fillna with 0


    # ----------------------------------------------------------------------
    # Filter flags

    x_reads['xi_srpm_gte_2'] <- as.integer(x_reads['srpm_pat'] >= 2)
    x_reads[is.na(x_reads['xi_srpm_gte_2']), 'xi_srpm_gte_2'] <- 0

    x_reads['rpkm_gt_1'] <- as.integer(x_reads['rpkm'] > 1)
    x_reads[is.na(x_reads['rpkm_gt_1']), 'rpkm_gt_1'] <- 0

    x_reads['lower_ci_gt_0'] <- as.integer(x_reads['lower_confidence_interval'] > 0)


    # ----------------------------------------------------------------------
    # Save

    # save
    if (!troubleshooting) {

        log_print(paste(Sys.time(), "Writing data..."))

        # ----------------------------------------------------------------------
        # X reads

        if (!dir.exists(file.path(out_dir, 'x_reads'))) {
            dir.create(file.path(out_dir, 'x_reads'), recursive=TRUE)
        }

        write.table(
            x_reads[intersect(x_read_output_cols, colnames(x_reads))],
            file = file.path(out_dir, 'x_reads', paste('x_reads-', mouse_id, '.csv', sep='')),
            row.names = FALSE,
            sep = ','
        )

        # ----------------------------------------------------------------------
        # Escape genes

        if (!dir.exists(file.path(out_dir, 'escape_genes'))) {
            dir.create(file.path(out_dir, 'escape_genes'), recursive=TRUE)
        }

        escape_genes = x_reads[
            (x_reads['rpkm_gt_1'] != 0) &
            (x_reads['xi_srpm_gte_2'] == 1) &
            (x_reads['lower_ci_gt_0'] == 1),
        ]

        write.table(
            escape_genes[order(escape_genes[, 'rpkm'], decreasing=TRUE),
                         intersect(escape_gene_cols, colnames(x_reads))],
            file = file.path(out_dir, 'escape_genes', paste('escape_genes-', mouse_id, '.csv', sep='')),
            row.names = FALSE,
            sep = ','
        )
    }

    loop_end_time = Sys.time()
    log_print(paste(Sys.time(), "Loop completed in:", difftime(loop_end_time, loop_start_time), '!'))
}

end_time = Sys.time()
log_print(paste('Script ended at:', Sys.time()))
log_print(paste("Script completed in:", difftime(end_time, start_time)))
log_close()
