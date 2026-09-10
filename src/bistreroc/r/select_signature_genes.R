suppressWarnings(suppressMessages(library(parallel)))

command_arguments <- commandArgs(trailingOnly = TRUE)
cell_type_profile_path <- command_arguments[1]
signature_output_path <- command_arguments[2]
minimum_markers_per_cell_type <- as.integer(command_arguments[3])
maximum_markers_per_cell_type <- as.integer(command_arguments[4])
q_value_threshold <- as.numeric(command_arguments[5])
configured_core_count <- as.integer(command_arguments[6])

stopifnot(
    minimum_markers_per_cell_type > 0,
    maximum_markers_per_cell_type >= minimum_markers_per_cell_type,
    q_value_threshold >= 0,
    q_value_threshold <= 1,
    configured_core_count > 0
)

reference_expression_table <- read.delim(
    cell_type_profile_path,
    header = TRUE,
    row.names = 1,
    check.names = FALSE,
    quote = "",
    comment.char = ""
)
reference_expression_matrix <- as.matrix(reference_expression_table)
storage.mode(reference_expression_matrix) <- "double"
stopifnot(ncol(reference_expression_matrix) %% 5 == 0)
cell_type_count <- ncol(reference_expression_matrix) / 5
output_cell_type_names <- colnames(reference_expression_matrix)[
    seq.int(1, ncol(reference_expression_matrix), by = 5)
]
sanitized_cell_type_names <- make.names(
    output_cell_type_names,
    unique = TRUE
)
cell_type_index_by_column <- rep(seq_len(cell_type_count), each = 5)

original_gene_names <- rownames(reference_expression_matrix)
sanitized_gene_names <- make.names(original_gene_names, unique = TRUE)
retained_gene_mask <- !(
    grepl("-", sanitized_gene_names, fixed = TRUE) &
        grepl("IG", sanitized_gene_names, fixed = TRUE)
)
reference_expression_matrix <- reference_expression_matrix[
    retained_gene_mask,
    ,
    drop = FALSE
]
sanitized_gene_names <- sanitized_gene_names[retained_gene_mask]
rownames(reference_expression_matrix) <- sanitized_gene_names
sorted_gene_indices <- order(sanitized_gene_names)
reference_expression_matrix <- reference_expression_matrix[
    sorted_gene_indices,
    ,
    drop = FALSE
]
sanitized_gene_names <- sanitized_gene_names[sorted_gene_indices]

if (max(reference_expression_matrix) > 50) {
    logarithm_offset <- as.double(min(reference_expression_matrix) <= 0)
    reference_expression_matrix[reference_expression_matrix < 0] <- 0
    log2_reference_values <- log10(
        reference_expression_matrix + logarithm_offset
    ) / log10(2)
} else {
    log2_reference_values <- reference_expression_matrix
}
linearized_reference_values <- 2 ^ log2_reference_values

calculate_q_values <- function(p_values) {
    # Apply the Storey-like marker-selection recurrence.

    p_value_count <- length(p_values)
    sorted_p_value_indices <- order(
        p_values,
        seq_along(p_values),
        method = "radix"
    )
    sorted_p_values <- p_values[sorted_p_value_indices]
    estimated_null_hypothesis_count <- 2 * sum(p_values > 0.5)
    sorted_q_values <- numeric(p_value_count)
    sorted_q_values[[p_value_count]] <- (
        sorted_p_values[[p_value_count]] *
            estimated_null_hypothesis_count / p_value_count
    )
    if (p_value_count > 1) {
        for (sorted_index in seq.int(p_value_count - 1, 1)) {
            sorted_q_values[[sorted_index]] <- min(
                1,
                sorted_q_values[[sorted_index + 1]],
                sorted_p_values[[sorted_index]] *
                    estimated_null_hypothesis_count / (sorted_index + 1)
            )
        }
    }
    q_values <- numeric(p_value_count)
    q_values[sorted_p_value_indices] <- sorted_q_values
    q_values
}

safe_divide <- function(numerator, denominator) {
    # Mark positive ratios with near-zero denominators as undefined.

    division_result <- numerator / denominator
    guarded_value_mask <- numerator > 1e-10 & denominator < 1e-10
    division_result[guarded_value_mask] <- NaN
    division_result
}

cell_type_profile_matrix <- matrix(
    nrow = nrow(reference_expression_matrix),
    ncol = length(sanitized_cell_type_names),
    dimnames = list(sanitized_gene_names, output_cell_type_names)
)
ranked_marker_genes_by_cell_type <- vector(
    "list",
    length(sanitized_cell_type_names)
)

for (cell_type_index in seq_along(sanitized_cell_type_names)) {
    target_column_mask <- cell_type_index_by_column == cell_type_index
    target_log2_values <- log2_reference_values[
        ,
        target_column_mask,
        drop = FALSE
    ]
    background_log2_values <- log2_reference_values[
        ,
        !target_column_mask,
        drop = FALSE
    ]
    target_replicate_count <- ncol(target_log2_values)
    background_replicate_count <- ncol(background_log2_values)
    target_log2_mean <- rowMeans(target_log2_values)
    background_log2_mean <- rowMeans(background_log2_values)
    log_fold_change <- target_log2_mean - background_log2_mean
    target_log2_variance <- rowSums(
        (target_log2_values - target_log2_mean) ^ 2
    ) / (target_replicate_count - 1)
    background_log2_variance <- rowSums(
        (background_log2_values - background_log2_mean) ^ 2
    ) / (background_replicate_count - 1)
    squared_standard_error <- (
        target_log2_variance / target_replicate_count +
            background_log2_variance / background_replicate_count
    )
    t_statistic <- safe_divide(
        log_fold_change,
        sqrt(squared_standard_error)
    )
    degrees_of_freedom <- safe_divide(squared_standard_error ^ 2, (
        (target_log2_variance / target_replicate_count) ^ 2 /
            (target_replicate_count - 1) +
            (background_log2_variance / background_replicate_count) ^ 2 /
                (background_replicate_count - 1)
    ))
    p_values <- (
        pt(-abs(t_statistic), degrees_of_freedom) +
            pt(abs(t_statistic), degrees_of_freedom, lower.tail = FALSE)
    )
    p_values[is.nan(p_values)] <- 1
    q_values <- calculate_q_values(p_values)
    p_value_threshold <- if (target_replicate_count == 2) 0.1 else 0.05
    eligible_gene_indices <- which(
        p_values < p_value_threshold &
            log_fold_change >= 2 &
            q_values <= q_value_threshold
    )
    marker_rank_order <- order(
        -log_fold_change[eligible_gene_indices],
        eligible_gene_indices,
        method = "radix"
    )
    ranked_marker_genes_by_cell_type[[cell_type_index]] <- sanitized_gene_names[
        eligible_gene_indices[marker_rank_order]
    ]
    cell_type_profile_matrix[, cell_type_index] <- apply(
        linearized_reference_values[, target_column_mask, drop = FALSE],
        1,
        median
    )
}

calculate_condition_number <- function(candidate_marker_count) {
    # Calculate exact matrix kappa for one per-type marker cutoff.

    candidate_marker_genes <- unique(unlist(lapply(
        ranked_marker_genes_by_cell_type,
        head,
        n = candidate_marker_count
    )))
    candidate_signature_row_mask <- sanitized_gene_names %in%
        candidate_marker_genes
    if (!any(candidate_signature_row_mask)) {
        return(Inf)
    }
    kappa(
        cell_type_profile_matrix[
            candidate_signature_row_mask,
            ,
            drop = FALSE
        ],
        exact = TRUE
    )
}

candidate_marker_counts <- minimum_markers_per_cell_type:
    maximum_markers_per_cell_type
active_worker_count <- min(configured_core_count, length(candidate_marker_counts))
if (active_worker_count == 1) {
    condition_numbers <- lapply(
        candidate_marker_counts,
        calculate_condition_number
    )
} else {
    condition_numbers <- mclapply(
        candidate_marker_counts,
        calculate_condition_number,
        mc.cores = active_worker_count,
        mc.set.seed = FALSE
    )
}
stopifnot(!any(vapply(condition_numbers, inherits, logical(1), "try-error")))
selected_marker_count <- candidate_marker_counts[
    which.min(unlist(condition_numbers))
]
selected_marker_genes <- unique(unlist(lapply(
    ranked_marker_genes_by_cell_type,
    head,
    n = selected_marker_count
)))
selected_signature_row_mask <- sanitized_gene_names %in% selected_marker_genes
signature_matrix <- cell_type_profile_matrix[
    selected_signature_row_mask,
    ,
    drop = FALSE
]

formatted_signature_values <- matrix(
    sprintf("%.17f", as.vector(t(signature_matrix))),
    nrow = nrow(signature_matrix),
    byrow = TRUE
)
output_text_lines <- c(
    paste(c("NAME", colnames(signature_matrix)), collapse = "\t"),
    apply(
        cbind(rownames(signature_matrix), formatted_signature_values),
        1,
        paste,
        collapse = "\t"
    )
)
writeLines(output_text_lines, signature_output_path, useBytes = TRUE)
message(
    "selected ",
    selected_marker_count,
    " markers per cell type (",
    nrow(signature_matrix),
    " unique genes)"
)
