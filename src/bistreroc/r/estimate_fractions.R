suppressWarnings(suppressMessages(library(data.table)))
suppressWarnings(suppressMessages(library(e1071)))
suppressWarnings(suppressMessages(library(parallel)))

command_arguments <- commandArgs(trailingOnly = TRUE)
signature_path <- command_arguments[1]
mixture_path <- command_arguments[2]
output_path <- command_arguments[3]
requested_core_count <- as.integer(command_arguments[4])

stopifnot(is.finite(requested_core_count), requested_core_count > 0)

signature_table <- fread(
    cmd = paste("cat", shQuote(signature_path)),
    check.names = FALSE
)
mixture_table <- fread(
    cmd = paste("cat", shQuote(mixture_path)),
    check.names = FALSE
)
signature_gene_column_name <- names(signature_table)[1]
mixture_gene_column_name <- names(mixture_table)[1]

stopifnot(
    ncol(signature_table) > 1,
    ncol(mixture_table) > 1,
    !anyDuplicated(signature_table[[signature_gene_column_name]]),
    !anyDuplicated(mixture_table[[mixture_gene_column_name]])
)

mixture_expression_column_names <- names(mixture_table)[-1]
if (max(as.matrix(mixture_table[, -1, with = FALSE])) < 50) {
    mixture_table[
        ,
        (mixture_expression_column_names) := lapply(
            .SD,
            function(expression_values) 2 ^ expression_values
        ),
        .SDcols = mixture_expression_column_names
    ]
}

shared_gene_names <- intersect(
    signature_table[[signature_gene_column_name]],
    mixture_table[[mixture_gene_column_name]]
)
stopifnot(length(shared_gene_names) > 1)
signature_table <- signature_table[match(
    shared_gene_names,
    signature_table[[signature_gene_column_name]]
)]
mixture_table <- mixture_table[match(
    shared_gene_names,
    mixture_table[[mixture_gene_column_name]]
)]
stopifnot(identical(
    signature_table[[signature_gene_column_name]],
    mixture_table[[mixture_gene_column_name]]
))

cell_type_names <- names(signature_table)[-1]
mixture_sample_names <- names(mixture_table)[-1]
signature_matrix <- data.matrix(signature_table[, -1, with = FALSE])
mixture_matrix <- data.matrix(mixture_table[, -1, with = FALSE])
stopifnot(
    all(is.finite(signature_matrix)),
    all(is.finite(mixture_matrix)),
    sd(signature_matrix) > 0,
    all(apply(mixture_matrix, 2, sd) > 0)
)

row_major_signature_values <- c(t(signature_matrix))
signature_mean <- Reduce(`+`, row_major_signature_values) /
    length(row_major_signature_values)
signature_sd <- sqrt(
    Reduce(`+`, (row_major_signature_values - signature_mean) ^ 2) /
        (length(row_major_signature_values) - 1)
)
standardized_signature_matrix <- (
    signature_matrix - signature_mean
) / signature_sd
standardized_mixture_matrix <- mixture_matrix
for (mixture_index in seq_len(ncol(mixture_matrix))) {
    mixture_values <- mixture_matrix[, mixture_index]
    mixture_mean <- Reduce(`+`, mixture_values) / length(mixture_values)
    mixture_sd <- sqrt(
        Reduce(`+`, (mixture_values - mixture_mean) ^ 2) /
            (length(mixture_values) - 1)
    )
    standardized_mixture_matrix[, mixture_index] <- (
        mixture_values - mixture_mean
    ) / mixture_sd
}
rounded_signature_matrix <- round(standardized_signature_matrix, 10)
rounded_mixture_matrix <- round(standardized_mixture_matrix, 10)
nu_candidates <- c(0.25, 0.50, 0.75)
model_fit_count <- ncol(standardized_mixture_matrix) * length(nu_candidates)
parallel_worker_count <- min(requested_core_count, model_fit_count)

fit_svr_model <- function(model_fit_index) {
    # Fit one linear nu-SVR configuration and return its scored fractions.

    mixture_index <- ((model_fit_index - 1) %/% length(nu_candidates)) + 1
    nu_index <- ((model_fit_index - 1) %% length(nu_candidates)) + 1
    svr_model <- svm(
        rounded_signature_matrix,
        rounded_mixture_matrix[, mixture_index],
        type = "nu-regression",
        kernel = "linear",
        nu = nu_candidates[nu_index],
        cost = 1,
        scale = FALSE
    )
    fraction_weights <- drop(crossprod(svr_model$coefs, svr_model$SV))
    fraction_weights[fraction_weights < 0] <- 0
    stopifnot(sum(fraction_weights) > 0)
    fractions <- fraction_weights / sum(fraction_weights)
    estimated_mixture <- drop(standardized_signature_matrix %*% fractions)
    c(
        fractions,
        p_value = 9999,
        correlation = cor(
            estimated_mixture,
            standardized_mixture_matrix[, mixture_index]
        ),
        rmse = sqrt(mean(
            (
                estimated_mixture -
                    standardized_mixture_matrix[, mixture_index]
            ) ^ 2
        )),
        selected_nu = nu_candidates[nu_index]
    )
}

if (parallel_worker_count == 1) {
    model_fit_rows <- lapply(seq_len(model_fit_count), fit_svr_model)
} else {
    model_fit_rows <- mclapply(
        seq_len(model_fit_count),
        fit_svr_model,
        mc.cores = parallel_worker_count,
        mc.set.seed = FALSE
    )
}
stopifnot(!any(vapply(model_fit_rows, inherits, logical(1), "try-error")))
model_fit_matrix <- do.call(rbind, model_fit_rows)

selected_model_rows <- lapply(
    seq_along(mixture_sample_names),
    function(mixture_index) {
        # Choose the minimum-RMSE model configuration for one mixture.

        first_model_fit_index <- (
            (mixture_index - 1) * length(nu_candidates) + 1
        )
        model_fit_indices <- first_model_fit_index:(
            first_model_fit_index + length(nu_candidates) - 1
        )
        mixture_model_fits <- model_fit_matrix[
            model_fit_indices,
            ,
            drop = FALSE
        ]
        mixture_model_fits[
            which.min(mixture_model_fits[, "rmse"]),
            ,
            drop = FALSE
        ]
    }
)
selected_result_matrix <- do.call(rbind, selected_model_rows)

output_connection <- file(output_path, open = "wt")
writeLines(
    paste(
        c("Mixture", cell_type_names, "P-value", "Correlation", "RMSE"),
        collapse = "\t"
    ),
    output_connection
)
for (mixture_index in seq_along(mixture_sample_names)) {
    reported_values <- selected_result_matrix[
        mixture_index,
        seq_len(ncol(selected_result_matrix) - 1)
    ]
    writeLines(
        paste(
            c(
                mixture_sample_names[mixture_index],
                sprintf("%.17f", reported_values)
            ),
            collapse = "\t"
        ),
        output_connection
    )
}
close(output_connection)
