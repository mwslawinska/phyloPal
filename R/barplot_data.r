#' Process barplot data at specified taxonomic level
#'
#' Prepares long-format data for barplot creation at a specified taxonomic level.
#' Handles unknown taxonomy, groups low-abundance taxa, and aggregates relative abundances (RA).
#' Works with data that already has RA (RA in the 0-1 range) calculated at ASV level.
#'
#' @param data A data frame with RA (relative abundance) already calculated per-ASV.
#' @param tax_level Character string specifying the taxonomic level to plot (e.g., "Family", "Genus").
#' @param group_vars Character vector of grouping variables (e.g. SampleID, SampleType).
#' @param SampleID_col Character string specifying the sample ID column (default: "SampleID")
#' @param low_abundance_threshold Numeric threshold below which taxa are considered low abundance (default: 0.01).
#' @param low_abundance_basis When to identify low-abundance taxa:
#'   \describe{
#'     \item{"per_sample"}{Mark taxa as low abundant at individual sample level BEFORE aggregation}
#'     \item{"post_aggregation"}{Mark taxa as low abundant AFTER aggregating across samples (default)}
#'   }
#' @param keep_ratype Controls whether low-abundance taxa are collapsed or kept separate:
#'   \describe{
#'     \item{"collapse"}{Relabel taxa below threshold as "low abundant" and collapse into one bin (default)}
#'     \item{"separate"}{Create \code{<tax_level>_original}; keep low/unknown as plot labels but do not merge originals}
#'   }
#' @param tax_original_suffix Suffix used to create the original-taxon column when \code{keep_ratype="separate"}
#'   (default: "_original"). Example: tax_level="Class" -> "Class_original".
#' @param preserve_higher_taxonomy Logical, whether to keep higher taxonomic levels (default: FALSE).
#' @param clean_taxonomy Logical, whether to clean taxonomy using replace_incertae_sedis_NAs (default: TRUE).
#' @param hierarchy Taxonomic hierarchy for cleaning (default: c("Domain", "Phylum", "Class", "Order", "Family", "Genus")).
#' @param clean_levels Levels to clean (default: c("Phylum", "Class", "Order", "Family", "Genus")). Uses \code{replace_incertae_sedis_NAs}
#' @param agg_fun Aggregation function for RA when summarizing groups: \code{"sum"} or \code{"mean"} (default: "sum").
#' @param normalize_by Character vector (or single string) giving grouping columns within which RA should sum to 1.
#'   If NULL: normalize within \code{SampleID_col} if present, else within \code{group_vars} excluding taxonomy columns.
#' @param drop_zero Logical; if TRUE, drop rows with \code{RA <= 0} at the end (default: TRUE). Useful to avoid inflating
#'   plotting keys with structural zeros.
#' @param unknown_label Label used for unknown taxonomy (default: "unknown").
#' @param low_label Label used for low-abundance bin (default: "low abundant").
#'
#' @section Low abundance handling (controlled by \code{keep_ratype}):
#' \describe{
#'   \item{collapse}{Relabel taxa below the threshold as "low abundant" and collapse into a single bin.}
#'   \item{separate}{Keep each original taxon in \code{<tax_level>_original}, while \code{<tax_level>}
#'   is replaced by \code{low_label}/\code{unknown_label} for flagged taxa.}
#' }
#'
#' @section Higher taxonomy handling (when \code{preserve_higher_taxonomy = TRUE}):
#' \describe{
#'   \item{collapse}{Higher taxonomy columns are also relabeled to "low abundant"/"unknown".}
#'   \item{separate}{True higher taxonomy is re-attached via \code{<tax_level>_original} into
#'   \code{<Higher>_true} columns (e.g. \code{Phylum_true}).}
#' }
#'
#' @return A data frame with processed relative abundances at specified taxonomic level.
#' \describe{
#'   \item{RA}{Relative abundance after aggregation and normalization (0–1).}
#'   \item{\code{<tax_level>}}{Taxon labels used for plotting (may include \code{low_label}/\code{unknown_label}).}
#'   \item{\code{<tax_level>_original}}{Only when \code{keep_ratype="separate"}: original taxon identity.}
#'   \item{\code{<Higher>_true}}{Only when \code{preserve_higher_taxonomy=TRUE} and \code{keep_ratype="separate"}: true higher ranks.}
#' }
#'
#' @examples
#' library(dplyr)
#'
#' # Minimal toy dataset (already aggregated at ASV level with RA in 0–1 range)
#' toy <- tibble::tibble(
#'   SampleID = c("S1","S1","S2","S2"),
#'   Group    = c("A","A","A","A"),
#'   Phylum   = c("P1","P2","P1","P2"),
#'   Class    = c("C1","C2","C1","C2"),
#'   RA       = c(0.9, 0.1, 0.6, 0.4)
#' )
#'
#' # Collapse low-abundance taxa into a single bin
#' out_collapse <- process_barplot_data(
#'   data = toy,
#'   tax_level = "Class",
#'   group_vars = "SampleID",
#'   low_abundance_threshold = 0.2,
#'   keep_ratype = "collapse",
#'   clean_taxonomy = FALSE
#' )
#'
#' # Keep low-abundance taxa separate but flagged
#' out_separate <- process_barplot_data(
#'   data = toy,
#'   tax_level = "Class",
#'   group_vars = "SampleID",
#'   low_abundance_threshold = 0.2,
#'   keep_ratype = "separate",
#'   clean_taxonomy = FALSE
#' )
#'
#' # In separate mode, an additional column "Class_original" appears
#' names(out_separate)
#' 
#' @export
process_barplot_data <- function(data,
                                 tax_level,
                                 group_vars,
                                 SampleID_col = "SampleID",
                                 low_abundance_threshold = 0.01,
                                 low_abundance_basis = c("post_aggregation", "per_sample"),
                                 keep_ratype = c("collapse", "separate"),
                                 tax_original_suffix = "_original",
                                 preserve_higher_taxonomy = FALSE,
                                 clean_taxonomy = TRUE,
                                 hierarchy = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
                                 clean_levels = c("Phylum", "Class", "Order", "Family", "Genus"),
                                 agg_fun = c("sum", "mean"),
                                 normalize_by = NULL,
                                 drop_zero = TRUE,
                                 unknown_label = "unknown",
                                 low_label = "low abundant") {
    # -----------------------------
    # Helpers / args
    # -----------------------------
    low_abundance_basis <- match.arg(low_abundance_basis)
    keep_ratype <- match.arg(keep_ratype)
    agg_fun <- match.arg(agg_fun)

    agg_apply <- function(x) {
        if (agg_fun == "sum") sum(x, na.rm = TRUE) else mean(x, na.rm = TRUE)
    }

    label_low_abundant_two_cols <- function(df, tax_level, threshold, tax_original_suffix,
                                            unknown_label, low_label) {
        orig_col <- paste0(tax_level, tax_original_suffix)

        if (!orig_col %in% names(df)) {
            df[[orig_col]] <- as.character(df[[tax_level]])
        } else {
            df[[orig_col]] <- as.character(df[[orig_col]])
        }

        is_unknown <- is.na(df[[orig_col]]) | df[[orig_col]] == unknown_label
        is_low <- (!is_unknown) & (df$RA < threshold)

        df[[tax_level]] <- dplyr::case_when(
            is_unknown ~ unknown_label,
            is_low ~ low_label,
            TRUE ~ df[[orig_col]]
        )
        df
    }

    propagate_higher <- function(df, tax_cols_to_keep, tax_level, unknown_label, low_label) {
        higher_cols <- setdiff(tax_cols_to_keep, tax_level)
        if (length(higher_cols) == 0) {
            return(df)
        }

        for (hc in higher_cols) {
            df[[hc]] <- dplyr::case_when(
                df[[tax_level]] == low_label ~ low_label,
                df[[tax_level]] == unknown_label ~ unknown_label,
                TRUE ~ df[[hc]]
            )
        }
        df
    }

    # -----------------------------
    # Validation
    # -----------------------------
    if (!"RA" %in% colnames(data)) stop("Data must contain 'RA' column (relative abundance).")
    if (!tax_level %in% colnames(data)) stop("tax_level '", tax_level, "' not found in data.")
    missing_group_vars <- setdiff(group_vars, colnames(data))
    if (length(missing_group_vars) > 0) {
        stop("Grouping variables not found in data: ", paste(missing_group_vars, collapse = ", "))
    }
    if (!is.numeric(low_abundance_threshold) || low_abundance_threshold < 0 || low_abundance_threshold > 1) {
        stop("low_abundance_threshold must be between 0 and 1.")
    }

    is_tbl <- inherits(data, "tbl_df")
    orig_col <- paste0(tax_level, tax_original_suffix)

    # -----------------------------
    # Step 1: Clean taxonomy
    # -----------------------------
    if (clean_taxonomy) {
        data_processed <- replace_incertae_sedis_NAs(
            df = data,
            hierarchy = hierarchy,
            clean_levels = clean_levels
        )
    } else {
        data_processed <- data
    }

    if (!SampleID_col %in% colnames(data_processed)) {
        stop("SampleID_col '", SampleID_col, "' not found in data.")
    }

    hierarchy_in_data <- intersect(hierarchy, colnames(data_processed))
    if (!tax_level %in% hierarchy_in_data) {
        stop(
            "tax_level must be one of taxonomy columns present in data. Found: ",
            paste(hierarchy_in_data, collapse = ", ")
        )
    }

    # Determine taxonomy columns to keep/append
    if (preserve_higher_taxonomy) {
        tax_idx <- which(hierarchy == tax_level)[1]
        tax_cols_to_keep <- intersect(hierarchy[1:tax_idx], colnames(data_processed))
    } else {
        tax_cols_to_keep <- tax_level
    }

    higher_cols <- setdiff(tax_cols_to_keep, tax_level)

    # -----------------------------
    # Build TRUE taxonomy map keyed by original taxon (only if needed)
    # -----------------------------
    tax_map <- NULL
    if (preserve_higher_taxonomy && keep_ratype == "separate" && length(higher_cols) > 0) {
        tax_map <- data_processed %>%
            dplyr::mutate(.orig = as.character(.data[[tax_level]])) %>%
            dplyr::select(dplyr::all_of(c(higher_cols)), .orig) %>%
            dplyr::distinct() %>%
            dplyr::rename(!!orig_col := .orig) %>%
            dplyr::group_by(.data[[orig_col]]) %>%
            dplyr::summarise(
                dplyr::across(dplyr::all_of(higher_cols), dplyr::first),
                .groups = "drop"
            )
    }

    # -----------------------------
    # Step 3: normalize_by
    # -----------------------------
    if (is.null(normalize_by)) {
        if (SampleID_col %in% colnames(data_processed)) {
            normalize_by <- SampleID_col
        } else {
            normalize_by <- setdiff(group_vars, hierarchy_in_data)
            if (length(normalize_by) == 0) stop("Could not infer normalize_by. Provide normalize_by explicitly.")
        }
    } else {
        missing_norm <- setdiff(normalize_by, colnames(data_processed))
        if (length(missing_norm) > 0) stop("normalize_by not found in data: ", paste(missing_norm, collapse = ", "))
    }

    # -----------------------------
    # Step 4: per-sample basis (optional)
    # -----------------------------
    sample_id_col <- SampleID_col

    if (low_abundance_basis == "per_sample" && low_abundance_threshold > 0) {
        if (is.null(sample_id_col)) stop("low_abundance_basis='per_sample' requires a SampleID column in the data.")

        grouping_sample_tax <- unique(c(group_vars, sample_id_col, tax_cols_to_keep))

        data_processed <- data_processed %>%
            dplyr::ungroup() %>%
            dplyr::group_by(dplyr::across(dplyr::all_of(grouping_sample_tax))) %>%
            dplyr::summarise(RA = agg_apply(.data$RA), .groups = "drop") %>%
            dplyr::group_by(.data[[sample_id_col]]) %>%
            dplyr::mutate(
                .denom = sum(.data$RA, na.rm = TRUE),
                RA = ifelse(.denom > 0, .data$RA / .denom, .data$RA),
                .denom = NULL
            ) %>%
            dplyr::ungroup()

        data_processed <- label_low_abundant_two_cols(
            df = data_processed,
            tax_level = tax_level,
            threshold = low_abundance_threshold,
            tax_original_suffix = tax_original_suffix,
            unknown_label = unknown_label,
            low_label = low_label
        )

        # Only propagate higher labels in collapse mode
        if (preserve_higher_taxonomy && keep_ratype == "collapse") {
            data_processed <- propagate_higher(
                data_processed, tax_cols_to_keep, tax_level,
                unknown_label = unknown_label, low_label = low_label
            )
        }

        if (keep_ratype == "collapse") {
            data_processed <- data_processed %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(grouping_sample_tax))) %>%
                dplyr::summarise(RA = sum(.data$RA, na.rm = TRUE), .groups = "drop")
        }
    }

    # -----------------------------
    # Step 5: aggregate across group_vars
    # -----------------------------
    grouping_cols <- unique(c(group_vars, tax_cols_to_keep))

    if (keep_ratype == "separate") {
        if (!(orig_col %in% colnames(data_processed))) {
            data_processed[[orig_col]] <- as.character(data_processed[[tax_level]])
        }
        grouping_cols <- unique(c(grouping_cols, orig_col))
    }

    data_aggregated <- data_processed %>%
        dplyr::ungroup() %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
        dplyr::summarise(RA = agg_apply(.data$RA), .groups = "drop")

    # -----------------------------
    # Step 6: normalize
    # -----------------------------
    data_aggregated <- data_aggregated %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(normalize_by))) %>%
        dplyr::mutate(
            .denom = sum(.data$RA, na.rm = TRUE),
            RA = ifelse(.denom > 0, .data$RA / .denom, .data$RA),
            .denom = NULL
        ) %>%
        dplyr::ungroup()

    # -----------------------------
    # Step 7: post-aggregation labeling (optional)
    # -----------------------------
    if (low_abundance_basis == "post_aggregation" && low_abundance_threshold > 0) {
        if (keep_ratype == "separate" && !(orig_col %in% colnames(data_aggregated))) {
            data_aggregated[[orig_col]] <- as.character(data_aggregated[[tax_level]])
        }

        data_aggregated <- label_low_abundant_two_cols(
            df = data_aggregated,
            tax_level = tax_level,
            threshold = low_abundance_threshold,
            tax_original_suffix = tax_original_suffix,
            unknown_label = unknown_label,
            low_label = low_label
        )

        if (preserve_higher_taxonomy && keep_ratype == "collapse") {
            data_aggregated <- propagate_higher(
                data_aggregated, tax_cols_to_keep, tax_level,
                unknown_label = unknown_label, low_label = low_label
            )
        }

        if (keep_ratype == "collapse") {
            grouping_cols_collapse <- unique(c(group_vars, tax_cols_to_keep))
            data_aggregated <- data_aggregated %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols_collapse))) %>%
                dplyr::summarise(RA = sum(.data$RA, na.rm = TRUE), .groups = "drop") %>%
                dplyr::group_by(dplyr::across(dplyr::all_of(normalize_by))) %>%
                dplyr::mutate(
                    .denom = sum(.data$RA, na.rm = TRUE),
                    RA = ifelse(.denom > 0, .data$RA / .denom, .data$RA),
                    .denom = NULL
                ) %>%
                dplyr::ungroup()
        }
    }

    # -----------------------------
    # Reattach TRUE higher taxonomy in separate mode + set plot-level higher taxonomy labels
    # -----------------------------
    if (preserve_higher_taxonomy && keep_ratype == "separate" && length(higher_cols) > 0) {
        # Join true taxonomy as *_true
        tax_map_true <- tax_map %>%
            dplyr::rename_with(~ paste0(.x, "_true"), dplyr::all_of(higher_cols))

        data_aggregated <- data_aggregated %>%
            dplyr::select(-dplyr::any_of(higher_cols)) %>%
            dplyr::left_join(tax_map_true, by = orig_col)

        # Create plot-level higher taxonomy columns (Domain/Phylum/...) from true ones
        for (hc in higher_cols) {
            data_aggregated[[hc]] <- data_aggregated[[paste0(hc, "_true")]]
        }

        # Force plot-level higher taxonomy to low/unknown when plotted taxon is low/unknown
        is_low <- data_aggregated[[tax_level]] == low_label
        is_unk <- data_aggregated[[tax_level]] == unknown_label

        for (hc in higher_cols) {
            data_aggregated[[hc]] <- dplyr::case_when(
                is_low ~ low_label,
                is_unk ~ unknown_label,
                TRUE ~ data_aggregated[[hc]]
            )
        }
    }

    # Drop structural zeros at the end (recommended for barplots)
    if (isTRUE(drop_zero)) {
        data_aggregated <- data_aggregated %>%
            dplyr::filter(.data$RA > 0)
    }

    if (is_tbl) data_aggregated <- tibble::as_tibble(data_aggregated)
    data_aggregated
}
