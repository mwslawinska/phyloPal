#' Aggregate microbiome data for alluvial plot visualization
#'
#' Aggregates ASV-level data with pre-calculated relative abundance (\code{RA}) to a chosen
#' taxonomic level and computes mean group compositions (sample-weighted). Optionally
#' completes missing taxa with zeros to ensure a rectangular structure required by
#' alluvial/Sankey plots.
#'
#' @param data Data frame with \code{RA} already calculated at ASV level.
#' @param tax_level Character string specifying the taxonomic level to aggregate to (e.g. \code{"Family"}).
#' @param group_col Character string specifying the grouping column (e.g. \code{"SampleType"}).
#' @param SampleID_col Character string specifying the sample ID column (default: \code{"SampleID"}).
#' @param groups Character vector of groups to include (NULL = all groups present in \code{group_col}).
#' @param complete_zero Logical; if TRUE, completes missing taxa with \code{RA = 0} within each sample
#'   before summarizing (default: TRUE). 
#' @param clean_taxonomy Logical; if TRUE, clean taxonomy using \code{replace_incertae_sedis_NAs} (default: TRUE).
#' @param preserve_higher_taxonomy Logical; if TRUE, keep higher taxonomic ranks up to \code{tax_level}
#'   (default: FALSE).
#' @param hierarchy Taxonomic hierarchy for cleaning (default:
#'   \code{c("Domain","Phylum","Class","Order","Family","Genus")}).
#' @param clean_levels Taxonomic levels to clean (default: \code{c("Phylum","Class","Order","Family","Genus")}).
#'
#' @section \code{complete_zero} behavior:
#' \describe{
#'   \item{\code{TRUE}}{Completes missing taxa within each sample with \code{RA = 0}, so taxa are averaged over all samples.}
#'   \item{\code{FALSE}}{Taxa absent from some samples are averaged only over samples where they appear (can inflate sporadic taxa).}
#' }
#' 
#' @section \code{preserve_higher_taxonomy} behavior:
#' \describe{
#'   \item{\code{TRUE}}{Keeps higher ranks up to \code{tax_level}. Note: rows created by the final completion
#'   may have missing higher-rank labels unless those columns are also completed/filled.}
#'   \item{\code{FALSE}}{Only returns \code{tax_level} plus grouping columns.}
#' }
#' 
#' @section Data processing in the function step-by-step:
#' \enumerate{
#'   \item Optionally clean taxonomy using \code{replace_incertae_sedis_NAs}.
#'   \item Aggregate ASV-level \code{RA} to \code{tax_level} within each sample using \code{sum()}.
#'   \item If \code{complete_zero = TRUE}, add missing taxa per sample with \code{RA = 0}
#'         (taxa are taken from all values observed in \code{data[[tax_level]]}).
#'   \item Compute the mean \code{RA} per \code{group_col} (mean across samples), then normalize within each group so
#'         \code{sum(RA) = 1}.
#'   \item Ensure a complete \code{group_col} \eqn{\times} \code{tax_level} grid (missing combinations are set to \code{RA = 0}).
#' }
#'
#' @return A data frame with mean relative abundance per \code{group_col} and \code{tax_level}.
#' \describe{
#'   \item{RA}{Mean relative abundance per group after normalization (0–1).}
#'   \item{\code{<group_col>}}{Grouping variable used on the alluvial x-axis (the column named by \code{group_col}).}
#'   \item{\code{<tax_level>}}{Taxon identifier at the aggregation level (the column named by \code{tax_level}).}
#'   \item{Higher taxonomy columns}{Only when \code{preserve_higher_taxonomy = TRUE}: higher ranks up to \code{tax_level}.}
#' }
#'
#' @examples
#' library(dplyr)
#' library(tidyr)
#'
#' # Toy example showing how complete_zero changes group means:
#' # Taxon B is absent from sample S1 (no row), but present in S2.
#' toy <- tibble::tibble(
#'   SampleID = c("S1","S2","S2"),
#'   Group    = c("G1","G1","G1"),
#'   Class    = c("A","A","B"),
#'   RA       = c(1.0, 0.5, 0.5)
#' )
#'
#' # Without completion: B is averaged only over samples where it appears (inflated)
#' out_no0 <- prepare_alluvial_data(
#'   data = toy,
#'   tax_level = "Class",
#'   group_col = "Group",
#'   SampleID_col = "SampleID",
#'   complete_zero = FALSE,
#'   clean_taxonomy = FALSE
#' )
#'
#' # With completion: missing taxa count as 0 in those samples
#' out_0 <- prepare_alluvial_data(
#'   data = toy,
#'   tax_level = "Class",
#'   group_col = "Group",
#'   SampleID_col = "SampleID",
#'   complete_zero = TRUE,
#'   clean_taxonomy = FALSE
#' )
#'
#' out_no0 %>% arrange(Class)
#' out_0   %>% arrange(Class)
#' 
#' @export
prepare_alluvial_data <- function(data,
                                  tax_level,
                                  group_col = "Group",
                                  SampleID_col = "SampleID",
                                  groups = NULL,
                                  complete_zero = TRUE,
                                  clean_taxonomy = TRUE,
                                  preserve_higher_taxonomy = FALSE,
                                  hierarchy = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
                                  clean_levels = c("Phylum", "Class", "Order", "Family", "Genus")) {
  # Input validation
  required_cols <- c(tax_level, group_col, SampleID_col, "RA")
  missing_cols <- setdiff(required_cols, colnames(data))
  if (length(missing_cols) > 0) {
    stop("Missing required columns: ", paste(missing_cols, collapse = ", "))
  }

  # If groups not specified, use all unique groups
  if (is.null(groups)) {
    groups <- unique(data[[group_col]])
  }

  # Step 1: Clean taxonomy if requested
  if (clean_taxonomy) {
    data <- data %>%
      replace_incertae_sedis_NAs(
        df = .,
        hierarchy = hierarchy,
        clean_levels = clean_levels
      )
  }

  # Step 2: Determine taxonomic columns to keep
  possible_tax_levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_cols_in_data <- intersect(possible_tax_levels, colnames(data))

  if (preserve_higher_taxonomy) {
    tax_level_index <- which(possible_tax_levels == tax_level)
    if (length(tax_level_index) == 0) {
      stop("tax_level must be one of: ", paste(possible_tax_levels, collapse = ", "))
    }
    higher_and_current_tax <- possible_tax_levels[1:tax_level_index]
    tax_cols_to_keep <- intersect(higher_and_current_tax, tax_cols_in_data)
  } else {
    tax_cols_to_keep <- tax_level
  }

  # Convert to symbols for NSE
  tax_sym <- rlang::sym(tax_level)
  group_sym <- rlang::sym(group_col)
  sample_sym <- rlang::sym(SampleID_col)

  # Step 3: Aggregate to taxonomic level per sample
  if (complete_zero) {
    if (group_col == SampleID_col) {
      # When group and sample are the same
      grouping_cols <- unique(c(SampleID_col, tax_cols_to_keep))

      data_agg <- data %>%
        dplyr::ungroup() %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
        dplyr::summarize(RA = sum(RA), .groups = "drop") %>%
        dplyr::group_by(!!sample_sym) %>%
        tidyr::complete(!!tax_sym := unique(dplyr::pull(data, !!tax_sym)),
          fill = list(RA = 0)
        ) %>%
        dplyr::ungroup()

      # Fill higher taxonomy columns
      if (preserve_higher_taxonomy) {
        for (tax_col in setdiff(tax_cols_to_keep, tax_level)) {
          data_agg <- data_agg %>%
            tidyr::fill(!!rlang::sym(tax_col), .direction = "downup")
        }
      }
    } else {
      # When group and sample are different
      grouping_cols <- unique(c(SampleID_col, group_col, tax_cols_to_keep))

      data_agg <- data %>%
        dplyr::ungroup() %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
        dplyr::summarize(RA = sum(RA), .groups = "drop") %>%
        dplyr::group_by(!!sample_sym) %>%
        tidyr::complete(!!tax_sym := unique(dplyr::pull(data, !!tax_sym)),
          fill = list(RA = 0)
        ) %>%
        tidyr::fill(!!group_sym, .direction = "downup") %>%
        dplyr::ungroup()

      # Fill higher taxonomy columns
      if (preserve_higher_taxonomy) {
        for (tax_col in setdiff(tax_cols_to_keep, tax_level)) {
          data_agg <- data_agg %>%
            tidyr::fill(!!rlang::sym(tax_col), .direction = "downup")
        }
      }
    }
  } else {
    # No zero completion
    if (group_col == SampleID_col) {
      grouping_cols <- unique(c(SampleID_col, tax_cols_to_keep))

      data_agg <- data %>%
        dplyr::ungroup() %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
        dplyr::summarize(RA = sum(RA), .groups = "drop")
    } else {
      grouping_cols <- unique(c(SampleID_col, group_col, tax_cols_to_keep))

      data_agg <- data %>%
        dplyr::ungroup() %>%
        dplyr::group_by(dplyr::across(dplyr::all_of(grouping_cols))) %>%
        dplyr::summarize(RA = sum(RA), .groups = "drop")
    }
  }

  # Step 4: Calculate mean RA per group and normalize
  if (preserve_higher_taxonomy) {
    summary_grouping <- unique(c(group_col, tax_cols_to_keep))
  } else {
    summary_grouping <- c(group_col, tax_level)
  }

  data_summary <- data_agg %>%
    dplyr::group_by(dplyr::across(dplyr::all_of(summary_grouping))) %>%
    dplyr::summarize(RA = mean(RA), .groups = "drop") %>%
    dplyr::group_by(!!group_sym) %>%
    dplyr::mutate(RA = RA / sum(RA)) %>%
    dplyr::ungroup() %>%
    dplyr::filter(!!group_sym %in% groups)

  # Ensure complete alluvial structure
  all_taxa <- unique(data_summary[[tax_level]])
  all_groups <- unique(data_summary[[group_col]])

  data_summary <- data_summary %>%
    tidyr::complete(
      !!rlang::sym(group_col) := all_groups,
      !!rlang::sym(tax_level) := all_taxa,
      fill = list(RA = 0)
    )

  # Step 5: Set higher taxonomy for unknown
  if (preserve_higher_taxonomy) {
    for (tax_col in setdiff(tax_cols_to_keep, tax_level)) {
      data_summary <- data_summary %>%
        dplyr::mutate(
          !!tax_col := dplyr::case_when(
            .data[[tax_level]] == "unknown" ~ "unknown",
            TRUE ~ .data[[tax_col]]
          )
        )
    }
  }

  return(data_summary)
}

#' Classify taxa by abundance patterns across groups
#'
#' Classifies taxa as shared/unique and abundant/low abundant based on their
#' relative abundances across multiple groups. Optionally detects “mixed abundance”
#' taxa that are abundant in some groups but low abundant in others.
#'
#' @param data A data frame containing at least \code{tax_level}, \code{group_col}, and \code{RA}
#'   (typically the output of \code{prepare_alluvial_data()}).
#' @param tax_level Character string specifying the taxonomic level (e.g. \code{"Family"}, \code{"Class"}).
#' @param group_col Character string specifying the grouping column (default: \code{"Group"}).
#' @param groups Character vector of groups to compare. If NULL, uses all unique values in \code{data[[group_col]]}.
#' @param low_abundance_threshold Numeric threshold below which taxa are considered low abundant (default: 0.01).
#' @param special_taxa Character vector of taxa that should never be marked as low abundant
#'   (default: \code{character(0)}).
#' @param detect_mixed_abundance Logical; if TRUE, detect taxa that are abundant in some groups and low abundant in others
#'   and label as \code{"shared mixed abundance"} (default: TRUE). These taxa keep their own palette key
#'   (\code{tax_color = <taxon>}) rather than being mapped to a low-abundance bin.
#' @param unknown_label Label used for unknown taxa (default: \code{"unknown"}).
#'
#' @details
#' This function does not set plotting order (legend or stratum stacking). In the
#' palette-driven workflow, ordering is handled in \code{plot_alluvial()} via the
#' order of \code{names(custom_palette)}.
#'
#' @return The input data with additional classification columns.
#' \describe{
#'   \item{tax_type}{Primary abundance pattern classification (e.g. shared/unique; abundant/low; and optionally mixed).}
#'   \item{category}{Final category used for summaries and plotting (may match \code{tax_type} depending on workflow).}
#'   \item{tax_color}{A label used for color mapping (e.g. low-abundant bins may share a common color label).}
#' }
#' \item{tax_val}{Character copy of \code{<tax_level>} used internally for classification.}
#' Other columns (including higher taxonomy ranks) are preserved if present in \code{data}.
#'
#' @examples
#' library(dplyr)
#'
#' # Two groups, three taxa:
#' # - A is abundant in both groups (shared abundant)
#' # - B is abundant in G1 but low in G2 (mixed abundance if enabled)
#' # - C appears only in G2 (unique)
#' toy <- tibble::tibble(
#'   Group = c("G1", "G1", "G1", "G2", "G2", "G2"),
#'   Class = c("A", "B", "unknown", "A", "B", "C"),
#'   RA    = c(0.80, 0.15, 0.05, 0.80, 0.005, 0.195)
#' )
#'
#' out <- classify_taxa_patterns(
#'   data = toy,
#'   tax_level = "Class",
#'   group_col = "Group",
#'   low_abundance_threshold = 0.01,
#'   detect_mixed_abundance = TRUE
#' )
#'
#' out %>% arrange(Group, Class)
#'
#' out_mixed <- classify_taxa_patterns(toy, "Class", "Group", detect_mixed_abundance = TRUE)
#' out_nomix <- classify_taxa_patterns(toy, "Class", "Group", detect_mixed_abundance = FALSE)
#' out_mixed %>%
#'   filter(Class == "B") %>%
#'   dplyr::select(Group, Class, RA, tax_type, tax_color)
#' out_nomix %>%
#'   filter(Class == "B") %>%
#'   dplyr::select(Group, Class, RA, tax_type, tax_color)
#'
#' @export
classify_taxa_patterns <- function(data,
                                   tax_level,
                                   group_col = "Group",
                                   groups = NULL,
                                   low_abundance_threshold = 0.01,
                                   special_taxa = character(0),
                                   detect_mixed_abundance = TRUE,
                                   unknown_label = "unknown") {
  # Input validation
  if (!all(c(tax_level, group_col, "RA") %in% colnames(data))) {
    stop("Data must contain columns: ", tax_level, ", ", group_col, ", and RA")
  }
  if (low_abundance_threshold < 0 || low_abundance_threshold > 1) {
    stop("low_abundance_threshold must be between 0 and 1")
  }

  if (is.null(groups)) {
    groups <- unique(data[[group_col]])
  }
  if (length(groups) == 0) stop("No groups found/provided in 'groups' or 'group_col'.")

  tax_sym <- rlang::sym(tax_level)
  group_sym <- rlang::sym(group_col)

  # Create tax_val for easier handling
  data <- data %>%
    dplyr::mutate(tax_val = as.character(!!tax_sym))

  # Step 1: Classify each taxon in each group
  low_abund_list <- vector("list", length(groups))
  high_abund_list <- vector("list", length(groups))
  absent_list <- vector("list", length(groups))

  for (i in seq_along(groups)) {
    group_data <- data %>%
      dplyr::filter(!!group_sym == groups[i]) %>%
      dplyr::mutate(
        abundance_status = dplyr::case_when(
          RA == 0 ~ "absent",
          RA < low_abundance_threshold &
            tax_val != unknown_label &
            !(tax_val %in% special_taxa) ~ "low",
          TRUE ~ "high"
        )
      )

    low_abund_list[[i]] <- group_data %>%
      dplyr::filter(abundance_status == "low") %>%
      dplyr::pull(tax_val) %>%
      unique()

    high_abund_list[[i]] <- group_data %>%
      dplyr::filter(abundance_status == "high") %>%
      dplyr::pull(tax_val) %>%
      unique()

    absent_list[[i]] <- group_data %>%
      dplyr::filter(abundance_status == "absent") %>%
      dplyr::pull(tax_val) %>%
      unique()
  }

  # Step 2: Identify patterns across groups
  shared_low <- setdiff(Reduce(intersect, low_abund_list), special_taxa)
present_set <- lapply(seq_along(groups), function(i) {
  c(low_abund_list[[i]], high_abund_list[[i]])
})
# Shared = abundant (high) in at least 2 groups
high_in_n_groups <- table(unlist(high_abund_list))
shared_taxa <- names(high_in_n_groups[high_in_n_groups >= 2])
unique_taxa  <- setdiff(unlist(present_set), shared_taxa)

  # Mixed abundance: high in one group, low in another
  if (isTRUE(detect_mixed_abundance) && length(groups) >= 2) {
    mixed_abundance <- unique(unlist(
      apply(combn(seq_along(groups), 2), 2, function(idx) {
        c(
          intersect(high_abund_list[[idx[1]]], low_abund_list[[idx[2]]]),
          intersect(high_abund_list[[idx[2]]], low_abund_list[[idx[1]]])
        )
      })
    ))
  } else {
    mixed_abundance <- character(0)
  }

  # Step 3: Assign tax_type
  data_classified <- data %>%
    dplyr::mutate(
  tax_type = dplyr::case_when(
  tax_val == unknown_label ~ "unknown",
  tax_val %in% special_taxa ~ "shared abundant",
  tax_val %in% mixed_abundance ~ "shared mixed abundance",
  tax_val %in% shared_low ~ "shared low abundant",
  tax_val %in% unique_taxa & RA >= low_abundance_threshold ~ "unique abundant",
  tax_val %in% unique_taxa & RA > 0 & RA < low_abundance_threshold ~ "unique low abundant",
  RA < low_abundance_threshold &
    !(tax_val %in% shared_low) &
    !(tax_val %in% mixed_abundance) ~ "unique low abundant",
  TRUE ~ "shared abundant"
)
    )

  # Step 4: Create category and tax_color
  abundant_shared <- Reduce(
    intersect,
    lapply(groups, function(g) {
      data_classified %>%
        dplyr::filter(
          !!group_sym == g,
          tax_type == "shared abundant"
        ) %>%
        dplyr::pull(tax_val) %>%
        unique()
    })
  )

  data_classified <- data_classified %>%
    dplyr::mutate(
      category = dplyr::case_when(
        tax_type == "shared mixed abundance" ~ "shared mixed abundance",
        tax_type == "shared abundant" ~ "shared abundant",
        tax_type == "shared low abundant" ~ "shared low abundant",
        tax_type == "unique abundant" & !(tax_val %in% abundant_shared) ~ "unique abundant",
        tax_type == "unique low abundant" ~ "unique low abundant",
        tax_type == "unknown" ~ "unknown",
        TRUE ~ "unclassified"
      ),
      # tax_color is a palette KEY
      tax_color = dplyr::case_when(
        tax_type == "shared low abundant" ~ "shared low abundant",
        tax_type == "unique low abundant" ~ "unique low abundant",
        tax_type == "shared mixed abundance" ~ tax_val,
        TRUE ~ tax_val
      )
    )

  # Step 5: Propagate unknown/low-abundant categories to higher taxonomy columns if present
  possible_tax_levels <- c("Domain", "Kingdom", "Phylum", "Class", "Order", "Family", "Genus", "Species")
  tax_cols_in_data <- intersect(possible_tax_levels, colnames(data_classified))

  if (length(tax_cols_in_data) > 1) {
    for (tax_col in tax_cols_in_data) {
      if (tax_col != tax_level) {
        data_classified <- data_classified %>%
          dplyr::mutate(
            !!tax_col := dplyr::case_when(
              tax_val == unknown_label ~ unknown_label,
              category == "shared low abundant" ~ "shared low abundant",
              category == "unique low abundant" ~ "unique low abundant",
              TRUE ~ .data[[tax_col]]
            )
          )
      }
    }
  }

  return(data_classified)
}
