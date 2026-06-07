#' Normalize Incertae sedis and fill missing taxonomy ranks
#'
#' Cleans hierarchical taxonomy columns by standardizing common "Incertae sedis"
#' variants, filling missing child ranks from stable parent ranks, and replacing
#' empty/uninformative entries with `"unknown"`.
#'
#' The function is designed for rank-ordered taxonomy (e.g., Domain → Phylum → Class → ...),
#' and performs a downward pass through the hierarchy to propagate parent information.
#'
#' @param df A data.frame or tibble containing taxonomy columns.
#' @param hierarchy Character vector of column names in rank order from highest
#'   (parent) to lowest (child). Only columns present in `df` are used.
#' @param clean_levels Character vector of ranks to normalize and fill. Only ranks
#'   present in both `df` and `hierarchy` are used.
#'
#' @details
#' The following rules are applied:
#' \itemize{
#'   \item "Incertae sedis" variants (case-insensitive; underscores/spaces tolerated)
#'         are converted to `"inc. sed."`.
#'   \item Empty strings (`""`) in `clean_levels` are treated as missing (`NA`).
#'   \item A parent is considered \emph{unstable} if it is missing (`NA`) or matches
#'         `"unknown"`, `"unclassified"`, or `"inc. sed."` (case-insensitive).
#'   \item If a child is `"inc. sed."` and the parent is stable, the child becomes
#'         `"<parent> inc. sed."`; otherwise it becomes `"unknown"`.
#'   \item If a child is missing and the parent is stable, the child becomes
#'         `"<parent>, unclassified"`; otherwise it becomes `"unknown"`.
#'   \item Underscores are converted to spaces and values are trimmed.
#' }
#'
#' @return An object of the same class as `df` (tibble in, tibble out) with cleaned
#'   taxonomy columns.
#'
#' @examples
#' tax_df <- data.frame(
#'   Domain = c("Bacteria", "Bacteria", NA),
#'   Phylum = c("Proteobacteria", "IncertaE_Sedis", "Actinobacteriota"),
#'   Class  = c(NA, "Alphaproteobacteria", ""),
#'   Order  = c("", "Rhizobiales", "inc. sed."),
#'   Family = c("Rhizobiaceae", NA, "unknown"),
#'   Genus  = c("", "inc. sed.", NA)
#' )
#'
#' replace_incertae_sedis_NAs(
#'   df = tax_df,
#'   hierarchy = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
#'   clean_levels = c("Phylum", "Class", "Order", "Family", "Genus")
#' )
#'
#' @export
replace_incertae_sedis_NAs <- function(df,
                                       hierarchy = c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
                                       clean_levels = c("Phylum", "Class", "Order", "Family", "Genus")) {
  is_tbl <- inherits(df, "tbl_df")
  hierarchy <- intersect(hierarchy, colnames(df))
  clean_levels <- intersect(clean_levels, hierarchy)

  # Normalize Incertae Sedis and empty → NA
  for (lvl in clean_levels) {
    df[[lvl]] <- gsub(
      "(_cls_|_ord_|_fam_|_gen_)?Incertae[_ ]?Sedis",
      "inc. sed.", df[[lvl]],
      ignore.case = TRUE
    )
    # Add space before inc. sed.
    df[[lvl]] <- gsub("(?<=\\S)inc\\. sed\\.", " inc. sed.", df[[lvl]], perl = TRUE)

    df[[lvl]][df[[lvl]] == ""] <- NA_character_

    # QIIME placeholders like Unknown_Family, Unknown_Genus, etc.
    df[[lvl]] <- dplyr::if_else(
      grepl("^Unknown(_|\\s)", df[[lvl]], ignore.case = TRUE),
      NA_character_,
      df[[lvl]]
    )
  }

  # Downward propagation
  for (i in 2:length(hierarchy)) {
    parent <- hierarchy[i - 1]
    child <- hierarchy[i]

    parent_val <- df[[parent]]
    child_val <- df[[child]]

    # Unstable parent: NA, unknown, unclassified, inc. sed.
    unstable_parent <- (parent %in% clean_levels) &
      (is.na(parent_val) | grepl("inc\\. sed\\.|unknown|unclassified", parent_val, ignore.case = TRUE))

    # Child inc. sed. → parent inc. sed. if parent stable, else unknown
    child_inc_sedis <- !is.na(child_val) & child_val == "inc. sed."
    df[[child]][child_inc_sedis & !unstable_parent] <- paste0(parent_val[child_inc_sedis & !unstable_parent], " inc. sed.")
    df[[child]][child_inc_sedis & unstable_parent] <- "unknown"

    # Child missing → propagate if parent stable, else unknown
    child_missing <- is.na(child_val)
    df[[child]][child_missing & !unstable_parent] <- paste0(parent_val[child_missing & !unstable_parent], ", unclassified")
    df[[child]][child_missing & unstable_parent] <- "unknown"
  }

  # Cleanup
  for (lvl in hierarchy) {
    df[[lvl]] <- gsub("_", " ", df[[lvl]])
    df[[lvl]] <- trimws(df[[lvl]])
    df[[lvl]] <- gsub("^unknown, unclassified$", "unknown", df[[lvl]])
    df[[lvl]] <- gsub("^, unclassified$", "unknown", df[[lvl]])
    df[[lvl]] <- gsub("  inc. sed.$", " inc. sed.", df[[lvl]])
  }

  if (is_tbl) df <- tibble::as_tibble(df)
  return(df)
}
