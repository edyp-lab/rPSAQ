

#' Deconvolute Peptide Abundances from MS2 Data
#'
#' Deconvolute peptide abundances by separating endogenous and PSAQ (Protein Standard Absolute Quantification)
#' contributions from MS2 isotopic data. This function processes experimental data to calculate corrected
#' intensities and areas for endogenous and PSAQ peptides.
#'
#' @param ms2_data A data frame containing MS2 fragments first and second isotopes abundances.
#'
#' @return A data frame containing deconvoluted peptide abundances with the following additional columns:
#'         \describe{
#'           \item{i1_endo_intensity}{Calculated intensity of the endogenous second isotope.}
#'           \item{I0_psaq_intensity}{Calculated intensity of the PSAQ first isotope.}
#'           \item{I0_endo_intensity}{Calculated intensity of the endogenous first isotope.}
#'           \item{experimental_ratio_intensity}{Ratio of PSAQ to endogenous intensities.}
#'           \item{i1_endo_area}{Calculated area of the endogenous second isotope.}
#'           \item{I0_psaq_area}{Calculated area of the PSAQ first isotope.}
#'           \item{I0_endo_area}{Calculated area of the endogenous first isotope.}
#'           \item{experimental_ratio_area}{Ratio of PSAQ to endogenous areas.}
#'         }
#'         The function returns NA for infinite values.
#'
#' @details
#' This function performs the following steps:
#' 1. Validates the required columns.
#' 2. Extracts and renames columns for first and second isotope data.
#' 3. Joins first and second isotope data.
#' 4. Generates theoretical isotopic distributions for the provided sequences.
#' 5. Merges experimental data with theoretical distributions.
#' 6. Calculates corrected intensities and areas for endogenous and PSAQ peptides.
#' 7. Computes experimental ratios for both intensities and areas.
#'
#' The input data frame must contain the following columns:
#' \describe{
#'   \item{seq}{Peptide sequence.}
#'   \item{mz}{Precursor m/z.}
#'   \item{frag_mz}{Fragment m/z.}
#'   \item{rt}{Retention time.}
#'   \item{m/z}{m/z value.}
#'   \item{frag_isotope}{Isotope label ("first" or "second").}
#'   \item{frag_label}{Fragment label (e.g., "y3").}
#'   \item{intensity}{One or more columns with "intensity" in their names.}
#'   \item{area}{One or more columns with "area" in their names.}
#' }
#'
#' @examples
#' # Deconvolute peptide abundances from MS2 data
#' data = deconvoluate_peptides_abundances(sample_10ng_r4)
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr full_join left_join select arrange mutate filter
#' @importFrom stringr str_split
#'
#' @seealso
#' \code{\link{compute_all_fragments_distributions}} for generating theoretical fragment distributions.
#'
#' @export
#'
deconvoluate_peptides_abundances = function(ms2_data) {

  required_columns = c("seq", "mz", "frag_mz", "rt", "m/z", "frag_isotope", "frag_label")

  #
  # Data check
  #
  intensity_column_idx = grep("intensity", names(ms2_data))
  area_column_idx = grep("area", names(ms2_data))

  data_check <- all(required_columns %in% names(ms2_data)) & length(intensity_column_idx) > 0 & length(area_column_idx) > 0

  if (!data_check) {
    print("The dataset must contain columns named (seq, mz, frag_mz, rt, m/z, frag_isotope, frag_label) and two columns suffixed by intensity and area")
    return()
  }

  # Extract sequences
  sequences = unique(ms2_data$seq)

  # Extract isotope 1 abundances and rename columns
  data_second_isotope = ms2_data %>% filter(tolower(.data$frag_isotope) == "second" ) %>% mutate(seq = stringr::str_split(.data$seq, " ", simplify = TRUE)[ ,1])
  data_second_isotope = data_second_isotope %>% rename(i1_frag_mz = .data$frag_mz)
  colnames(data_second_isotope)[intensity_column_idx] =  "i1_intensity"
  colnames(data_second_isotope)[area_column_idx] = "i1_area"

  # Extract isotope 0 abundances, rename columns and concatenate columns with initial data
  data_first_isotope = ms2_data %>% filter(tolower(.data$frag_isotope) == "first") %>% mutate(seq = stringr::str_split(.data$seq, " ", simplify = TRUE)[ ,1])
  data_first_isotope = data_first_isotope %>% rename(i0_frag_mz = .data$frag_mz)
  colnames(data_first_isotope)[intensity_column_idx] =  "i0_intensity"
  colnames(data_first_isotope)[area_column_idx] = "i0_area"

  data_ms2 = dplyr::full_join(data_first_isotope, data_second_isotope, by = c("seq", "frag_label"), suffix = c("_i0", "_i1"))
  data_ms2 = data_ms2 %>% arrange(seq, desc(.data$i0_intensity))

  # Generate the theoretical isotopic distribution of fragments when the selection window targets the second isotope
  theory = compute_all_fragments_distributions(sequences, fragment_type = c("y"))

  # Join by peptide sequence and fragment mz
  merge_ms2 = dplyr::left_join(data_ms2, theory, by = c("seq" = "precursor_seq", "frag_label" = "fragment_label"))

  # Compute ratios based on intensities
  merge_ms2 = merge_ms2 %>% mutate(i1_endo_intensity = .data$i0_intensity * .data$ms2_isotopic_ratio)
  merge_ms2 = merge_ms2 %>% mutate(I0_psaq_intensity = .data$i1_intensity - .data$i1_endo_intensity)
  merge_ms2 = merge_ms2 %>% mutate(I0_endo_intensity = (.data$i0_intensity + .data$i1_endo_intensity) / .data$ms1_isotopic_ratio)
  merge_ms2 = merge_ms2 %>% mutate(experimental_ratio_intensity =  .data$I0_psaq_intensity / .data$I0_endo_intensity)
  # Compute ratios based on area
  merge_ms2 = merge_ms2 %>% mutate(i1_endo_area = .data$i0_area * .data$ms2_isotopic_ratio)
  merge_ms2 = merge_ms2 %>% mutate(I0_psaq_area = .data$i1_area - .data$i1_endo_area)
  merge_ms2 = merge_ms2 %>% mutate(I0_endo_area = (.data$i0_area + .data$i1_endo_area) / .data$ms1_isotopic_ratio)
  merge_ms2 = merge_ms2 %>% mutate(experimental_ratio_area =  .data$I0_psaq_area / .data$I0_endo_area)

  # Replace Inf by NA
  merge_ms2[sapply(merge_ms2, is.infinite)] = NA

  merge_ms2 = merge_ms2 %>% select(all_of(names(data_ms2)), .data$experimental_ratio_intensity, .data$experimental_ratio_area) %>% arrange(.data$seq, .data$i0_frag_mz)

  return(merge_ms2)
}



#' Perform Batch PSAQ Analysis on Multiple Excel Files
#'
#' Perform batch processing of PSAQ (Protein Standard Absolute Quantification) analysis
#' on multiple Excel files matching a specific pattern. This function reads each file,
#' deconvolutes peptide abundances, and saves the results to new Excel files.
#'
#' @param dir A character string specifying the directory where the Excel files must be searched.
#'                  These sequences are used to generate theoretical isotopic distributions.
#' @param xlsx_files_pattern A character string specifying the pattern to match Excel files.
#'                           Example pattern : "^784.*\\\.xlsx$", which matches files starting with "784"
#'                           and ending with ".xlsx". This pattern is passed to \code{\link{list.files}}.
#'
#' @return This function does not return a value. Instead, it writes the results of the deconvolution
#'         to new Excel files in the current working directory. Each output file is named
#'         "calculated_ratios_[original_filename]".
#'
#' @details
#' This function performs the following steps for each Excel file matching the pattern:
#' 1. Identifies all Excel files in the current working directory that match the specified pattern.
#' 2. For each file, reads the file and call \code{\link{deconvoluate_peptides_abundances}} to process the data.
#' 3. Writes the results to a new Excel file with the prefix "calculated_ratios_" followed by the original filename.
#'
#' The function uses \code{\link{deconvoluate_peptides_abundances}} to deconvolute peptide abundances
#' by separating endogenous and PSAQ contributions from MS2 isotopic data.
#'
#' @examples
#' # Perform batch PSAQ analysis on files starting with "784" and ending with ".xlsx"
#' batch_psaq_analysis("inst/extdata", "^784.*\\.xlsx$")
#'
#' # Perform batch PSAQ analysis on files starting with "sample" and ending with ".xlsx"
#' batch_psaq_analysis("inst/extdata", "^sample.*\\.xlsx$")
#'
#' @importFrom openxlsx write.xlsx
#'
#' @seealso
#' \code{\link{deconvoluate_peptides_abundances}} for deconvoluting peptide abundances from a single Excel file.
#' \code{\link{list.files}} for listing files matching a pattern.
#'
#' @export
#'
batch_psaq_analysis = function(dir, xlsx_files_pattern) {
  files <- list.files(dir, pattern = xlsx_files_pattern)
  print(paste0("Found ", length(files), " files matching pattern ",xlsx_files_pattern, " in ", dir))
  for(filename in files) {
    print(paste0("Processing file : ", filename, " ..."))
    inputfile = file.path(dir, filename)
    ms2_data = openxlsx::read.xlsx(inputfile, sheet = 1, colNames = TRUE)
    psaq_deconvolution = deconvoluate_peptides_abundances(ms2_data)
    outputfile = file.path(dir, paste0("calculated_ratios_", filename))
    openxlsx::write.xlsx(psaq_deconvolution, outputfile)
  }
}

#' Plot PSAQ Experimental Ratios
#'
#' Create a scatter plot to visualize experimental PSAQ (Protein Standard Absolute Quantification) ratios
#' from deconvoluted peptide abundance data.
#'
#' @param psaq_deconvolution A data frame containing deconvoluted peptide abundance data.
#'                           This data frame should include the following columns:
#'                           \describe{
#'                             \item{experimental_ratio_area}{Experimental ratio of PSAQ to endogenous absed on area abundances.}
#'                             \item{experimental_ratio_intensity}{Experimental ratio of PSAQ to endogenous based on intensity abundances.}
#'                             \item{seq}{Peptide sequence for coloring points by sequence.}
#'                           }
#'
#' @return A ggplot object displaying experimental PSAQ ratios. The plot includes:
#'         \describe{
#'           \item{Circular points}{Representing experimental_ratio_area values.}
#'           \item{X-shaped points}{Representing experimental_ratio_intensity values (if present).}
#'         }
#'         The y-axis is limited to the range [-0.5, 2.5].
#'
#' @details
#' This function creates a scatter plot where:
#' \enumerate{
#'   \item The x-axis represents the row indices of the input data frame.
#'   \item The y-axis represents the experimental ratio values.
#'   \item Points are colored by peptide sequence.
#'   \item Circular points represent experimental_ratio_area values.
#'   \item X-shaped points represent experimental_ratio_intensity values (if this column exists in the data).
#'   \item The y-axis is constrained between -0.5 and 2.5 to focus on relevant ratio values.
#' }
#'
#' @examples
#' # Deconvoluate then plot PSAQ example data
#' psaq_data = deconvoluate_peptides_abundances(sample_10ng_r4)
#' plot = plot_psaq_ratios(psaq_data)
#'
#'
#' @importFrom ggplot2 ggplot geom_point ylim aes labs
#'
#' @seealso
#' \code{\link{deconvoluate_peptides_abundances}} for generating the input data.
#' \code{\link[ggplot2]{ggplot}} for creating the base plot.
#' \code{\link[ggplot2]{geom_point}} for adding points to the plot.
#'
#' @export
#'
plot_psaq_ratios = function(psaq_deconvolution) {
  plot = ggplot(psaq_deconvolution) +
    geom_point(aes(x = as.numeric(row.names(psaq_deconvolution)), y = .data$experimental_ratio_intensity, color = seq), shape = 4) +
    geom_point(aes(x = as.numeric(row.names(psaq_deconvolution)), y = .data$experimental_ratio_area, color = seq)) +
    ylim(-0.5, 2.5) + labs(x = "peptides/fragments", y = "ratio PSAQ/endo")
  return(plot)
}
