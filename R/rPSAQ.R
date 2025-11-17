

#' Deconvolute Peptide Abundances from MS2 Data
#'
#' Deconvolute peptide abundances by separating endogenous and PSAQ (Protein Standard Absolute Quantification)
#' contributions from MS2 isotopic data. This function processes experimental data to calculate
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
#'   \item{frag_isotope}{Isotope label ("first" or "second").}
#'   \item{frag_label}{Fragment label (e.g., "y3").}
#'   \item{intensity}{One or more columns with "intensity" in their names.}
#'   \item{area}{One or more columns with "area" in their names.}
#' }
#'
#' @examples
#' # Deconvolute peptide abundances from MS2 data
#' data = deconvolute_peptides_abundances(sample_10ng_r4)
#'
#' @importFrom openxlsx read.xlsx
#' @importFrom dplyr full_join left_join select arrange mutate filter
#' @importFrom stringr str_split
#'
#' @seealso
#' \code{\link{compute_all_fragments_distributions}} for generating theoretical fragment distributions.
#' \code{\link{reshape_data}} for reshaping abundance.
#'
#' @export
#'
deconvolute_peptides_abundances = function(ms2_data) {

  # Extract sequences
  sequences = unique(ms2_data$seq)

  data_ms2 = reshape_data(ms2_data)

  # Generate the theoretical isotopic distribution of fragments when the selection window targets the second isotope
  theory = compute_all_fragments_distributions(sequences, fragment_type = c("y", "b"))

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



#' Reshape MS2 Abundance Data
#'
#' Reshape MS2 data by separating first and second isotope measurements and renaming columns
#' for downstream analysis.
#'
#' @param ms2_data A data frame containing MS2 isotopic data with the following required columns:
#'        \describe{
#'          \item{seq}{Peptide sequence (may contain additional information separated by spaces).}
#'          \item{mz}{Precursor m/z value.}
#'          \item{frag_mz}{Fragment m/z value.}
#'          \item{rt}{Retention time.}
#'          \item{frag_isotope}{Isotope label, must contain values "first" or "second" (case-insensitive).}
#'          \item{frag_label}{Fragment label (e.g., "y3", "b5").}
#'        }
#'        Additionally, the data frame must contain at least one column with "intensity" in its name
#'        and at least one column with "area" in its name.
#'
#' @return A data frame with reshaped isotopic data containing:
#'         \describe{
#'           \item{seq}{Cleaned peptide sequence (first part before space if present).}
#'           \item{mz}{Precursor m/z value.}
#'           \item{i0_frag_mz}{Fragment m/z for first isotope measurements.}
#'           \item{i1_frag_mz}{Fragment m/z for second isotope measurements.}
#'           \item{i0_intensity}{Intensity for first isotope measurements.}
#'           \item{i1_intensity}{Intensity for second isotope measurements.}
#'           \item{i0_area}{Area for first isotope measurements.}
#'           \item{i1_area}{Area for second isotope measurements.}
#'           \item{frag_label}{Fragment label.}
#'           \item{rt}{Retention time.}
#'         }
#'         The returned data frame is sorted by peptide sequence and descending first isotope intensity.
#'         If input data doesn't meet requirements, the function returns NULL and prints an error message.
#'
#' @details
#' This function performs the following operations:
#' 1. Validates that the input data contains all required columns
#' 2. Extracts and processes first isotope measurements (frag_isotope = "first")
#' 3. Extracts and processes second isotope measurements (frag_isotope = "second")
#' 4. Cleans peptide sequences by taking the first part before any space
#' 5. Renames columns to indicate isotope origin (i0_ for first isotope, i1_ for second isotope)
#' 6. Joins first and second isotope data by peptide sequence and fragment label
#' 7. Sorts the resulting data frame by peptide sequence and descending first isotope intensity
#'
#' @examples
#' # Example data frame with required columns
#' data <- data.frame(
#'   seq = c("PEPTIDE1 first", "PEPTIDE1 second", "PEPTIDE2 first", "PEPTIDE2 second"),
#'   mz = c(500, 500, 600, 600),
#'   frag_mz = c(200, 200.5, 300, 300.5),
#'   rt = c(10, 10, 15, 15),
#'   frag_isotope = c("first", "second", "first", "second"),
#'   frag_label = c("y3", "y3", "y4", "y4"),
#'   peptide_intensity = c(1000, 500, 1500, 750),
#'   peptide_area = c(5000, 2500, 7500, 3750)
#' )
#'
#' # Reshape the data
#' reshaped_data <- reshape_data(data)
#' head(reshaped_data)
#'
#' @importFrom dplyr filter mutate rename full_join arrange
#' @importFrom stringr str_split
#'
#'
#' @export
#'
reshape_data = function(ms2_data) {
  required_columns = c("seq", "mz", "frag_mz", "rt", "frag_isotope", "frag_label")

  #
  # Data check
  #
  intensity_column_idx = grep("intensity", names(ms2_data))
  area_column_idx = grep("area", names(ms2_data))

  data_check <- all(required_columns %in% names(ms2_data)) & length(intensity_column_idx) > 0 & length(area_column_idx) > 0

  if (!data_check) {
    print("The dataset must contain columns named (seq, mz, frag_mz, rt, frag_isotope, frag_label) and two columns suffixed by intensity and area")
    return()
  }


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

  return(data_ms2)
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
#' @param output_plot A boolean specifying if graphical representation must be generated. If TRUE
#'                        two plot (one based on area abundances, one based on intensity abundances) will
#'                        be generated for each input file
#'
#' @return This function does not return a value. Instead, it writes the results of the deconvolution
#'         to new Excel files in the current working directory. Each output file is named
#'         "calculated_ratios_[original_filename]".
#'
#' @details
#' This function performs the following steps for each Excel file matching the pattern:
#' 1. Identifies all Excel files in the current working directory that match the specified pattern.
#' 2. For each file, reads the file and call \code{\link{deconvolute_peptides_abundances}} to process the data.
#' 3. Writes the results to a new Excel file with the prefix "calculated_ratios_" followed by the original filename.
#'
#' The function uses \code{\link{deconvolute_peptides_abundances}} to deconvolute peptide abundances
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
#' @importFrom ggplot2 labs
#' @importFrom grDevices dev.off pdf
#' @importFrom utils, read.table
#'
#' @seealso
#' \code{\link{deconvolute_peptides_abundances}} for deconvoluting peptide abundances from a single Excel file.
#' \code{\link{list.files}} for listing files matching a pattern.
#'
#' @export
#'
batch_psaq_analysis = function(dir, xlsx_files_pattern, output_plot = FALSE) {
  files <- list.files(dir, pattern = xlsx_files_pattern)
  print(paste0("Found ", length(files), " files matching pattern ",xlsx_files_pattern, " in ", dir))
  if (output_plot) {
    grDevices::pdf(file.path(dir, "plots.pdf"))
  }
  for(filename in files) {
    print(paste0("Processing file : ", filename, " ..."))
    inputfile = file.path(dir, filename)
    if (endsWith(filename, "xlsx")) {
      ms2_data = openxlsx::read.xlsx(inputfile, sheet = 1, colNames = TRUE)
    } else if (endsWith(filename, "csv")) {
      ms2_data = utils::read.table(file = inputfile, sep = ";", dec = ".", header = TRUE)
    }
    psaq_deconvolution = deconvolute_peptides_abundances(ms2_data)
    if (output_plot) {
      plot_area = plot_psaq_ratios(psaq_deconvolution, "area")
      plot_area = plot_area + labs(title = paste0(filename, " (Area) "))
      print(plot_area)
      plot_intensity = plot_psaq_ratios(psaq_deconvolution, "intensity")
      plot_intensity = plot_intensity + labs(title = paste0(filename, " (Intensity) "))
      print(plot_intensity)
    }
    if (output_plot) {
      grDevices::dev.off()
    }
    outputfile = file.path(dir, paste0("calculated_ratios_", filename))
    openxlsx::write.xlsx(psaq_deconvolution, outputfile)
  }
}

#' Plot PSAQ Experimental Ratios
#'
#' Creates a scatter plot to visualize experimental PSAQ (Protein Standard Absolute Quantification) ratios
#' from deconvoluted peptide abundance data.
#'
#' @param psaq_deconvolution A data frame containing deconvoluted peptide abundance data.
#'                           This data frame should include the following columns:
#'                           \describe{
#'                             \item{experimental_ratio_area}{Experimental ratio of PSAQ to endogenous absed on area abundances.}
#'                             \item{experimental_ratio_intensity}{Experimental ratio of PSAQ to endogenous based on intensity abundances.}
#'                             \item{seq}{Peptide sequence for coloring points by sequence.}
#'                           }
#' @param base Plot intensity or area based ratios. Possible values are "intensity" or "area".
#'
#' @return A ggplot object displaying experimental PSAQ ratios. The plot includess the median experimental ratio
#' represented as an horizontal line. The y-axis is limited to the range [-0.5, 2.5].
#'
#' @details
#' This function creates a scatter plot where:
#' \enumerate{
#'   \item The x-axis represents the row indices of the input data frame.
#'   \item The y-axis represents the experimental ratio values.
#'   \item Points are colored by peptide sequence.
#'   \item The y-axis is constrained between -0.5 and 2.5 to focus on relevant ratio values.
#' }
#'
#' @examples
#' # Deconvolute then plot PSAQ example data
#' psaq_data = deconvolute_peptides_abundances(sample_10ng_r4)
#' plot_psaq_ratios(psaq_data, "area")
#'
#'
#' @importFrom ggplot2 ggplot geom_point ylim aes labs geom_segment annotate
#' @importFrom stats median
#'
#' @seealso
#' \code{\link{deconvolute_peptides_abundances}} for generating the input data.
#' \code{\link[ggplot2]{ggplot}} for creating the base plot.
#' \code{\link[ggplot2]{geom_point}} for adding points to the plot.
#'
#' @export
#'
plot_psaq_ratios = function(psaq_deconvolution, base = "area") {
  if (tolower(base) == "area") {
    med = median(psaq_deconvolution$experimental_ratio_area, na.rm = T)
    sd = sd(psaq_deconvolution$experimental_ratio_area, na.rm = T)
    plot = ggplot(psaq_deconvolution) +
      geom_point(aes(x = as.numeric(row.names(psaq_deconvolution)), y = .data$experimental_ratio_area, color = seq)) +
      ylim(-0.5, 2.5) + labs(x = "peptides/fragments", y = "ratio PSAQ/endo (area)") +
      geom_segment(x = 0, y = med, xend = nrow(psaq_deconvolution)+1, yend = med) +
      annotate("text", x=0, y=2.5, label= paste0("median ratio = ", round(med,3)), size = 3.5, hjust = 0) +
      annotate("text", x=0, y=2.4, label= paste0("std devation = ", round(sd,4)), size = 3.5, hjust = 0)
  } else{
    med = median(psaq_deconvolution$experimental_ratio_intensity, na.rm = T)
    sd = sd(psaq_deconvolution$experimental_ratio_intensity, na.rm = T)
    plot = ggplot(psaq_deconvolution) +
      geom_point(aes(x = as.numeric(row.names(psaq_deconvolution)), y = .data$experimental_ratio_intensity, color = seq)) +
      ylim(-0.5, 2.5) + labs(x = "peptides/fragments", y = "ratio PSAQ/endo (intensity)") +
      geom_segment(x = 0, y = med, xend = nrow(psaq_deconvolution)+1, yend = med) +
      annotate("text", x=0, y=2.5, label= paste0("median ratio = ", round(med,3)), size = 3.5, hjust = 0) +
      annotate("text", x=0, y=2.4, label= paste0("std devation = ", round(sd,4)), size = 3.5, hjust = 0)

  }

  return(plot)
}
