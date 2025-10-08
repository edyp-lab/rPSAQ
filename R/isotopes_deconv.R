#' Calculate the theoretical isotopic distribution of a peptide
#'
#' This function computes the theoretical isotopic distribution of a precursor from its amino acids sequence.
#'  The distribution represents the relative abundance of isotopologues (molecules with different isotopic compositions)
#'  by taking into account the natural abundance  of heavy isotopes (^13C, ^2H, ^15N, ^18O, and ^33S).
#' It is useful for predicting isotope patterns in mass spectrometry.
#'
#' @param sequence A character string representing the amino acid sequence of the peptide (e.g., "ALQASALAAWGGK").
#' @param size An integer specifying the number of isotopologues to compute (default = 4).
#'             This corresponds to the number of isotope peaks to model, starting from the monoisotopic peak.
#'
#' @return A data frame with two columns:
#'         \describe{
#'           \item{proba}{Raw probabilities of each isotopologue.}
#'           \item{percent}{Relative abundance of each isotopologue as a percentage of the monoisotopic peak.}
#'         }
#'
#' @details
#' The function uses the \code{\link[OrgMassSpecR]{ConvertPeptide}} function to convert the peptide sequence
#' into its elemental composition, and the \code{\link[sinib]{dsinib}} function to compute the isotopic distribution.
#'
#' @examples
#' # Example usage
#' isotopic_distribution("ALQASALAAWGGK", size = 4)
#'
#' @importFrom OrgMassSpecR ConvertPeptide
#' @importFrom sinib dsinib
#' @import dplyr
#'
#' @seealso
#' \code{\link[OrgMassSpecR]{ConvertPeptide}} for converting peptide sequences to elemental compositions.
#' \code{\link[sinib]{dsinib}} for computing isotopic distributions.
#'
#' @export
#'
isotopic_distribution = function(sequence, size = 4) {

  proba_isotop = c(0.0107, 0.000184, 0.00364, 0.000380781, 0.007948959) # occurrence of the second isotope of each atom

  formula = OrgMassSpecR::ConvertPeptide(sequence, output = "elements", IAA = FALSE)
  n_espece = c(formula$C, formula$H, formula$N, formula$O, formula$S)
  isotopes = c(0:size-1)
  proba = sinib::dsinib(isotopes, size = n_espece, prob = proba_isotop)
  dist = data.frame(proba)
  dist = dist %>% mutate(percent = 100*proba/proba[1])

  return(dist)

}


#' Enumerate y and/or b fragments of a peptide
#'
#' This function enumerates y and/or b fragments of a peptide sequence, names each fragment,
#' computes their theoretical mass and atomic composition.
#'
#' @param sequence A character string representing the amino acid sequence of the peptide (e.g., `"ALQASALAAWGGK"`).
#' @param fragment_type A character vector specifying the types of fragments to generate.
#'                      Possible values are `"y"`, `"b"`, or `c("y", "b")` (default is `"y"`).
#'                      If `"y"`, y-ion fragments are generated. If `"b"`, b-ion fragments are generated.
#'                      If both are specified, both y-ion and b-ion fragments are generated.
#'
#' @return A data frame with the following columns:
#'         \describe{
#'           \item{id}{Integer representing the fragment number.}
#'           \item{precursor_seq}{Character string of the precursor peptide sequence.}
#'           \item{precursor_mass}{Double representing the monoisotopic mass of the precursor peptide.}
#'           \item{fragment_label}{Character string representing the label of the fragment (e.g., `"y3"`, `"b5"`).}
#'           \item{fragment_seq}{Character string of the fragment sequence.}
#'           \item{C}{Integer representing the number of carbon atoms in the fragment.}
#'           \item{H}{Integer representing the number of hydrogen atoms in the fragment.}
#'           \item{N}{Integer representing the number of nitrogen atoms in the fragment.}
#'           \item{O}{Integer representing the number of oxygen atoms in the fragment.}
#'           \item{S}{Integer representing the number of sulfur atoms in the fragment.}
#'           \item{fragment_mass}{Double representing the monoisotopic mass of the fragment.}
#'         }
#'
#' @details
#' The function uses the `ConvertPeptide` function to convert peptide sequences into their elemental compositions
#' and the `MonoisotopicMass` function to compute the monoisotopic mass of each fragment.
#' The precursor mass is computed with a charge of 2, while the fragment masses are computed with a charge of 1.
#' The function does not consider IAA (iodoacetamide) modification.
#'
#' @examples
#' # Generate y-ion fragments for a peptide
#' generate_peptide_fragments("ALQASALAAWGGK")
#'
#' # Generate both y-ion and b-ion fragments for a peptide
#' generate_peptide_fragments("ALQASALAAWGGK", fragment_type = c("y", "b"))
#'
#' @importFrom OrgMassSpecR ConvertPeptide MonoisotopicMass
#'
#' @seealso
#' \code{\link[OrgMassSpecR]{ConvertPeptide}} for converting peptide sequences to elemental compositions.
#' \code{\link[OrgMassSpecR]{MonoisotopicMass}} for computing monoisotopic masses.
#'
#' @export
#'
generate_peptide_fragments = function(sequence, fragment_type = c("y")) {

fragments = data.frame(id = integer(),
                       precursor_seq = character(),
                       precursor_mass = double(),
                       fragment_label = character(),
                       fragment_seq = character(),
                       C = integer(),
                       H = integer(),
                       N = integer(),
                       O = integer(),
                       S  = integer(),
                       fragment_mass = double())

prec_mass = OrgMassSpecR::MonoisotopicMass(ConvertPeptide(sequence, output = "elements", IAA = FALSE), charge = 2)

if (is.element("y", fragment_type)) {
  for(i in nchar(sequence):2) {
    f_sequence = substr(sequence, i, nchar(sequence))
    f_number = nchar(sequence)-i+1
    f_formula =  OrgMassSpecR::ConvertPeptide(f_sequence, output = "elements", IAA = FALSE)
    f_monoMass = OrgMassSpecR::MonoisotopicMass(f_formula, charge = 1)
    fragments = rbind(fragments,
                      data.frame(id = f_number,
                                 precursor_seq = sequence,
                                 precursor_mass = prec_mass,
                                 fragment_label = paste0("y", f_number),
                                 fragment_seq = f_sequence,
                                 C = f_formula$C,
                                 H = f_formula$H,
                                 N = f_formula$N,
                                 O = f_formula$O,
                                 S = f_formula$S,
                                 fragment_mass = f_monoMass))
  }
}

if (is.element("b", fragment_type)) {
  for(i in 2:nchar(sequence)-1) {
    f_sequence = substr(sequence, 1, i)
    f_number = i
    f_formula =  OrgMassSpecR::ConvertPeptide(f_sequence, output = "elements", IAA = FALSE)
    f_monoMass = OrgMassSpecR::MonoisotopicMass(f_formula, charge = 1)
    fragments = rbind(fragments,
                      data.frame(id = f_number,
                                 precursor_seq = sequence,
                                 precursor_mass = prec_mass,
                                 fragment_label = paste0("b", f_number),
                                 fragment_seq = f_sequence,
                                 C = f_formula$C,
                                 H = f_formula$H,
                                 N = f_formula$N,
                                 O = f_formula$O,
                                 S = f_formula$S,
                                 fragment_mass = f_monoMass))

  }
}
return(fragments)

}

#' Compute the Isotopic Distribution of an Fragment
#'
#' Compute the isotopic distribution of the specified fragment, taking into account that only
#' the second isotope of the precursor have been included in the selection window.
#' This function uses a hypergeometric distribution to model the contribution of precursor isotopes
#' to the fragment isotopic distribution.
#'
#' @param fragment_seq A character string representing the amino acid sequence of the fragment.
#' @param precursor_seq A character string representing the amino acid sequence of the precursor peptide.
#' @param size An integer specifying the number of isotopes to consider.
#'
#' @return A data frame with the following columns:
#'         \describe{
#'           \item{proba}{A numeric vector containing the probability associated with each MS2 isotope.}
#'         }
#'         The probabilities are computed using a hypergeometric distribution, where the probability
#'         of each isotope is based on the number of carbon-13 atoms included in the fragment.
#'
#' @details
#' The function uses a hypergeometric distribution to model the contribution of precursor isotopes
#' to the fragment isotopic distribution. The hypergeometric distribution is defined as follows:
#' \deqn{dhyper(x, m, n, k)}
#' where:
#' \describe{
#'   \item{x}{Number of "white balls" drawn (number of carbon-13 atoms in the fragment isotope).}
#'   \item{m}{Number of "white balls" in the urn (1, representing the carbon-13 atoms in the precursor).}
#'   \item{n}{Number of "black balls" in the urn (number of carbon-12 atoms in the precursor minus 1).}
#'   \item{k}{Number of balls drawn from the urn (number of carbon atoms in the fragment).}
#' }
#' In this context:
#' \describe{
#'   \item{White balls}{Number of carbon-13 atoms (MS1 isotope rank).}
#'   \item{Black balls}{Number of carbon-12 atoms (total carbon in the precursor minus the isotope rank).}
#' }
#'
#' @examples
#' # Compute the isotopic distribution for a fragment of the precursor peptide
#' ms2_isotopic_distribution("LAAWGGK", "ALQASALAAWGGK", 4)
#' ms2_isotopic_distribution("WGGK", "ALQASALAAWGGK", 4)
#' ms2_isotopic_distribution("ALQASALAAWGGK", "ALQASALAAWGGK", 4)
#'
#' @importFrom OrgMassSpecR ConvertPeptide
#' @importFrom stats dhyper
#'
#' @seealso
#' \code{\link[OrgMassSpecR]{ConvertPeptide}} for converting peptide sequences to elemental compositions.
#'
#' @export
#'
ms2_isotopic_distribution = function(fragment_seq, precursor_seq, size) {

  fragment_formula = OrgMassSpecR::ConvertPeptide(fragment_seq, output = "elements", IAA = FALSE)
  precursor_formula = OrgMassSpecR::ConvertPeptide(precursor_seq, output = "elements", IAA = FALSE)
  indexes = c(1:size)
  h = dhyper((indexes-1), 1, precursor_formula$C-1, fragment_formula$C)
  isotopes = data.frame(proba = h)
  return(isotopes)

}

#' Compute the Isotopic Distribution of MS2 Fragments of a Precursor
#'
#' Compute the isotopic distribution of y (and/or b) fragments of the specified precursor peptide sequence.
#' This function calculates the isotopic ratios for both MS1 and MS2 levels and appends them to the fragment data.
#'
#' @param sequence A character string representing the amino acid sequence of the peptide (e.g., "ALQASALAAWGGK").
#' @param fragment_type A character vector specifying the types of fragments to generate.
#'                      Possible values are "y", "b", or c("y", "b") (default is "y").
#'                      If "y", y-ion fragments are generated. If "b", b-ion fragments are generated.
#'                      If both are specified, both y-ion and b-ion fragments are generated.
#'
#' @return A data frame containing the following columns:
#'         \describe{
#'           \item{id}{Integer representing the fragment number.}
#'           \item{precursor_seq}{Character string of the precursor peptide sequence.}
#'           \item{precursor_mass}{Double representing the monoisotopic mass of the precursor peptide.}
#'           \item{fragment_label}{Character string representing the label of the fragment (e.g., "y3", "b5").}
#'           \item{fragment_seq}{Character string of the fragment sequence.}
#'           \item{C, H, N, O, S}{Integers representing the number of each type of atom in the fragment.}
#'           \item{fragment_mass}{Double representing the monoisotopic mass of the fragment.}
#'           \item{ms1_isotopic_ratio}{Double representing the isotopic ratio (M+1/M) at the MS1 level.}
#'           \item{ms2_isotopic_ratio}{Double representing the isotopic ratio (M+1/M) at the MS2 level for each fragment.}
#'         }
#'
#' @details
#' This function performs the following steps:
#' 1. Generates the specified fragments (y and/or b ions) using the `generate_peptide_fragments` function.
#' 2. Computes the isotopic distribution of the precursor peptide using the `isotopic_distribution` function.
#' 3. Calculates the MS1 isotopic ratio (M+1/M) and adds it to the fragment data.
#' 4. For each fragment, computes the MS2 isotopic distribution using the `ms2_isotopic_distribution` function.
#' 5. Calculates the MS2 isotopic ratio (M+1/M) for each fragment and appends it to the fragment data.
#'
#' The isotopic ratio is defined as the ratio of the probability of the first isotopologue (M+1) to the monoisotopic peak (M).
#'
#' @examples
#' # Compute the isotopic distribution for y-ion fragments of a peptide
#' compute_fragments_isotopic_distribution("ALQASALAAWGGK", c("y"))
#'
#' # Compute the isotopic distribution for b-ion fragments of a peptide
#' compute_fragments_isotopic_distribution("ALQASALAAWGGK", c("b"))
#'
#' @seealso
#' \code{\link{generate_peptide_fragments}} for generating peptide fragments.
#' \code{\link{isotopic_distribution}} for computing the isotopic distribution of a peptide.
#' \code{\link{ms2_isotopic_distribution}} for computing the isotopic distribution of MS2 fragments.
#'
#'
#' @export
#'
compute_fragments_isotopic_distribution = function(sequence, fragment_type = c("y")) {

  added = NULL
  fragments = generate_peptide_fragments(sequence, fragment_type = fragment_type)
  ms1_isotopic_dist = isotopic_distribution(sequence)
  fragments$ms1_isotopic_ratio = ms1_isotopic_dist$percent[2] / ms1_isotopic_dist$percent[1]

  for (i in 1 : nrow(fragments)) {
    ms2_isotopic_dist = ms2_isotopic_distribution(fragments[i,"fragment_seq"], sequence, 4)
    ms2_isotopic_ratio = ms2_isotopic_dist$proba[2] / ms2_isotopic_dist$proba[1]

    added = rbind(added, data.frame(ms2_isotopic_ratio))
  }

  fragments = cbind(fragments, added)

  return(fragments)
}

#' Generate Isotopic Distributions of Fragments from Peptide Sequences
#'
#' Generate theoretical fragments (y and/or b ions) and their isotopic distributions
#' for a list of peptide sequences.
#'
#' @param sequences A character vector containing one or more peptide sequences
#'                  (e.g., c("ALQASALAAWGGK", "PEPTIDEK")).
#' @param fragment_type A character vector specifying the types of fragments to generate.
#'                      Possible values are "y", "b", or c("y", "b") (default is "y").
#'                      If "y", y-ion fragments are generated. If "b", b-ion fragments are generated.
#'                      If both are specified, both y-ion and b-ion fragments are generated.
#'
#' @return A data frame containing the combined theoretical fragments and their isotopic distributions
#'         for all input sequences. The data frame includes the following columns:
#'         \describe{
#'           \item{id}{Integer representing the fragment number.}
#'           \item{precursor_seq}{Character string of the precursor peptide sequence.}
#'           \item{precursor_mass}{Double representing the monoisotopic mass of the precursor peptide.}
#'           \item{fragment_label}{Character string representing the label of the fragment (e.g., "y3", "b5").}
#'           \item{fragment_seq}{Character string of the fragment sequence.}
#'           \item{C, H, N, O, S}{Integers representing the number of each type of atom in the fragment.}
#'           \item{fragment_mass}{Double representing the monoisotopic mass of the fragment.}
#'           \item{ms1_isotopic_ratio}{Double representing the isotopic ratio (M+1/M) at the MS1 level.}
#'           \item{ms2_isotopic_ratio}{Double representing the isotopic ratio (M+1/M) at the MS2 level for each fragment.}
#'         }
#'
#' @details
#' This function iterates through a list of peptide sequences and generates theoretical fragments
#' with their isotopic distributions for each sequence. It combines the results into a single data frame.
#' The function uses \code{\link{compute_fragments_isotopic_distribution}} to calculate the fragments
#' and their isotopic distributions for each individual sequence.
#'
#' @examples
#' # Generate y-ion fragments and their isotopic distributions for multiple peptides
#' compute_all_fragments_distributions(c("ALQASALAAWGGK", "PEPTIDEK"))
#'
#' # Generate both y-ion and b-ion fragments for multiple peptides
#' compute_all_fragments_distributions(c("ALQASALAAWGGK", "PEPTIDEK"), fragment_type = c("y", "b"))
#'
#' @seealso
#' \code{\link{compute_fragments_isotopic_distribution}} for computing fragments and isotopic distributions
#' of a single peptide sequence.
#'
#'
#' @export
#'
compute_all_fragments_distributions  <- function(sequences, fragment_type = c("y")) {
  theoretical <- NULL
  for(sequence in sequences) {
    fragments_theory <- compute_fragments_isotopic_distribution(sequence, fragment_type = fragment_type)
    theoretical <- rbind(theoretical, fragments_theory)
  }

  return(theoretical)
}


