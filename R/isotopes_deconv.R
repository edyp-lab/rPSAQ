

#' The Atoms Natural Light and Heavy Isotopes
#'
#' @returns a named list of isotopic distribution of light and heavy (light+1) isotopes
#' @export
#'
isotopes_abundances = function() {
  return(list(C = c(0.9893, 0.0107), H = c(0.999885, 0.000115), N = c(0.996360, 0.003640), O = c(0.997570, 0.000380), S = c(0.9499, 0.0075)))
}

#' Calculate Isotope Probabilities Using Binomial Distribution
#'
#' Calculate the probabilities of observing different numbers of heavy isotopes in a molecule
#' using the probabilities distribution.
#'
#' @param x A vector of integers representing the number of heavy isotopes (k) to calculate probabilities for.
#'         Typically a sequence from 0 to the maximum number of heavy isotopes of interest.
#' @param n An integer representing the total number of atoms that can have heavy isotopes
#'         (e.g., number of carbon atoms in a molecule when calculating ^13C probabilities).
#' @param prob A numeric vector of length 2 containing:
#'        \describe{
#'          \item{prob[1]}{Probability of the light isotope (e.g., probability of ^12C).}
#'          \item{prob[2]}{Probability of the heavy isotope (e.g., probability of ^13C).}
#'        }
#'
#'
#' @return A numeric vector of probabilities corresponding to each value in x, calculated using the binomial formula:
#'         P(X = k) = C(n, k) * p^k * (q)^(n-k), where p is the probability of the heavy isotope and q the probability
#'         of the light isotope.
#'
#' @details
#' This function calculates isotope probabilities using the binomial probability mass function:
#' P(X = k) = C(n, k) * p^k * (q)^(n-k)
#' where:
#' \describe{
#'   \item{C(n, k)}{Binomial coefficient "n choose k"}
#'   \item{p}{prob[2], probability of the heavy isotope}
#'   \item{q}{prob[1], probability of the light isotope}
#'   \item{n}{Total number of atoms}
#'   \item{k}{Number of heavy isotopes (values in x)}
#' }
#'
#' The function is particularly useful for modeling isotopic distributions in mass spectrometry,
#' where it can predict the relative abundances of different isotopologues.
#'
#' @examples
#' # Calculate probabilities for 0 to 3 ^13C atoms in a molecule with 10 carbon atoms
#' # Natural abundance of ^13C is ~1.07%
#' x <- 0:3
#' n <- 10
#' prob <- c(0.9893, 0.0107)  # Probabilities for ^12C and ^13C
#' isotopes_probabilities(x, n, prob)
#'
#' # Calculate probabilities for 0 to 5 ^15N atoms in a molecule with 20 nitrogen atoms
#' # Natural abundance of ^15N is ~0.364%
#' x <- 0:5
#' n <- 20
#' prob <- c(0.99636, 0.00364)  # Probabilities for ^14N and ^15N
#' isotopes_probabilities(x, n, prob)
#'
#' @seealso
#' \code{\link{choose}} for calculating binomial coefficients.
#'
#'
#' @export
#'
isotopes_probabilities = function(x, n, prob) {
  probabilities = sapply(x, function(k) { choose(n, k)*prob[2]^k*prob[1]^(n-k) })
  return (probabilities)
}


#' Calculate the theoretical isotopic distribution of the first two isotopes of a peptide
#'
#' This function computes the theoretical isotopic distribution of a precursor from its amino acids sequence.
#'  The distribution represents the relative abundance of isotopologues (molecules with different isotopic compositions)
#'  by taking into account the natural abundance of +1 isotopes (^13C, ^2H, ^15N, ^17O, and ^33S).
#'
#' @param sequence A character string representing the amino acid sequence of the peptide (e.g., "ALQASALAAWGGK").
#'
#' @return A data frame with two columns:
#'         \describe{
#'           \item{proba}{Raw probabilities of each isotopologue.}
#'         }
#'
#' @details
#' The function uses the \code{\link[OrgMassSpecR]{ConvertPeptide}} function to convert the peptide sequence
#' into its elemental composition, and the \code{\link[sinib]{dsinib}} function to compute the isotopic distribution.
#'
#' @examples
#' # Example usage
#' isotopic_distribution("ALQASALAAWGGK")
#' isotopic_distribution("GILAADESVGTMGNR")
#' isotopic_distribution("LSFSYGR")
#'
#' @importFrom OrgMassSpecR ConvertPeptide
#' @import dplyr
#'
#' @seealso
#' \code{\link[OrgMassSpecR]{ConvertPeptide}} for converting peptide sequences to elemental compositions.
#'
#' @export
#'
isotopic_distribution = function(sequence) {

  proba_isotopes = isotopes_abundances()
  atoms= names(proba_isotopes)

  formula = OrgMassSpecR::ConvertPeptide(sequence, output = "elements", IAA = FALSE)
  formula[["H"]] = formula[["H"]] + 2 #doubly charged ion
  d = list()

  for (a in atoms) {
    if (formula[[a]] > 0) {
#      d[[a]] = dbinom(c(0:1), formula[[a]], proba_isotop[[a]])
      d[[a]] = isotopes_probabilities(c(0:1), formula[[a]], proba_isotopes[[a]])
    }
  }

  dist = data.frame(matrix(unlist(d), nrow=length(d), byrow=TRUE))

  i1 = prod(dist$X1)
  i2 = 0
  for (i in 1:length(atoms)) {
    if (formula[[i]] > 0) {
      i2 = i2+dist[i,2]*prod(dist[-i, 1], na.rm = TRUE)
    }
  }

  return(data.frame(proba = c(i1, i2)))

}

#' Compute the Isotopic Distribution of the first two isotopes of a Fragment
#'
#' Compute the isotopic distribution of the specified fragment, taking into account that only
#' the second isotope of the precursor have been included in the selection window.
#'
#' @param fragment_seq A character string representing the amino acid sequence of the fragment.
#' @param precursor_seq A character string representing the amino acid sequence of the precursor peptide.
#'
#' @return A data frame with the following columns:
#'         \describe{
#'           \item{proba}{A numeric vector containing the probability associated with each MS2 isotope.}
#'         }
#'         The probabilities of each isotope is based on the number of heavy atoms included in the fragment.
#'
#'
#' @examples
#' # Compute the isotopic distribution for a fragment of the precursor peptide
#' ms2_isotopic_distribution("LAAWGGK", "ALQASALAAWGGK")
#' ms2_isotopic_distribution("WGGK", "ALQASALAAWGGK")
#' ms2_isotopic_distribution("ALQASALAAWGGK", "ALQASALAAWGGK")
#' ms2_isotopic_distribution("SYGR", "LSFSYGR")
#'
#' @importFrom OrgMassSpecR ConvertPeptide
#'
#' @seealso
#' \code{\link[OrgMassSpecR]{ConvertPeptide}} for converting peptide sequences to elemental compositions.
#'
#' @export
#'
ms2_isotopic_distribution = function(fragment_seq, precursor_seq) {

  fragment_formula = OrgMassSpecR::ConvertPeptide(fragment_seq, output = "elements", IAA = FALSE)
  precursor_formula = OrgMassSpecR::ConvertPeptide(precursor_seq, output = "elements", IAA = FALSE)
  precursor_formula[["H"]] = precursor_formula[["H"]] + 2 # doubly charged precursor ion
  fragment_formula[["H"]] = fragment_formula[["H"]] + 1 # singly charged fragment ion

  proba_isotopes = isotopes_abundances()
  atoms= names(proba_isotopes)
  p = list()
  d = list()
  for (a in atoms) {
    if (precursor_formula[[a]] > 0) {
      d[[a]] = isotopes_probabilities(c(0:1), precursor_formula[[a]], proba_isotopes[[a]])
      p[[a]] =  c(1-(fragment_formula[[a]]/precursor_formula[[a]]), (fragment_formula[[a]]/precursor_formula[[a]]))
    }
  }

  dist = data.frame(matrix(unlist(d), nrow=length(d), byrow=TRUE))
  d = list()
  for (i in 1:length(atoms)) {
    if (precursor_formula[[i]] > 0) {
      d[[atoms[i]]] = dist[i,2]*prod(dist[-i, 1], na.rm = TRUE)
    }
  }

  d = as.list(unlist(d)/sum(unlist(d)))

  h1 = 0
  h2 = 0
  for (a in atoms) {
    if (precursor_formula[[a]] > 0) {
      h1 = h1 + d[[a]]*p[[a]][1]
      h2 = h2 + d[[a]]*p[[a]][2]
    }
  }

  isotopes = data.frame(proba = c(h1, h2))
  return(isotopes)

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
  fragments$ms1_isotopic_ratio = ms1_isotopic_dist$proba[2] / ms1_isotopic_dist$proba[1]
  fragments$ms1_isotopic_norm_ratio = ms1_isotopic_dist$proba[2] / (ms1_isotopic_dist$proba[1]+ms1_isotopic_dist$proba[2])

    for (i in 1 : nrow(fragments)) {

    ms2_isotopic_dist = ms2_isotopic_distribution(fragments[i,"fragment_seq"], sequence)
    ms2_isotopic_ratio = ms2_isotopic_dist$proba[2] / ms2_isotopic_dist$proba[1]
    ms2_isotopic_norm_ratio = ms2_isotopic_dist$proba[2] / (ms2_isotopic_dist$proba[1] + ms2_isotopic_dist$proba[2])

    added = rbind(added, data.frame(ms2_isotopic_ratio, ms2_isotopic_norm_ratio))
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


