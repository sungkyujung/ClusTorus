#' 6VXX: Structure of the SARS-CoV-2 spike glycoprotein(closed state)
#'
#' The torsion angle dataset of the SARS-CoV-2 spike glycopreotein.
#'
#' @format This list contains the following components:
#' \describe{
#'   \item{\code{psi}}{protein backbone chain angle
#'     for atoms C, N, C_alpha and C, in arc-degree.}
#'   \item{\code{phi}}{protein backbone chain angle
#'     for atoms N, C_alpha, C and N, in arc-degree.}
#'   \item{\code{omega}}{protein backbone chain angle
#'     for atoms C_alpha, C, N and C-alpha, in arc-degree.}
#'   \item{\code{chi1}}{side chain torsion angle
#'     for atoms N, C_alpha, C_beta and *G, in arc-degree.}
#'   \item{\code{chi2}}{side chain torsion angle
#'     for atoms C_alpha, C_beta, *G and *D, in arc-degree.}
#'   \item{\code{chi3}}{side chain torsion angle
#'     for atomsC_beta, *G, *D and *E, in arc-degree.}
#'   \item{\code{chi4}}{side chain torsion angle
#'     for atoms *G, *D, *E and *Z, in arc-degree.}
#'   \item{\code{chi5}}{side chain torsion angle
#'     for atoms *D, *E, *Z and NH1, in arc-degree.}
#'   \item{\code{alpha}}{virtual torsion angle between
#'     consecutive C_alpha atoms.}
#'   \item{\code{coords}}{numeric matrix of 'justified' coordinates.}
#'   \item{\code{tbl}}{a numeric matrix of psi, phi, and chi torsion
#'     angles.}
#' }
#' @source This data can be downloaded in
#'   \url{https://www.rcsb.org/structure/6VXX}, or with using R package
#'   \code{bio3d}.
#' @references Walls, A.C., et al. (2020), "Structure of the SARS-CoV-2
#'   spike glycoprotein (closed state)" Cell 181: 281, DOI:10.2210/pdb6vxx/pdb.
#'   Retrived from \url{https://www.wwpdb.org/pdb?id=pdb_00006vxx}
#' @seealso Description of the angluar information is from the 'value'
#'   part of \code{torsion.pdb} in the package \code{bio3d}.
"SARS_CoV_2"

#' ILE: Structure of the Isoleucine
#'
#' An isomer of leucine, essential branched-chain aliphatic amino acid found in many proteins.
#'
#' @format This list contains the following components:
#' \describe{
#'   \item{\code{tbl}}{a numeric matrix of psi, phi, and chi torsion
#'     angles.}
#' }
#' @details
#'   ILE data is generated with collection of different pdb files. To select adequate protein
#'   data, we use PISCES server. (the method is introduced in articles of references.)
#'   To select high-quality protein data, we use several benchmarks:
#'   resolution : 1.6A(angstrom) or better,
#'   R-factor : 0.22 or better,
#'   Sequence percentage identity: <= 25%.
#'   Then, we select ILE only angular data for each protein data. To see the detail code, visit
#'   \url{https://github.com/sungkyujung/ClusTorus}
#' @source This data is extracted from PISCES server \url{http://dunbrack.fccc.edu/pisces/}
#' @references Data description is from \url{https://www.rcsb.org/ligand/ILE}.
#'
#'   The data extracting method is from Harder, T., Boomsma, W., Paluszewski, M. et al.(2010)
#'   "Beyond rotamers: a generative, probabilistic model of side chains in proteins". BMC Bioinformatics 11, 306
#'   \doi{https://doi.org/10.1186/1471-2105-11-306} and
#'
#'   Kanti V. Mardia , John T. Kent , Zhengzheng Zhang , Charles C. Taylor & Thomas Hamelryck (2012)
#'   "Mixtures of concentrated multivariate sine distributions with applications to bioinformatics",
#'   Journal of Applied Statistics, 39:11, 2475-2492, DOI: 10.1080/02664763.2012.719221
#'
#' @seealso Description of the angluar information is from the 'value'
#'   part of \code{torsion.pdb} in the package \code{bio3d}.
"ILE"

#' toydata1: Labelled Data for 5 Clusters
#'
#' Artificially generated data on the 2 dimensional torus
#'
#' @format This \code{data.frame} contains the following components:
#' \describe{
#'   \item{\code{phi}}{column for the first angle}
#'   \item{\code{psi}}{column for the second angle}
#'   \item{\code{label}}{column for the clustering membership}
#' }
#' @details
#'   toydata1 is an artificial data generated from a mixture of 5 clusters,
#'   where three clusters are sampled from bivariate normal distributions
#'   and the other two are each sampled from the uniform distribution on a rectangle.
#'
#' @references This simulation data is from S. Jung, K. Park, B. Kim (2021)
#'   "Clustering on the torus by conformal prediction". Annals of Applied Statistics
"toydata1"

#' toydata2: Labelled Data for 3 Clusters
#'
#' Artificially generated data on the 2 dimensional torus
#'
#' @format This \code{data.frame} contains the following components:
#' \describe{
#'   \item{\code{phi}}{column for the first angle}
#'   \item{\code{psi}}{column for the second angle}
#'   \item{\code{label}}{column for the clustering membership}
#' }
#' @details
#'   toydata2 is an artificial data generated from a mixture of 3 clusters,
#'   where the first cluster is sampled from a spherical normal distribution,
#'   the second cluster is from the uniform distribution on a large “L”-shaped region,
#'   and the third cluster of size 50 is sampled from the uniform distribution on the
#'   entire 2-dimensional torus.
#'
#' @references This simulation data is from S. Jung, K. Park, B. Kim (2021)
#'   "Clustering on the torus by conformal prediction". Annals of Applied Statistics
"toydata2"
