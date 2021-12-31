#' 6VXX: Structure of the SARS-CoV-2 spike glycoprotein(closed state)
#'
#' The torsion angle dataset of SARS-CoV-2 spike glycopreotein.
#'
#' @format \code{data_6VXX} consists of following informations:
#' \describe{
#'   \item{\code{phi}}{main chain torsion angle for atoms C,N,CA,C.}
#'   \item{\code{psi}}{main chain torsion angle for atoms N,CA,C,N.}
#'   \item{\code{omega}}{main chain torsion angle for atoms CA,C,N,CA.}
#'   \item{\code{alpha}}{virtual torsion angle between consecutive C-alpha atoms.}
#'   \item{\code{chi1}}{side chain torsion angle for atoms N,CA,CB,*G.}
#'   \item{\code{chi2}}{side chain torsion angle for atoms CA,CB,*G,*D.}
#'   \item{\code{chi3}}{side chain torsion angle for atoms CB,*G,*D,*E.}
#'   \item{\code{chi4}}{side chain torsion angle for atoms *G,*D,*E,*Z.}
#'   \item{\code{chi5}}{side chain torsion angle for atoms *D,*E,*Z, NH1.}
#'   \item{\code{coords}}{numeric matrix of ‘justified’ coordinates.}
#'   \item{\code{tbl}}{a numeric matrix of psi, phi and chi torsion angles.}
#' }
#' @source This data can be downloaded in
#'   \url{https://www.rcsb.org/structure/6VXX}, or with using R package
#'   \code{bio3d}. Precisely, we use the code: \code{bio3d::torsion.pdb(bio3d::read.pdb("6vxx"))}
#' @references Walls, A. C., Park, Y. J., Tortorici, M. A., Wall, A., McGuire, A. T., & Veesler, D. (2020). Structure, function, and antigenicity of the SARS-CoV-2 spike glycoprotein. \emph{Cell}, 181(2), 281-292.
#'   Retrived from \url{https://www.wwpdb.org/pdb?id=pdb_00006vxx}
#' @seealso Description of the angluar information is from the 'value'
#'   part of \code{torsion.pdb} in the package \code{bio3d}.
"data_6VXX"

#' SARS-CoV-2: chain B of Structure of the SARS-CoV-2 spike glycoprotein(closed state)
#'
#' The torsion angle dataset of the chain B of SARS-CoV-2 spike glycopreotein. This data is originally
#' from first two main torsion angles of \code{data_6VXX}.
#'
#' This data is obtained with following codes:
#'
#' @format This data.frame contains the following columns:
#' \describe{
#'   \item{\code{phi}}{main chain torsion angle for atoms C,N,CA,C.}
#'   \item{\code{psi}}{main chain torsion angle for atoms N,CA,C,N.}
#' }
#' @source This data can be downloaded in
#'   \url{https://www.rcsb.org/structure/6VXX}, or with using R package
#'   \code{bio3d}. To see the precise extracting code, visit \url{https://github.com/sungkyujung/ClusTorus/tree/master/data-raw}
#' @references Walls, A. C., Park, Y. J., Tortorici, M. A., Wall, A., McGuire, A. T., & Veesler, D. (2020). Structure, function, and antigenicity of the SARS-CoV-2 spike glycoprotein. \emph{Cell}, 181(2), 281-292.
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
#'   \item{\code{phi}}{main chain torsion angle for atoms C,N,CA,C.}
#'   \item{\code{psi}}{main chain torsion angle for atoms N,CA,C,N.}
#'   \item{\code{chi1}}{side chain torsion angle for atoms N,CA,CB,*G.}
#'   \item{\code{chi2}}{side chain torsion angle for atoms CA,CB,*G,*D.}
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
#'   The data extracting method is from Harder, T., Boomsma, W., Paluszewski, M., Frellsen, J., Johansson, K. E., & Hamelryck, T. (2010). Beyond rotamers: a generative, probabilistic model of side chains in proteins. \emph{BMC bioinformatics}, 11(1), 1-13.
#'
#'   Mardia, K. V., Kent, J. T., Zhang, Z., Taylor, C. C., & Hamelryck, T. (2012). Mixtures of concentrated multivariate sine distributions with applications to bioinformatics. \emph{Journal of Applied Statistics}, 39(11), 2475-2492.
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
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
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
#' @references Jung, S., Park, K., & Kim, B. (2021). Clustering on the torus by conformal prediction. \emph{The Annals of Applied Statistics}, 15(4), 1583-1603.
"toydata2"
