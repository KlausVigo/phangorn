#' srr_stats
#'
#' All of the following standards initially have `@srrstatsTODO` tags.
#' These may be moved at any time to any other locations in your code.
#' Once addressed, please modify the tag from `@srrstatsTODO` to `@srrstats`,
#' or `@srrstatsNA`, ensuring that references to every one of the following
#' standards remain somewhere within your code.
#' (These comments may be deleted at any time.)
#'
#' @srrstatsVerbose TRUE
#'
#' @noRd
NULL

#' NA_standards
#'
#' Any non-applicable standards can have their tags changed from `@srrstatsTODO`
#' to `@srrstatsNA`, and placed together in this block, along with explanations
#' for why each of these standards have been deemed not applicable.
#' (These comments may also be deleted at any time.)
#'
#' @srrstatsNA {G1.4a} *All internal (non-exported) functions should also be documented in standard [`roxygen2`](https://roxygen2.r-lib.org/) format, along with a final `@noRd` tag to suppress automatic generation of `.Rd` files.*
#' This seems disproportionate work and might lead to poor programming style
#'
#' The following 2 standards will be in an appendix with the paper as comparison
#' are too time consuming for CRAN and require often different software installed.
#' @srrstatsNA {G1.5} *Software should include all code necessary to reproduce results which form the basis of performance claims made in associated publications.*
#' See also comments below. The package is pretty old so I wonder if those claims
#' should be more relevant to packages offering an alternative implementation to phangorn?
#' @srrstatsNA {G1.6} *Software should include code necessary to compare performance claims with alternative implementations in other R packages.*
#' @srrstatsNA {G2.2} *Appropriately prohibit or restrict submission of multivariate input to parameters expected to be univariate.*
#' @srrstatsNA {G2.5} *Where inputs are expected to be of `factor` type, secondary documentation should explicitly state whether these should be `ordered` or not, and those inputs should provide appropriate error or other routines to ensure inputs follow these expectations.*
#' @srrstatsNA {G2.4a} *explicit conversion to `integer` via `as.integer()`*
#' @srrstatsNA {G2.4b} *explicit conversion to continuous via `as.numeric()`*
#' @srrstatsNA {G2.4c} *explicit conversion to character via `as.character()` (and not `paste` or `paste0`)*
#' @srrstatsNA {G2.4d} *explicit conversion to factor via `as.factor()`*
#' @srrstatsNA {G2.6} *Software which accepts one-dimensional input should ensure values are appropriately pre-processed regardless of class structures.*
#' @srrstatsNA {G2.10} *Software should ensure that extraction or filtering of single columns from tabular inputs should not presume any particular default behaviour, and should ensure all column-extraction operations behave consistently regardless of the class of tabular data used as input.*
#' @srrstatsNA {G2.11} *Software should ensure that `data.frame`-like tabular objects which have columns which do not themselves have standard class attributes (typically, `vector`) are appropriately processed, and do not error without reason. This behaviour should be tested. Again, columns created by the [`units` package](https://github.com/r-quantities/units/) provide a good test case.*
#' @srrstatsNA {G2.12} *Software should ensure that `data.frame`-like tabular objects which have list columns should ensure that those columns are appropriately pre-processed either through being removed, converted to equivalent vector columns where appropriate, or some other appropriate treatment such as an informative error. This behaviour should be tested.*
#' @srrstatsNA {G2.14a} *error on missing data*
#' @srrstatsNA {G2.14c} *replace missing data with appropriately imputed values*
#' @srrstatsNA {G3.1} *Statistical software which relies on covariance calculations should enable users to choose between different algorithms for calculating covariances, and should not rely solely on covariances from the `stats::cov` function.*
#' @srrstatsNA {G3.1a} *The ability to use arbitrarily specified covariance methods should be documented (typically in examples or vignettes).*
#' @srrstatsNA {G5.6b} *Parameter recovery tests should be run with multiple random seeds when either data simulation or the algorithm contains a random component. (When long-running, such tests may be part of an extended, rather than regular, test suite; see G5.10-4.12, below).*
#' @srrstatsNA {G5.11} *Where extended tests require large data sets or other assets, these should be provided for downloading and fetched as part of the testing workflow.*
#' @srrstatsNA {G5.11a} *When any downloads of additional data necessary for extended tests fail, the tests themselves should not fail, rather be skipped and implicitly succeed with an appropriate diagnostic message.*
#' @srrstatsNA {EA2.0} *EDA Software which accepts standard tabular data and implements or relies upon extensive table filter and join operations should utilise an **index column** system*
#' @srrstatsNA {EA2.1} *All values in an index column must be unique, and this uniqueness should be affirmed as a pre-processing step for all input data.*
#' @srrstatsNA {EA2.2} *Index columns should be explicitly identified, either:*
#' @srrstatsNA {EA2.2a} *by using an appropriate class system, or*
#' @srrstatsNA {EA2.2b} *through setting an `attribute` on a table, `x`, of `attr(x, "index") <- <index_col_name>`.*
#' @srrstatsNA {EA2.3} *Table join operations should not be based on any assumed variable or column names*
#' @srrstatsNA {EA2.4} *Use and demand an explicit class system for such input (for example, via the [`DM` package](https://github.com/krlmlr/dm)).*
#' @srrstatsNA {EA2.5} *Ensure all individual tables follow the above standards for Index Columns*
#' @srrstatsNA {EA5.3} *Column-based summary statistics should always indicate the `storage.mode`, `class`, or equivalent defining attribute of each column.*
#' @srrstatsNA {EA5.5} *All visualisations should include units on all axes where such are specified or otherwise obtainable from input data or other routines.*
#' @srrstatsNA {EA5.6} *Any packages which internally bundle libraries used for dynamic visualization and which are also bundled in other, pre-existing R packages, should explain the necessity and advantage of re-bundling that library.*
#' @srrstatsNA {EA6.0b} *Dimensions of tabular objects*
#' @srrstatsNA {EA6.0c} *Column names (or equivalent) of tabular objects*
#' @srrstatsNA {EA6.0d} *Classes or types of all columns contained within `data.frame`-type tabular objects *
#' @srrstatsNA {EA6.0e} *Values of single-valued objects; for `numeric` values either using `testthat::expect_equal()` or equivalent with a defined value for the `tolerance` parameter, or using `round(..., digits = x)` with some defined value of `x` prior to testing equality.*
#' @srrstatsNA {UL1.2} *Unsupervised learning which uses row or column names to label output objects should assert that input data have non-default row or column names, and issue an informative message when these are not provided.*
#' @srrstatsNA {UL2.1} *Unsupervised Learning Software should document any transformations applied to input data, for example conversion of label-values to `factor`, and should provide ways to explicitly avoid any default transformations (with error or warning conditions where appropriate).*
#' @srrstatsNA {UL2.2} *Unsupervised Learning Software which accepts missing values in input data should implement explicit parameters controlling the processing of missing values, ideally distinguishing `NA` or `NaN` values from `Inf` values.*
#' @srrstatsNA {UL2.3} *Unsupervised Learning Software should implement pre-processing routines to identify whether aspects of input data are perfectly collinear.*
#' @srrstatsNA {UL3.0} *Algorithms which apply sequential labels to input data (such as clustering or partitioning algorithms) should ensure that the sequence follows decreasing group sizes (so labels of "1", "a", or "A" describe the largest group, "2", "b", or "B" the second largest, and so on.)*
#' @srrstatsNA {UL3.1} *Dimensionality reduction or equivalent algorithms which label dimensions should ensure that that sequences of labels follows decreasing "importance" (for example, eigenvalues or variance contributions).*
#' @srrstatsNA {UL3.2} *Unsupervised Learning Software for which input data does not generally include labels (such as `array`-like data with no row names) should provide an additional parameter to enable cases to be labelled.*
#' @srrstatsNA {UL3.4} *Objects returned from Unsupervised Learning Software which labels, categorise, or partitions data into discrete groups should include, or provide immediate access to, quantitative information on intra-group variances or equivalent, as well as on inter-group relationships where applicable.*
#' @srrstatsNA {UL7.5} *Batch processing routines should be explicitly tested, commonly via extended tests (see G4.10–G4.12).*
#' @srrstatsNA {UL7.5a} * Tests of batch processing routines should demonstrate that equivalent results are obtained from direct (non-batch) processing.*
#' @noRd
NULL


#' General comments
#'
#' @srrstats {G1.1}
#' The package is on CRAN since 2008-03-25. So in many cases the phangorn will
#' have the first implementation within R of an algorithm,
#' which has previously been implemented in other languages or contexts,
#' e.g. for phylogenetic inference for maximum parsimony and maximum likelihood.
#' Plot methods like consensusNetworks, densiTree like plots, split networks.
#'
#' @srrstats {G1.2} CONTRIBUTING.md contains some information about the life cycle.
#' As phangorn is already over 10 years on CRAN many functions are mature.
#' Some newer functions or less used function might be more experimental.
#' This is usually stated in the manual.
#' @srrstats {G1.3} We try to be describe the terminology as clearly and unambiguously as possible.
#' *All statistical terminology should be clarified and unambiguously defined.*
#' @srrstats {G1.4}
#' All Rd.pages are generated with roxygen2.
#' *Software should use [`roxygen2`](https://roxygen2.r-lib.org/) to document all functions.*
#'
#' The package has some 172 and counting functions, many have some form of checking input arguments.
#' I started to do this more systematic using the checkmate package for all user facing functions.
#' Additionally there are functions in assert.R to test common objects in phangorn.
#' @srrstats {G2.0} *Implement assertions on lengths of inputs, particularly through asserting that inputs expected to be single- or multi-valued are indeed so.*
#' @srrstats {G2.0a} *Provide explicit secondary documentation of any expectations on lengths of inputs*
#' @srrstats {G2.1} *Implement assertions on types of inputs (see the initial point on nomenclature above).*
#' @srrstats {G2.1a} *Provide explicit secondary documentation of expectations on data types of all vector inputs.*
#'
#' phyDat stores categorial data similar to a factor (one assumes that the whole
#' matrix share the same states.
#' When transforming matrices to an phyDat object missing values are deleted and
#' a warning is issued
#'
#' @srrstats {G2.13, G2.14, G2.14b}
#' *Statistical Software should implement appropriate checks for missing data as part of initial pre-processing prior to passing data to analytic algorithms.*
#' *Where possible, all functions should provide options for users to specify how to handle missing (`NA`) data, with options minimally including:*
#'
#' General comments to EDA:
#'
#' The target audience are on the one hand (evolutionary) biologists wanting to
#' infer phylogenies using MP, ML, distance methods.
#' Several methods for tipdated or ultrametric tree exist.
#' Tree comparison and some of the distance based methods (e.g. UPGMA with NNI)
#' provides improve and generalisation on hierarchical clustering.
#' Additionally phangorn can be used in teaching as whole workflow can be done
#' from EDA, inference and interpreting results are possible in one environment.
#' @srrstats {EA1.0} *Identify one or more target audiences for whom the software is intended*
#' The input for MP, ML methods are multiple sequence alignments. There exist
#' several possibilities to read different file formats and transform msa from different formats.
#' ?read.phyDat methods(as.phyDat) connecting to bioconductor etc.
#' For distance based methods distances (objects of class dist).
#' @srrstats {EA1.1} *Identify the kinds of data the software is capable of analysing (see *Kinds of Data* below).*
#' @srrstats {EA1.2} *Identify the kinds of questions the software is intended to help explore.*
#' @srrstats {EA1.3} *Identify the kinds of data each function is intended to accept as input*
#' @srrstats {EA2.6} we use inherits, so phylo or multiPhylo object might have
#' additional attributes
#' @srrstats {EA4.1, EA4.2, EA5.2, EA5.4} densiTree and plotSeqLogo plotBS, round values,
#' use pretty for axes. (examples are densiTree.R. plotAnc.R, plotBS.R)
#'
#' @srrstats {EA5.0} We follow best practices to be as accessible as possible.
#' @srrstats {EA5.0a} We try our best, this is with trees which can differ in
#' number of labels quite challenging.
#' @srrstats {EA5.0b} in the image.phyDat we use color schemes which are
#' commonly used in bioinformatics, otherwise we often use color blind friendly
#' palettes.
#' @srrstats {EA5.1} Graphical parameter are often from called via '...' argument,
#' relying on R standards, otherwise these explicitly stated.
#'
#' General comments to UL (or not yet put into different files):
#'
#' @srrstats {UL2.0} That's why there are the EDA functions available.
#' One major point is to hierarchical clustering is that in phylogenetics the
#' assumption that trees / distances are ultrametric is often relaxed
#' (additive distances).
#' @srrstats {UL3.3} Ancestral reconstruction is an example where ancestral
#' states can fast extracted from a model.
#' @srrstats {UL4.0, UL4.1, UL4.2, UL4.3, UL4.3a, UL4.4} pml object has print
#' methods, contains a call object, a print method (and several more generics).
#' @srrstats {UL6.0, UL6.1, UL6.2} There are several generic plot methods
#' plot.pml(), plot.networx().
#' Labels (e.g.) bootstrap support are placed for readability. The plot functions
#' mainly build on the ape package.
#'
NULL

