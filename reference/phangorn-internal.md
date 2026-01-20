# Internal phangorn Functions

Internal phangorn functions.

## Usage

``` r
threshStateC(x, thresholds)

base_freq(obj)

assert_phyDat(x, contains_label = !is.null(label), label = NULL)

assert_phylo(x, has_edge_length = FALSE, is_rooted = FALSE,
  is_ultrametric = FALSE, is_binary = FALSE, unique_tiplabel = TRUE)

assert_multiPhylo(x, has_edge_length = FALSE, is_rooted = FALSE,
  is_ultrametric = FALSE, is_binary = FALSE)

assert_treeish(x, null.ok = FALSE)

assert_phyDat(x, contains_label = !is.null(label), label = NULL)

assert_pml(x)

clean_phylo(x, ..., unroot = FALSE, collapse.singles = FALSE,
  reorder = FALSE, order = "postorder", multi2di = FALSE,
  di2multi = FALSE, tol = 1e-07, compress = FALSE)

candidate_tree(x, method = c("unrooted", "ultrametric", "tipdated"),
  eps = 1e-08, tip.dates = NULL, ...)

hash(x, ...)

hash_topo(x, rooted = FALSE, ...)

binomial_mk(k = 4)

map_duplicates(x, dist = length(x) < 500, ...)

edQt(Q = c(1, 1, 1, 1, 1, 1), bf = c(0.25, 0.25, 0.25, 0.25))

pml.free()

pml.init(data, k = 1L)

coords(obj, dim = "3D")

pmlPen(object, lambda, ...)

relabel(x, ref)

edge_length_hclust(x, dm, method = "average")
```
