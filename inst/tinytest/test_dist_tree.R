data(Laurasiatherian)
dm <- dist.ml(Laurasiatherian)

# check nnls functions
    tree_nj <- NJ(dm)
    tree_unj <- UNJ(dm)
    tree_nnls_unj <- nnls.phylo(tree_unj, dm)
    tree_nnls_nj <- nnls.phylo(tree_nj, dm)

    tree_upgma <- upgma(dm)
    tree_wpgma <- wpgma(dm)
    tree_nnls_upgma <- nnls.phylo(tree_upgma, dm, rooted=TRUE)
    tree_nnls_wpgma <- nnls.phylo(tree_wpgma, dm, rooted=TRUE)

    expect_equal(tree_upgma, tree_nnls_upgma)
    expect_false(all.equal(tree_wpgma, tree_nnls_wpgma))
    expect_equal(tree_unj, tree_nnls_unj)
    expect_false(all.equal(tree_nj, tree_nnls_nj))
