context("simSeq")

test_that("compare to seq-gen", {

    skip_on_cran()
    # Hypothesis: 
    #   phangorn::simSeq and seq-gen generate alignments 
    #   with approximately the same number of mutations
    #
    # Method:
    #
    # - Create a phylogeny to work on
    # - Simulate alignments with seq-gen, count number of mutations
    # - Simulate alignments with phangorn::simSeq, count number of mutations
    # - Compare distributions of number of mutations
    # - Reject hypothesis if distributions are different

    set.seed(42)

    # Alignment length, in base pairs
    sequence_length <- 1000

    # Set the rigor of the experiment
    n_replicates <- 100

    # Create the simplest tree possible to simulate alignments on 
    tree <- ape::read.tree(text = "(t1:1,t2:1);")

    ############################################################################
    # seq-gen
    ############################################################################
    
    # Use a simple root sequence, note the uppercase adenine
    root_sequence <- rep("A", sequence_length)
    root_sequence_str <- paste0(root_sequence, collapse = "")
    
    # seq-gen 
    diffs_1 <- rep(NA, n_replicates)
    diffs_2 <- rep(NA, n_replicates)
    for (i in seq(1, n_replicates))
    {
    
      lines <- c(
        paste0("1 ", sequence_length),
        paste0("t0 ", root_sequence_str),
        "1",
        ape::write.tree(phy = tree)
      )
      filename <- tempfile()
      writeLines(text = lines, con = filename)
      # Put the root sequence at index 1
      out <- system2(command = "seq-gen", 
        args = c(paste0("-l", sequence_length), "-mHKY", "-k1"),
        stdin = filename,
        stdout = TRUE,
        stderr = FALSE
      )
      seq1 <- strsplit(out[2], split = " ")[[1]][9]
      seq2 <- strsplit(out[3], split = " ")[[1]][9]
      diffs_1[i] <- sum(strsplit(seq1, split = "")[[1]] != root_sequence)
      diffs_2[i] <- sum(strsplit(seq2, split = "")[[1]] != root_sequence)
    }
    seq_gen_diffs <- c(diffs_1, diffs_2)

    ############################################################################
    # phangorn::simSeq
    ############################################################################

    # Use a simple root sequence, note the lower-case adenine
    root_sequence <- rep("a", sequence_length)

    # Measure with `phangorn::simSeq`
    diffs_1 <- rep(NA, n_replicates)
    diffs_2 <- rep(NA, n_replicates)
    for (i in seq(1, n_replicates))
    {
      data <- simSeq(
        tree, 
        l = sequence_length, 
        rootseq = root_sequence, 
        type = "DNA", 
        rate = 1.0
      )
      diffs_1[i] <- sum(as.character(data)[1, ] != root_sequence)
      diffs_2[i] <- sum(as.character(data)[2, ] != root_sequence)
    }
    sim_seq_diffs <- c(diffs_1, diffs_2)

    ############################################################################
    # Compare
    ############################################################################
    
    # Use Kolmogorov-Smirnov test to determine if two distributions
    # are the same. Warning: p-value will be approximate in the presence of ties
    test <- suppressWarnings(
        ks.test(seq_gen_diffs, sim_seq_diffs)
    )
    # If p-value is less than 0.05, we assume the distributions to
    # be different. If the distributions are different, phangorn::simSeq
    # and seq-gen produce a different number of mutations. 
    p_value <- test$p.value
    
    # Expect distributions not to be similar
    expect_gt(p_value, 0.05)
    
    
    ############################################################################
    # test for site pattern
    ############################################################################
    
    tree <- read.tree(text="(t1:0.1,t2:0.1,t3:0.2);")
    dat <- allSitePattern(3)
    prob <- as.vector(exp(pml(tree, dat)$site))
    sitePattern <- function(x){
        attr(x, "index") <- NULL
        attr(x, "weight") <- rep(1, length(attr(x, "weight")))
        x <- as.data.frame(x, stringsAsFactors = FALSE)
        do.call("paste", c(x, sep = ""))  
    }
    ref <- sitePattern(dat)
    attr(dat, "weight") <- 1000 * prob
    fit <- pml(tree, dat)
    
    p_value1 <- numeric(100)
    p_value2 <- numeric(100)
    for(i in 1:100){
        x <- simSeq(tree)
        y <- sitePattern(x)
        obs <- numeric(length(ref))
        obs[match(y, ref)] <- attr(x, "weight")
        
        test <- suppressWarnings(
            chisq.test(obs, p=prob)
        )
        # If p-value is less than 0.05, we assume the distributions to
        # be different. 
        p_value1[i] <- test$p.value
    }
    expect_gt(mean(p_value1), 0.05)
    
    for(i in 1:100){
        x <- simSeq(fit)
        y <- sitePattern(x)
        obs <- numeric(length(ref))
        obs[match(y, ref)] <- attr(x, "weight")
        
        test <- suppressWarnings(
            chisq.test(obs, p=prob)
        )
        # If p-value is less than 0.05, we assume the distributions to
        # be different. 
        p_value2[i] <- test$p.value
    }
    expect_gt(mean(p_value2), 0.05)
    
})
