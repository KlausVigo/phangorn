---
layout: page
title: Using binary or discrete data with phangorn
description: Using binary or discrete data with phangorn
---



If data are not sequence alignment an `phyDat` object then there are generic functions `as.phyDat()` in phangorn to transform matrices and data.frames into phyDat objects.
For example you can read in your data with `read.table()` or `read.csv()`, but you might need to transpose your data. For matrices `as.phyDat()` assumes that the entries each row belongs to one individual (taxa), but for data.frame each column. For binary data you can transform these with a command like (depending how you coded them):

    as.phyDat(data, type="USER", levels = c(0, 1))
    as.phyDat(data, type="USER", levels = c(TRUE, FALSE))


