---
title: "Exploratory analysis using aRchaic"
author: "Kushal K Dey"
date: "4/5/2017"
output: html_document
---

`aRchaic` provides various functions and tools to perform exploratory analysis of ancient DNA
sample using the mutational features recorded in its MFF file. In this section,
we shall focus on a few of these tools.

## Read length distribution

aRchaic provides a function `read_length_distribution ()` to track the density plot of the read
lengths for all reads, as well as reads containing C->T mutations and those containing C->T
mutations near the ends of the reads, latter two being potential representatives of
damaged reads. We demonstrate by applying this function on an MFF file


