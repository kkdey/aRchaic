# aRchaic

a R/python software for exploration, clustering, visualization and classification of DNA damage patterns 

## Authors

[Kushal K Dey*](http://kkdey.github.io/), [Hussein al Asadi
*](https://halasadi.wordpress.com/), [John Novembre](http://jnpopgen.org/), [Matthew Stephens](http://stephenslab.uchicago.edu/)


## Installation

The user can install the `aRchaic` package in R from Github via devtools.

```
devtools::install_github("kkdey/aRchaic")
```

To load the package in R

```
library(aRchaic)
```

## Methods

We propose a model based approach to cluster ancient (UG/non UDG) and modern samples based on known DNA damage patterns like type of mismatch, the strand breaks and position of mismatch on the read. Our model is based on a fast and more efficient version of the Grade of Membership model and we present a novel way to interpret and visualize the clusters. Additionally we also provide functions to visualize DNA damage patterns for a single BAM file, methods to classify each read in the BAM file as modern or ancient, along with model based classification techniques for moderns and ancient samples.


Check our [Project Webpage](https://kkdey.github.io/aRchaic/)

## Contact

For any inquiries or questions related to the package, please open an issue in this repository. You can also contact us at [kkdey@uchicago.edu](kkdey@uchicago.edu) or [halasadi@gmail.com](halasadi@uchicago.edu)


Also users are welcome to contribute to the package by submitting pull request. 

## Acknowledgements

The authors would like to acknowledge Anna Di Rienzo, Choongwon Jeong, Anna Gosling, John Lindo, David Witonsky, Joseph Marcus, John Blischak and members of Stephens Lab and Novembre Lab for helpful discussions.



