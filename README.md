[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.3936045.svg)](https://doi.org/10.5281/zenodo.3936045)

# longCombat: Longitudinal ComBat R Package

Longitudinal ComBat uses an empirical Bayes method to harmonize means and variances of the residuals across batches in a linear mixed effects model framework. Detailed methods are described in the manuscript: 

Beer JC, Tustison NJ, Cook PA, Davatzikos C, Sheline YI, Shinohara RT, Linn KA. (2020) Longitudinal ComBat: A method for harmonizing longitudinal multi-scanner imaging data. NeuroImage. In press. https://doi.org/10.1016/j.neuroimage.2020.117129.

Install package with: 
```{r, include=FALSE}
install.packages("devtools")
devtools::install_github("jcbeer/longCombat")
```

Note: longCombat currently will not run if tidyverse suite is loaded. This may be fixed in the future. For now, please run longCombat before loading tidyverse.

[Contact Joanne Beer](mailto:joanne.beer@pennmedicine.upenn.edu?subject=[GitHub]%20longCombat) with any questions, comments, or suggestions. Feedback is appreciated. 