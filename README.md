Riborex
=======

Riborex is a R package for identification of differential translation from Ribo-seq data.

DEPENDENCIES
------------
* DESeq2
* edgeR

INSTALLATION
------------
First, please make sure you have DESeq2 and edgeR installed.
Start R and enter:
```r
  ## try http:// if https:// URLs are not supported
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  biocLite("edgeR")
```
To install Riborex, download "riborex-1.2.3.tar.gz" from releases
at https://github.com/smithlabcode/riborex, start
a terminal and CD into the directory where you downloaded Riborex,
start R and enter
```r
  install.packages("riborex-1.2.3.tar.gz", repos=NULL, type="source")
```
Alternatively, you could also install `devtools` package and then install
`riborex` to get latest changes :

```r
  install.packages('devtools')
  library(devtools)
  options(unzip='internal')
  devtools::install_github('smithlabcode/riborex')
```


DOCUMENTATION
-------------
Please refer to vignettes/riborex.pdf for how to use riborex.

Contacts and bug reports
------------------------
Andrew D. Smith
andrewds@usc.edu

Wenzheng Li
wenzhenl@usc.edu

Weili Wang
weiliw@usc.edu

If you found a bug or mistake in this project, we would like to know about it.
Before you send us the bug report though, please check the following:

1. Are you using the latest version? The bug you found may already have been
   fixed.
2. Check that your input is in the correct format and you have selected the
   correct options.
3. Please reduce your input to the smallest possible size that still produces
   the bug; we will need your input data to reproduce the problem, and the
   smaller you can make it, the easier it will be.


Copyright and License Information
---------------------------------
Copyright (C) 2017 University of Southern California, Wenzheng Li, Weili Wang
and Andrew D. Smith

Authors: Wenzheng Li, Weili Wang, Philip J. Uren, Luiz OF Penalva, Andrew D. Smith

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see http://www.gnu.org/licenses/.
