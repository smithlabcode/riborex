# Riborex

Riborex is a R package for identification of differential translation from Ribo-seq data.

--------------------------------------------------------

[Online Paper](https://academic.oup.com/bioinformatics/article/33/11/1735/2964727) | [PDF](https://watermark.silverchair.com/btx047.pdf?token=AQECAHi208BE49Ooan9kkhW_Ercy7Dm3ZL_9Cf3qfKAc485ysgAAAr0wggK5BgkqhkiG9w0BBwagggKqMIICpgIBADCCAp8GCSqGSIb3DQEHATAeBglghkgBZQMEAS4wEQQMEeGzK92KGEWtk0FoAgEQgIICcFUmPpYJo1tGujsBoP4lYse6EnU5MCy_Hp5z2ae7rC9jSUgL4I_UscKpEadF2YIwltYzPiN63ZHToP6ec-l8MPEDHwrFISgfjFkYiz7CdWwTLFo-1zKx449YprMc7gPtXGUjF4koXDC0VUKEAz5Df6litXeErwErQfN1MvTfapzv6e2uocdeaIb_gVHP8_SOgzivRHJ5iFDunIjS_43s_G5N42L4lIOAaEl_5gKQ7MmBoH04Nzd8yc8GX9Q6D5gBri-Lqsu3HM5x4_EuGRrAp8dP40rtwCBF7yoPj_VHhirTgN1XP82JZJO7QWzqnsZXJ6n9w4lOTgxpY24DONTXGAVAAf-CwkuIK58cvqvrulZwOwj9cikYExj5sXpYY9J4u1y-wVufDIFBLa912qaVRJV9h7k3nGihq7MJkpNKLwTzqj35VbtmGDgqoDJMQxjs9K5WdSCW4Ecc6uPFQk9XwC7PzYAiJIRn7mjVUjyciBKe9XRvcpJbK4C3nl__XDRhM9sFmNlA7mwPhIePhOhzhLWACUSCltLBPiRNDNv34VB_Rj_--SPxG5gnRmfqGew7GeDYUo_W-rjvCC0PeBASozZUe3h6sT4s7eIXwWLbzuJ_ujWcXM7jLaqCjJgl-SGUcEYUiWUvXK9YfGNshpgAf99cCpu4mNhGxPL672URz1_RhY2kS2zHbG84MOgdd-SqvarTGYugArsUN6AlMXvHTR5ZIWr4PE47FjOYPCn4NS9QV-42dpj1a7TuJSsDlbDdoHdXcxXTTFhO_wfKjVy_OCXR1jrqxSHTURi3NR6yF0W4mO-oCuTM7DS4RtI1LvnhDA) | [Supplementary File](https://oup.silverchair-cdn.com/oup/backfile/Content_public/Journal/bioinformatics/33/11/10.1093_bioinformatics_btx047/3/btx047_supp.pdf?Expires=1620962473&Signature=h3eYZdiu4ddVx5L-NXUsO123Ej50pLrumpFn74i1FjsaS7wbQNqrJqMZg4G9XKryOWUJXNBDLxeTnYBXOKPHW5O1QFuDS6k8r8gvjwSMpEtOMdR0wP1tsnhKKi~Oi8sBIPEdHalsQcJCZMZXkNuKpav7u~dgDZ2pUkhiaGoOyLomnBHfafjM79-Lprnk4NLFO3KzOOZ0-IPtSvVDYczqEo8HascrsO1B3Rbs-1aQdQv8~hoTAfV7Y8hXOA3XO2WSO1Kp16qrgmqUWCtqJUfd6QZPlwy~JbZ4OS9HPjTNyM4e0JVWLuWNsjN4~aQixjZ9S9JEzEmmnBHGndxaH6Xs~g__&Key-Pair-Id=APKAIE5G5CRDK6RD3PGA)

--------------------------------------------------------

## DEPENDENCIES
* DESeq2
* edgeR
* fdrtool
--------------------------------------------------------
## INSTALLATION
We strongly recommend that you install Riborex via conda:
```bash
   conda install -c bioconda r-riborex
```
To install locally, please make sure you have DESeq2 and edgeR installed.
Then start R and enter:
```r
  ## try http:// if https:// URLs are not supported
  source("https://bioconductor.org/biocLite.R")
  biocLite("DESeq2")
  biocLite("edgeR")
  biocLite('fdrtool')
```
To install Riborex, download the latest version
"riborex-x.x.x.tar.gz" from releases
at https://github.com/smithlabcode/riborex, start
a terminal and CD into the directory where you downloaded Riborex,
start R and enter
```r
  install.packages("riborex-x.x.x.tar.gz", repos=NULL, type="source")
```
Alternatively, you could also install `devtools` package and then install
`riborex` to get latest changes :

```r
  install.packages('devtools')
  library(devtools)
  options(unzip='internal')
  devtools::install_github('smithlabcode/riborex')
```
--------------------------------------------------------

## DOCUMENTATION
Please refer to vignettes/riborex.pdf for how to use riborex.

## Contacts and bug reports
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

--------------------------------------------------------
## Copyright and License Information
Copyright (C) 2017-2020 University of Southern California, Wenzheng Li, Weili Wang
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
