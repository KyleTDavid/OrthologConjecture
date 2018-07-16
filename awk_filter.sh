#!/bin/bash
#apply filters to masterdata matrix, for specifics please see text
awk -F '\t' '($4!~"Bilateria") && ($4!~"Caenorhabditis elegans") && ($4!~"Chordata") && ($4!~"Ciona") && ($4!~"Ciona intestinalis") && ($4!~"Ciona savignyi") && ($4!~"Drosophila melanogaster") && ($4!~"Ecdysozoa") && ($4!~"Opisthokonta") && ($4!~"Saccharomyces cerevisiae") && ($7 > 0) && ($7 <= 2) && ($8 >= 0.01) && ($8 <= 2) && ($12 > 0) && ($12 <= 2) && ($13 >= 0.01) && ($13 <= 2) && ($9 <= 10) && ($14 <= 10)' masterdata_unfiltered.txt > output.txt



