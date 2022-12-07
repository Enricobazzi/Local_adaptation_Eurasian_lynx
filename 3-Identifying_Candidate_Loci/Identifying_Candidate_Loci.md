---
title: "Identifying Candidate Loci"
author: "Enrico"
date: "2022-10-17"
output: html_document
---


```
cd /home/ebazzicalupo/Local_adaptation_Eurasian_lynx

for var in bio9 bio2 bio15 bio16 mean_snow_days
 do
  echo "analyzing association with ${var}"
  screen -dmS aux_${var}  sh -c "3-Identifying_Candidate_Loci/code/baypass_aux_v2.sh ${var}; exec /bin/bash"
done

for var in bio9 bio2 bio15 bio16 mean_snow_days
 do
  echo "generating genwin windows from baypass results of ${var}"
  Rscript 3-Identifying_Candidate_Loci/code/run_genwin.R ${var}
done
```

```
for i in 1 2 4 6 7 8
 do
  echo "running RDA with line $i of formula file"
  Rscript 3-Identifying_Candidate_Loci/code/run_rda_exploration.R $i
done
```

```

```

