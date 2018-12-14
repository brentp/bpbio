# bpbio: a single binary with a collection of sub-commands useful for genomics

+ plot-sv-vcf: create a plot that shows the number of large and small variants separated by DEL/DUP/BND/INV for a multi-sample SV VCF. (idea from @ernfrid)
+ variexpr: simple expression on variants for great good


## Libraries
+ countstats: get median and percentiles for integers in constant space and little time.
+ pedfile: ped (fam) file parsing and matching to VCF.
+ duko: ergonomic embedded javascript via [duktape](https://duktape.org/) and [duktape-nim](https://github.com/manguluka/duktape-nim)
+ variexpr: simple expressions for filtering and labeling VCF records e.g.: denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && kid.DP > 10


## Abandoned 

+ homsv: look for depth changes in self-chains or homologous regions
+ homsv-merge: merge output from homsv

## Variexpr

```
bpbio variexpr \
   --pass-only --vcf $vcf \
   --ped $ped \
   --load functions.js \
   --out-vcf annotated.bcf \
   --trio "denovo:INFO.AC == 1 && kid.alts == 1 && mom.alts == 0 && dad.alts == 0 && (mom.AD[1] + dad.AD[1]) == 0 \
                   && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 && kid.DP >= 12 && mom.DP >= 12 && dad.DP >= 12" \
   --trio "informative:kid.GQ > 20 && dad.GQ > 20 && mom.GQ > 20 && kid.alts == 1 && ((mom.alts == 1 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 1))" \
   --trio "recessive:recessive_func(kid, mom, dad)"
