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

variexpr finds all trios in a VCF, PED pair and let's the user specify an expression with indentifiers
of `kid`, `mom`, `dad` that is applied to each possible trio. samples that pass that filter have the id
of the kid added to the INFO field.

```
bpbio variexpr \
   --pass-only \ # output only variants that pass one of the filters.
   --vcf $vcf \
   --ped $ped \
   --load functions.js \ 
   --out-vcf annotated.bcf \
   --info "variant.call_rate > 0.9" \ # this filter is applied before the trio filters and can speed evaluation if it is stringent.
   --trio "denovo:kid.alts == 1 && mom.alts == 0 && dad.alts == 0 \
                   && kid.AB > 0.25 && kid.AB < 0.75 \
                   && (mom.AD[1] + dad.AD[1]) == 0 \
                   && kid.GQ >= 20 && mom.GQ >= 20 && dad.GQ >= 20 \
                   && kid.DP >= 12 && mom.DP >= 12 && dad.DP >= 12" \
   --trio "informative:kid.GQ > 20 && dad.GQ > 20 && mom.GQ > 20 && kid.alts == 1 && ((mom.alts == 1 && dad.alts == 0) || (mom.alts == 0 && dad.alts == 1))" \
   --trio "recessive:recessive_func(kid, mom, dad)"


Note that `variexpr` does not give direct access to the genotypes, instead exposing `alts` where 0 is homozygous reference, 1 is heterozygous, 2 is
homozygous alternate and -1 when the genotype is unknown. It is recommended to **decompose** a VCF before sending to `variexpr`

### Attributes

 + anything in the INFO is available as e.g. INFO.CSQ
 + if FORMAT.AB is not present, it is added so one can filter with kid.AB > 0.25 && kid.AB < 0.25
 + variant attributes are: `CHROM`, `POS`, `start`, `end`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`
 + calculated variant attributes include: `aaf`, `hwe_score`, `call_rate`, `num_hom_ref`, `num_het`, `num_hom_alt`, `num_unknown`

 + sample attributes (via `kid`, `mom`, `dad`) include in the FORMAT. available as e.g. kid.AD[1]
 + sample attributes from the ped for `affected`, `sex` are available as, e.g. kid.sex.
