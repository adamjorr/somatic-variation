library(vcfR)
library(dplyr)
library(stringr)
library(tidyr)

denovos = read.vcfR("deduped.vcf.gz")
phased = read.vcfR("../phasing/phased.vcf.gz")

db = vcfR2tidy(denovos)
pos = db$fix %>% select(CHROM,POS) %>% split(., .$CHROM)

a = extract_gt_tidy(phased, "PS")

dbPS = tibble(CHROM = getCHROM(phased), POS = getPOS(phased), PS = a$gt_PS) %>% filter(!is.na(PS))
rlePS = dbPS %>% pull(PS) %>% split(dbPS$CHROM) %>% lapply(rle)
cumPS = lapply(rlePS,function(x) cumsum(x$lengths))

hpw = lapply(setNames(names(pos),names(pos)), function(n) {
    dnmpos = pos[[n]] %>% pull(POS)
    pspos = dbPS %>% filter(CHROM==n) %>% pull(POS)
    hiz = cumPS[[n]]
    loz = (hiz-rlePS[[n]]$lengths+1)

    p = sapply(dnmpos, function(y) {
        lo = max(which(y > pspos))
        hi = min(which(y < pspos))
        lo = min(which(lo < hiz))
        hi = min(which(hi < hiz))
        if(lo != hi) {
            return(0)
        }
        hi = hiz[lo]
        lo = loz[lo]
        pspos[hi]-pspos[lo]+1
    })
    cbind(pos[[n]],HPLEN=p)
})

x = do.call(rbind,hpw)
write.table(x,file="hplen.tsv",row.names=F,quote=F,sep="\t")
