#!/usr/bin/Rscript --vanilla
library(RColorBrewer)
library(stringr)
library(data.table)
library(vcfR)

scriptname = commandArgs() %>% str_subset("^--file=") %>% str_replace("^--file=","")
setwd(dirname(scriptname))
pdfname = basename(scriptname) %>% str_replace("\\.[^.]*$", ".pdf")

denovos = read.vcfR("vcf/goodsites.vcf")
db = vcfR2tidy(denovos)

# analysis indicates that all het->hom mutations occur at the root
# therefore we will swap them to be hom->het
v = db$fix$DNT %>% str_detect("^(.)/\\1") %>% `!` %>% which

s = db$fix$DNT %>% str_split("->",simplify=TRUE)

x = apply(s,1,function(g){
    g = str_split(g,"/",simplify=TRUE)
    if(g[1,1] != g[1,2]) {
        g = g[c(2,1),]
    }
    o = g[1,] != g[2,]
    if(sum(o) != 1) {
        return(NULL)
    }
    g[,o]
})
x = t(x)

tr1 = function(x) {
    switch(x, "A" = "T", "C" = "G", "G" = "C", "T" = "A")
}

tr2 = function(x) {
    switch(x, "A" = "T:A", "C" = "G:C", "G" = "C:G", "T" = "A:T")
}


o = str_detect(x[,1], "^[AG]")
x[o,1] = str_replace_all(x[o,1], "[ACGT]", tr1)
x[o,2] = str_replace_all(x[o,2], "[ACGT]", tr1)

x[,1] = str_replace_all(x[,1], "[ACGT]", tr2)
x[,2] = str_replace_all(x[,2], "[ACGT]", tr2)

data = table(str_c(x[,1], " > ", x[,2]))
data = data[c(2,4,1,3,6,5)]

# 6.30 × 3.54

pdf(pdfname, family="Helvetica",
    width=2.25, height=2.25,pointsize=8)
par(mai=c(0.5,0.2,0.02,0))

cz <- brewer.pal(3, "Set1")
bp = barplot(data,col=cz[c(3,3,2,2,2,2)],
    axes=FALSE,ylim=c(0,42), xaxt="n")
par(mgp=c(1,0.5,0))
axis(2)
text(x=bp[,1]-0.5, y=-6, names(data), cex=1.2^-1, xpd=TRUE,srt=45)

# text(x=2,y=37, sprintf("%d mutations called",sum(data)), cex=1.2^1,pos=4)
# text(x=2,y=30, sprintf("%0.2f mutation called / meter", sum(data)/90.09),
#     cex=1.2^1, pos=4)

invisible(dev.off())
embedFonts(pdfname, options="-DPDFSETTINGS=/prepress")


# The width of figures, when printed, will usually be
# 5.5 cm (2.25 inches or 1 column) or 12.0 cm
# (4.75 inches or 2 columns). Bar graphs, simple line graphs,
# and gels may be reduced to a smaller width. Symbols and
# lettering should be large enough to be legible after
# reduction [a reduced size of about 7 points (2 mm) high,
# and not smaller than 5 points]. Avoid wide variation in
# type size within a single figure. In laying out information
# in a figure, the objective is to maximize the space given
# to presentation of the data. Avoid wasted white space and
# clutter.