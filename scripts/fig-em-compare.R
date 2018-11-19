library(ape)
library(phytools)
library(stringr)
library(data.table)
library(RColorBrewer)

add_alpha <- function(col, alpha=1){
    apply(sapply(col, col2rgb)/255, 2, function(x) {
        rgb(x[1], x[2], x[3], alpha=alpha)
    })
}

scriptname = commandArgs() %>% str_subset("^--file=") %>% str_replace("^--file=","")
pdfname = scriptname %>% str_replace("\\.[^.]*$", ".pdf")

the_tree = "((((M1:698,M2:638)MD:75,M3:796)MC:773,M4:844)MB:372,(((M7:978,M6:928)MG:111,M8:1029)MF:112,M5:1297)ME:358)MA:0;"
a = read.tree(text=the_tree)
par_tree = "((((M1:0.0411,M2:0.19178)MD:0.0411,M3:0.0274)MC:0,M4:0.10959)MB:0.08904,(((M7:0.0411,M6:0.0411)MG:0.0274,M8:0.16438)MF:0.0274,M5:0.17808)ME:0.08904)MA:0;"
b = read.tree(text=par_tree)

a$tip.label =a$tip.label %>% str_replace("^M","")
b$tip.label =b$tip.label %>% str_replace("^M","")

ab = cophylo(a,b,rotate.multi=TRUE)

#usable slide area 4.25279in x 3.56883in

pdf(pdfname, width=4.25, height=2.75,family="Helvetica")
par(mai=c(0.2,0.2,0.2,0.2),lwd=0.1)

cz <- brewer.pal(3, "Set1")
cz = add_alpha(cz,1.0)

plot(ab,lwd=3,
	edge.col=list(left=rep(cz[2],length(a$edge.length)),right=rep(cz[3],length(b$edge.length)))
)

invisible(dev.off())
embedFonts(pdfname, options="-DPDFSETTINGS=/prepress")
