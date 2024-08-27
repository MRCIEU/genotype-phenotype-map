library(dplyr)
library(data.table)
library(parallel)

fn <- function(f, suf="") {
    vc <- paste0(f, ".unphased.vcor1")
    if(!file.exists(vc)) {
        message("no ", vc)
        return(NULL)
    }
    bim <- paste0(f, ".bim")

    a <- fread(vc, header=FALSE) %>% as.matrix
    b <- fread(bim) %>% as_tibble()

    ind <- (b$V5 > b$V6)
    i <- ind
    ind <- as.numeric(ind)
    ind[ind == 0] <- -1

    fm <- ind %*% t(ind)
    # dim(fm)
    # fm[1:10,1:10]

    stopifnot(all(diag(fm) == 1))

    a <- a * fm
    temp <- b$V6[i]
    b$V6[i] <- b$V5[i]
    b$V5[i] <- temp

    fwrite(a, file=paste0(vc, suf, ".gz"), row=F, col=F, qu=F, sep=" ")
    system(paste0("cp ", bim, " ", bim, ".orig"))
    fwrite(b, file=paste0(bim, suf), row=F, col=F, qu=F, sep="\t")
    return(NULL)
}

eur <- fread("eur_ldregions.txt", header=FALSE) %>% as_tibble()
fp <- file.path("temp", eur$V1, paste0(eur$V2, "-", eur$V3))

mclapply(fp, fn)


# Check
s <- tibble(
    fp=fp,
    a=sapply(paste0(fp,".unphased.vcor1"), file.size),
    gz=sapply(paste0(fp,".unphased.vcor1.gz"), file.size),
    rat=gz/a
)

summary(s$a)
summary(s$rat)

todo <- subset(s, is.na(rat))
lapply(todo$fp, fn)
