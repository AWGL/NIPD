# Christopher Medway
# functions for NIPD analysis with Relative Haplotype Doseage

library(ggplot2)
library(reshape2)

plotSPRT <- function(h1, h2, q0, q1) {
    
    # generates a SPRT plot given:
    # h1: vector of allele counts supporting haplotype-1
    # h2: vector of allele counts supporting haplotype-2
    # q0: expected frequency of h1 in maternal plasma given that h2 is inherited by the fetus
    # q1: expected frequency of h1 in maternal plasma given that h1 is inherited by the fetus
    
    h1_cum <- cumsum(h1)
    h2_cum <- cumsum(h2)
    
    N <- h1_cum + h2_cum
    
    h1_prop <- h1_cum / N 
    
    # calculate confidence intervals
    d <- ( (1 - q1) / (1 - q0) )
    g <- (q1 * (1 - q0)) / (q0 * (1 - q1))
    
    range_reads <- seq(min(N),(max(N)))
    upper <- vector(mode = "numeric", length = length(range_reads))
    lower <- vector(mode = "numeric", length = length(range_reads))
    
    for (i in seq(range_reads)) {
        print(i)
        upper[i] <- ((log(1200) / range_reads[i]) - log(d)) / log(g)   
        lower[i] <- ((log(1/1200)/ range_reads[i]) - log(d)) / log(g) 
    }
    
    df <- data.frame(range_reads, upper, lower)
    df$h1_prop[match(N,df$range_reads)] <- h1_prop
    
    
    dfmelt <- reshape2::melt(
        data = df,
        measure.vars = c("upper","lower"),
        id.vars = c("range_reads","h1_prop")
    )
    
    p <- ggplot(data = dfmelt, aes(x = as.numeric(range_reads), y = value, split = variable)) +
        geom_line() +
        geom_point(data = dfmelt, aes(x = range_reads, y = h1_prop), shape = 4) +
        ylim(0.35, 0.6)
    
    return(p)
}



calculateQ <- function(fetalFractions, model = 1, assumption = 1) {
    
    # calculates q0 and q1 variables
    # fetal fraction = needs to be approximated from data 
    # model          = [1] x-linked recessive (currently only option)
    # assumption     = [1] male fetus inherits haplotype with mutation (currently only option)
    
    if (model == 1) {
        message("model is x-linked recessive")
        
        if (assumption == 1) {
            message("assuming fetus inherits haplotype with mutation")
            
            out <- lapply(fetalFractions, function(frac) {
                
                q1 <- (0.5 + frac)
                q0 <- (0.5 - frac)
                return(c("fetal_frac" = frac, "q1" = q1, "q0" = q0))
            })
            
            return(out)
        } else {print("unknown assumption")}
    } else {print("unknown model")}
}


calculateH1Reads <- function(q, depth, nsnps, perms) {
    
  # given x number of SNPs at y depth, calculates the probablility of
  # support for H1 at a given fetal fraction(q)
  
  out <- lapply(q, function(x) {
    
    snpPerm <- lapply(seq(nsnps), function(n) {
      
      readsH1 <- rbinom(prob = q, n = perms, size = depth)
      uCI <- readsH1[order(readsH1)][floor((perms / 100) * 95)]
      lCI <- readsH1[order(readsH1)][floor((perms / 100) * 5)]
      
      return(data.frame("snp" = paste0("rs",n) , lCI, uCI, stringsAsFactors = F))
    })
    df <- do.call(rbind, snpPerm)  # now calculate upper and low props
  })

 
  return(out)
}







