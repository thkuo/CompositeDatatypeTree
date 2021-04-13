plotTTree3<- function (ttree, showLabels = TRUE, sample.col= NA, sample.pch= NA, 
                       legend.col=NA, legend.pch= NA, legend.label= NA, 
                       plot_time_interval= NA, plot.title= NA,
                       hidden_counts.col= c(), showMissingLinks = 0) 
{
    nam = ttree$nam
    ttree = ttree$ttree
    ttree = cbind(ttree, rep(1, nrow(ttree)))
    if (showMissingLinks > 0) {
        i = which(is.na(ttree[, 2]))[1]
        while (i < nrow(ttree)) {
            w = which(ttree[, 3] == i)
            if (length(w) == 1) {
                ttree[w, 3] = ttree[i, 3]
                ttree[w, 4] = ttree[w, 4] + ttree[i, 4]
                ttree = ttree[-i, ]
                ttree[which(ttree[, 3] > i), 3] = ttree[which(ttree[, 
                                                                    3] > i), 3] - 1
            }
            else i = i + 1
        }
    }
    if (showMissingLinks == 2) {
        ttree[which(ttree[, 4] >= 2), 4] = 2
    }
    n = nrow(ttree)
    ys <- rep(0, n)
    scale <- rep(1, n)
    todo = c(which(ttree[, 3] == 0))
    while (length(todo) > 0) {
        f = which(ttree[, 3] == todo[1])
        o = rank(-ttree[f, 1])
        f[o] = f
        for (i in f) {
            ys[i] = ys[todo[1]] + scale[todo[1]] * which(f == 
                                                             i)/(length(f) + 1)
            scale[i] = scale[todo[1]]/(length(f) + 1)
            todo = c(todo, i)
        }
        todo = todo[-1]
    }
    ys = rank(ys)
    par(yaxt = "n", bty = "n")
    mi = min(ttree[which(!is.na(ttree[, 1])), 1])
    ma = max(ttree[which(!is.na(ttree[, 1])), 1])
    xlim= c(mi - (ma - mi) * 0.05, ma + (ma - mi) * 0.05)
    if ((!is.na(plot_time_interval)) && (length(plot_time_interval) ==2)){
        xlim<- (plot_time_interval)
    }
    plot(c(), c(), 
         xlim = xlim, 
         ylim = c(0, n + 1), xlab = "year", ylab = "", xaxt='n')
    if (xlim[2]-xlim[1] >= 2){
        axis(side = 1, at= seq(from= floor(xlim[1]), to= ceiling(xlim[2]), by=1))   
    }else{
        axis(side = 1, at= seq(from= floor(xlim[1]), to= ceiling(xlim[2]), by=.2))   
    }
    pal = gray.colors(max(ttree[, 4]))
    if((length(hidden_counts.col)>0) & (length(hidden_counts.col) >= length(pal))){
        pal<- hidden_counts.col
    }else if (length(hidden_counts.col) < length(pal)){
        sprintf('%d colors in hidden_counts.col but %d needed. Using default instead', length(hidden_counts.col) , length(pal))
    }
    for (i in 1:n) {
        if (ttree[i, 3] != 0) {
            dircol = pal[ttree[i, 4]]
            arrows(ttree[ttree[i, 3], 1], 
                   ys[ttree[i, 3]], 
                   ttree[i, 1], 
                   ys[i], 
                   length = 0, col = dircol)
        }
        if (showLabels && !is.na(ttree[i, 2])) 
            text(ttree[i, 1], ys[i], nam[i], pos = 4, cex = 0.5)
    }
    for (i in 1:n) {
        points(ttree[i, 1], 
               ys[i], 
               # pch = 21, 
               pch= ifelse(is.na(ttree[i, 2]), 21,
                           ifelse(is.na(sample.pch[i]),21, sample.pch[i])),
               bg = ifelse(is.na(ttree[i, 2]), "white", "black"),
               col = ifelse(is.na(ttree[i, 2]), "black",
                           sample.col[i]),
               cex = 1
               # cex = 0.5
               )
    }
    if (length(pal) >= 2) {
        legend("top", legend = 0:(length(pal) - 1), col = pal, 
               lty = 1, cex = 0.5, title = "Missing links")
    }
    if (all(c(!is.na(legend.col), !is.na(legend.pch), !is.na(legend.label)))){
        legend("topleft", 
               pch = legend.pch, 
               col = legend.col, 
               legend = legend.label)    
    }
    if (!is.na(plot.title)){
        legend("topleft", plot.title, bty="n")
    }
    
}