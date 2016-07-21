
## Need more work.


# a new method of associatin between genome tracks
# ~ LD50  CC50  (Coverage Concordance)  DC50 (Distance Concordance)   Cannot use GC50
# How much do the regions in track A need to be expanded to cover 50% of the regions in track B.

# repeatedly run findOverlaps with increasing maxgap parameter.
# Express this as % or function of the coverage of one or both tracks or % of the genome.

# Currently seems quite slow but does work. Probably slowed by search space? 
# Could add in heuristic based on likely values. (1000, 10000)  # not much quicker

# TODO
# also need to scale result (or add extra part) to account for coverage of original tracks.
#   and the size of the genome?  Concerned that a track covering most/all the genome should score highly.
#   Beginning to sound a bit like concordance when meant to be simple. 
cc50 <- function(track.A, track.B, cc=.5, interval = c(1, 1000000), lower = min(interval),
                 upper = max(interval))  {
    require(GenomicRanges)
    
    if(length(track.A) <= length(track.B)) {
        short.track <- track.A
        long.track <- track.B
    } else {
        short.track <- track.B
        long.track <- track.A
    }
    
    # need a function to minimise or maximise. (maximum=TRUE). Can give upper limit.
    f <- function (x, short.track, long.track, cc=.5) {
        overlap.count <- sum(overlapsAny(short.track, long.track, maxgap=x))# how many of the short track are overlapped by the long track.
        return(abs((length(short.track)*cc) - overlap.count))  # absolute difference from half the length of the shorter track.
    }
    xmin <- optimize(f, tol = 1, short.track=short.track, long.track=long.track, lower=lower, upper=upper, cc=cc)
    #xmin$cc <- xmin$minimum / max(sum(width(track.A)), sum(width(track.B)))
    xmin$cc <- log10(xmin$minimum * max(sum(width(track.A)), sum(width(track.B))))
    if(xmin$minimum >= (upper *0.95)) warning("cc50 near upper limit of search range, consider extending interval")
    return(xmin)
}


#cc50(emb.GR, gwas.GR)
#cc50(emb.GR, f5.GR)
#cc50(emb.GR, f5.GR, cc=.25)
#cc50(emb.GR, f5.GR, cc=.75)
#cc50(f5.GR, shift(f5.GR, 10000))

#cc50(emb.GR, gwas.GR, cc=.25)
#cc50(emb.GR, gwas.GR, cc=.75)
#system.time(cc50(emb.GR, gwas.GR, interval = c(0, 1000000)))
#system.time(cc50(emb.GR, gwas.GR, interval = c(1000, 10000)))

#sum(overlapsAny(emb.GR, gwas.GR))
#sum(overlapsAny(gwas.GR, emb.GR))
