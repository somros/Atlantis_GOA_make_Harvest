biom <- 1000
catch <- 300
exploit <- catch/biom
exploit

fmort <- -1 * log(1-exploit)
fmort

1-exp(-fmort) # this is exploit, the proportion of pop caught in a year

biom1 <- biom * exp(-fmort)
biom1

fday <- fmort / 365
biomday <- c(biom, rep(0,365))
for(i in 2:366){biomday[i] <- biomday[i-1] * exp(-fday)}
biomday

# test
biom * exp(-fmort)

mfc <- 1-exp(-fmort / 365)

biom2 <- c(biom, rep(0,365))
for(i in 2:366){biom2[i] <- biom2[i-1]-biom2[i-1]*mfc}
biom2
# 
# mfc <- 1-exp(-fmort / 365)
# 
# 1-mfc <- exp(-fmort / 365)
# 
# log(1 - mfc) <- (-fmort / 365)
# 
# fmort <- -365 * log(1 - 0.0001)

