pcavars <- function(vars, name, proportion = 0.95) {
    v <- values(vars)
    n <- sum(complete.cases(v))
    # Sample the study area with n-non.na and creates an environmental table
    sr <- sampleRandom(vars, n)
    pca <- prcomp(scale(sr))
    summary.pca <- summary(pca)
    axis.nb <- which(summary.pca$importance["Cumulative Proportion",] >= proportion)[1]
    pcavars <- predict(vars, pca, index = 1:axis.nb)
    return(pcavars)
    }
