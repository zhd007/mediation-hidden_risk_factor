library("tidyverse")
library("reshape2")
library("dplyr")
library("mixtools")
library("ggplot2")
library("SimDesign")
library(msm)
library(fitdistrplus)
library(Matrix)
library(magic)

# This is a script to simulate mediation analysis for continuous outcomes.


#set.seed(2020)
# this is a vector of snp names: e.g. snp1-snp50 
num_snp <- 500               # number of snp in the model
num_boot <- 2000
set <- 100
em.times <- 21

known_freq <- c(0.5, 0.6, 0.8, 0.2, 0.5, 0.8 )    # freq of m1, m2, m3, m4, m5, m_Hidden
known_beta <- c(0.4, 0.2, 0.3, 0.2, 0.4)            # beta1 and beta2, beta3, beat4, beta5

sample_size <- 10000

# c_overall <- c(0)
# a_overall <- c(0)

true_hidden <- runif(n = set, min = 0.02, max = 0.5)
#true_hidden <- rep(0, times = set)
# the snp names
snp_col <- c()
for (i in 1:num_snp) {
  snp_col[i] <- paste("snp", i, sep="")
}

# FiellerRatioCI_basic <- function(a,b,V,alpha=0.05){
#   theta <- a/b
#   v11 <- V[1,1]
#   v12 <- V[1,2]
#   v22 <- V[2,2]
#   
#   z <- qnorm(1-alpha/2)
#   g <- (z^2)*v22/b^2
#   C <- sqrt(v11 - 2*theta*v12 + theta^2 * v22 - g*(v11-v12^2/v22))
#   minS <- (1/(1-g))*(theta- g*v12/v22 - z/b * C)
#   maxS <- (1/(1-g))*(theta- g*v12/v22 + z/b * C)
#   return(c(ratio=theta,min=minS,max=maxS))
# }



# ------------------------------------- simulation -------------------------------------


record <- data.frame("freq" = -9, "b_hat" = -9, "est_a" = -9, "true_freq" = -9, 
                     "b_hat.null" = -9, "b_hat.alt" = -9, "b_hat.ml" = -9, 
                     "b.low" = -9, "b.high" = -9,
                     "b.low2" = -9, "b.high2" = -9,
                     "bt.alt.low" = -9, "bt.alt.high" = -9,
                     "bt.null.low" = -9, "bt.null.high" = -9,
                     "separate" = -9, "method" = -9, "true_hidden" = -9)

for (k in 1:set) {
  tryCatch({
    print(paste("------------------------ Iteration", k, "---------------------------"))
    
    df.snp <- data.frame(matrix(vector(), sample_size, length(snp_col),
                                dimnames=list(c(), snp_col)),
                         stringsAsFactors=F)
    # generate snp
    for (i in 1:length(snp_col)) {
      #### frequency of the SNPs
      #freq <- runif(n = 1, min = 0.1, max = 0.4)
      freq <- runif(n = 1, min = 0.1, max = 0.9)
      #freq <- runif(n = 1, min = 0.005, max = 0.015)
      cur_snp_v <- sample(c(0, 1, 2), size = sample_size, replace = TRUE, prob = c((1-freq)^2, 2*freq*(1-freq), freq^2))
      cur_snp_name <- paste(snp_col[i])
      assign(cur_snp_name, cur_snp_v)
      df.snp[[cur_snp_name]] <- cur_snp_v
    }
    
    
    
    
    #cor12 <- mvrnorm(n = sample_size, mu = c(50,5), Sigma = matrix(c(2,1.5,1.5,2), nrow = 2, ncol = 2))
    #cor34 <- mvrnorm(n = sample_size, mu = c(10,6), Sigma = matrix(c(3,1.2,1.2,3), nrow = 2, ncol = 2))
    
    #betaS1_2 <- mvrnorm(n = num_snp, mu = c(0.2,0.2), Sigma = matrix(c(0.05,0.005,0.005,0.05), nrow = 2, ncol = 2))
    
    upper_level_effect <- rnorm(n = length(known_beta)+1, 0.2, 0.01)
    
    ############# known: m1
  
    # a vector of indicator of presence
    seq1 <- sample(c(0, 1), size = num_snp, replace = TRUE, prob = c(1-known_freq[1], known_freq[1]))
    #seq1 <- rep(1, num_snp)
    # a vector of beta(snp->mediator)
    #betaS1 <- rnorm(n = num_snp, mean = 0.3, sd = 0.1)
    betaS1 <- rnorm(n = num_snp, mean = upper_level_effect[1], sd = 0.04)
    #betaS1 <- betaS1_2[,1]
    
    m1 <- data.matrix(df.snp) %*% (seq1 * betaS1)
    m1 <- m1[,1]
    #m1 <- m1 + 50 # constant
    m1 <- m1 + rnorm(sample_size, mean = 0, sd = 2) # add some noise: size of noise?
    #m1 <- m1 + cor12[,1]
    #m1 <- m1 + runif(sample_size, min = -0.3, max = 0.3)
    
    ############# known: m2
    seq2 <- sample(c(0, 1), size = num_snp, replace = TRUE, prob = c(1-known_freq[2], known_freq[2]))
    #seq2 <- rep(1, num_snp)
    # a vector of beta(snp->mediator)
    #betaS2 <- rnorm(n = num_snp, mean = 0.3, sd = 0.1)
    betaS2 <- rnorm(n = num_snp, mean = upper_level_effect[2], sd = 0.04)
    #betaS2 <- betaS1_2[,2]
    
    m2 <- data.matrix(df.snp) %*% (seq2 * betaS2)
    m2 <- m2[,1]
    #m2 <- m2 + 5 # constant
    m2 <- m2 + rnorm(sample_size, mean = 0, sd = 3) + 0.9 * m1
    #m2 <- m2 + cor12[,2]
    #m2 <- m2 + runif(sample_size, min = -0.4, max = 0.4)
    
    
    ############# known: m3
    seq3 <- sample(c(0, 1), size = num_snp, replace = TRUE, prob = c(1-known_freq[3], known_freq[3]))
    #seq3 <- rep(1, num_snp)
    # a vector of beta(snp->mediator)
    #betaS3 <- rnorm(n = num_snp, mean = 0.3, sd = 0.1)
    betaS3 <- rnorm(n = num_snp, mean = upper_level_effect[3], sd = 0.04)
    
    m3 <- data.matrix(df.snp) %*% (seq3 * betaS3)
    m3 <- m3[,1]
    #m3 <- m3 + 10 # constant
    m3 <- m3 + rnorm(sample_size, mean = 0, sd = 4)
    #m3 <- m3 + cor34[,1]
    #m3 <- m3 + runif(sample_size, min = -0.4, max = 0.4)
    
    ############# known: m4
    seq4 <- sample(c(0, 1), size = num_snp, replace = TRUE, prob = c(1-known_freq[4], known_freq[4]))
    #seq4 <- rep(1, num_snp)
    # a vector of beta(snp->mediator)
    #betaS4 <- rnorm(n = num_snp, mean = 0.3, sd = 0.1)
    betaS4 <- rnorm(n = num_snp, mean = upper_level_effect[4], sd = 0.04)
    
    m4 <- data.matrix(df.snp) %*% (seq4 * betaS4)
    m4 <- m4[,1]
    #m4 <- m4 + 6 # constant
    m4 <- m4 + rnorm(sample_size, mean = 0, sd = 3)
    #m4 <- m4 + cor34[,2]
    #m4 <- m4 + runif(sample_size, min = -0.4, max = 0.4)
    
    ############# known: m5
    seq5 <- sample(c(0, 1), size = num_snp, replace = TRUE, prob = c(1-known_freq[5], known_freq[5]))
    #seq5 <- rep(1, num_snp)
    # a vector of beta(snp->mediator)
    #betaS5 <- rnorm(n = num_snp, mean = 0.3, sd = 0.1)
    betaS5 <- rnorm(n = num_snp, mean = upper_level_effect[5], sd = 0.04)
    
    m5 <- data.matrix(df.snp) %*% (seq5 * betaS5)
    m5 <- m5[,1]
    m5 <- m5 + 15 # constant
    m5 <- m5 + rnorm(sample_size, mean = 0, sd = 4)
    #m5 <- m5 + runif(sample_size, min = -0.4, max = 0.4)
    
    ############# hidden: m6
    seq6 <- sample(c(0, 1), size = num_snp, replace = TRUE, prob = c(1-known_freq[6], known_freq[6]))
    # a vector of beta(snp->mediator)
    #betaS6 <- rnorm(n = num_snp, mean = 0.3, sd = 0.1)
    betaS6 <- rnorm(n = num_snp, mean = upper_level_effect[6], sd = 0.04)
    
    m6 <- data.matrix(df.snp) %*% (seq6 * betaS6)
    m6 <- m6[,1]
    m6 <- m6 + 20 # constant
    m6 <- m6 + rnorm(sample_size, mean = 0, sd = 2) 
    #m6 <- m6 + runif(sample_size, min = -0.3, max = 0.3)
    
    ############ covariates
    covariate1 <- rnorm(sample_size, mean = 7, sd = 0.5)
    covariate2 <- rnorm(sample_size, mean = 4, sd = 0.4)
    
    
    ############## mediator -> outcome
    y <- known_beta[1] * m1 + known_beta[2] * m2 + known_beta[3] * m3 + 
      known_beta[4] * m4 + known_beta[5] * m5 +
      true_hidden[k] * m6 + rnorm(sample_size, mean = 0, sd = 0.2) +
      0.2 * covariate1 - 0.1 * covariate2
    # y <- known_beta[1] * m1 + known_beta[2] * m2 + known_beta[3] * m3 + 
    #   known_beta[4] * m4 + known_beta[5] * m5 +
    #   true_hidden[k] * m6 + runif(sample_size, min = -0.3, max = 0.3)
    
    # combined
    df.snp$m1 <- m1
    df.snp$m2 <- m2
    df.snp$m3 <- m3
    df.snp$m4 <- m4
    df.snp$m5 <- m5
    df.snp$m6 <- m6
    df.snp$y <- y
    df.snp$covariate1 <- covariate1
    df.snp$covariate2 <- covariate2
    
    
    para_decision <- vector(mode = "logical", length = num_snp)
    beta_m1 <- vector(mode = "logical", length = num_snp)
    beta_m2 <- vector(mode = "logical", length = num_snp)
    beta_m3 <- vector(mode = "logical", length = num_snp)
    beta_m4 <- vector(mode = "logical", length = num_snp)
    beta_m5 <- vector(mode = "logical", length = num_snp)
    
    # regression
    # snp -> mediators
    # df.a <- data.frame("a1" = -9, "a2" = -9, "a3" = -9,
    #                    "a1.var" = -9, "a2.var" = -9, "a3.var" = -9)
    # for (i in 1:num_snp) {
    #   eq.m1 <- paste("m1~", snp_col[i], sep = "")
    #   eq.m2 <- paste("m2~", snp_col[i], sep = "")
    #   eq.m3 <- paste("m3~", snp_col[i], sep = "")
    #   
    #   row.a <- data.frame("a1" = summary(lm(eq.m1, data = df.snp))$coefficient[2, 1],
    #                       "a2" = summary(lm(eq.m2, data = df.snp))$coefficient[2, 1],
    #                       "a3" = summary(lm(eq.m3, data = df.snp))$coefficient[2, 1],
    #                       "a1.var" = summary(lm(eq.m1, data = df.snp))$coefficient[2, 2]^2,
    #                       "a2.var" = summary(lm(eq.m1, data = df.snp))$coefficient[2, 2]^2,
    #                       "a3.var" = summary(lm(eq.m1, data = df.snp))$coefficient[2, 2]^2)
    #   df.a <- rbind(df.a, row.a)
    # }
    # df.a <- df.a[df.a$a1 != -9,]
    
    # alternative way
    # 70 snps
    # eq.m1.all <- "m1~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+snp51+snp52+snp53+snp54+snp55+snp56+snp57+snp58+snp59+snp60+snp61+snp62+snp63+snp64+snp65+snp66+snp67+snp68+snp69+snp70"
    # eq.m2.all <- "m2~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+snp51+snp52+snp53+snp54+snp55+snp56+snp57+snp58+snp59+snp60+snp61+snp62+snp63+snp64+snp65+snp66+snp67+snp68+snp69+snp70"
    # eq.m3.all <- "m3~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+snp51+snp52+snp53+snp54+snp55+snp56+snp57+snp58+snp59+snp60+snp61+snp62+snp63+snp64+snp65+snp66+snp67+snp68+snp69+snp70"
    
    
    # 500 snps
    snp_sum <- paste("snp", 1:num_snp, sep = "", collpase = "")
    snp_sum.cat <- do.call(paste, c(as.list(snp_sum), sep = "+"))
    eq.m1.all <- paste("m1~", snp_sum.cat, sep = "")
    eq.m2.all <- paste("m2~", snp_sum.cat, sep = "")
    eq.m3.all <- paste("m3~", snp_sum.cat, sep = "")
    eq.m4.all <- paste("m4~", snp_sum.cat, sep = "")
    eq.m5.all <- paste("m5~", snp_sum.cat, sep = "")
    eq.m6.all <- paste("m6~", snp_sum.cat, sep = "")
    
    m1_model <- lm(eq.m1.all, data = df.snp)
    m2_model <- lm(eq.m2.all, data = df.snp)
    m3_model <- lm(eq.m3.all, data = df.snp)
    m4_model <- lm(eq.m4.all, data = df.snp)
    m5_model <- lm(eq.m5.all, data = df.snp)
    m6_model <- lm(eq.m6.all, data = df.snp)
    
    
    m1_coefs <- summary(m1_model)$coefficient
    m2_coefs <- summary(m2_model)$coefficient
    m3_coefs <- summary(m3_model)$coefficient
    m4_coefs <- summary(m4_model)$coefficient
    m5_coefs <- summary(m5_model)$coefficient
    m6_coefs <- summary(m6_model)$coefficient
    
    df.a <- data.frame("a1.alt" = m1_coefs[2:(num_snp+1),1],
                       "a2.alt" = m2_coefs[2:(num_snp+1),1], 
                       "a3.alt" = m3_coefs[2:(num_snp+1),1],
                       "a4.alt" = m4_coefs[2:(num_snp+1),1],
                       "a5.alt" = m5_coefs[2:(num_snp+1),1],
                       "a6.alt" = m6_coefs[2:(num_snp+1),1],
                       "a1.var.alt" = m1_coefs[2:(num_snp+1),2]^2,
                       "a2.var.alt" = m2_coefs[2:(num_snp+1),2]^2, 
                       "a3.var.alt" = m3_coefs[2:(num_snp+1),2]^2,
                       "a4.var.alt" = m4_coefs[2:(num_snp+1),2]^2,
                       "a5.var.alt" = m5_coefs[2:(num_snp+1),2]^2,
                       "a6.var.alt" = m6_coefs[2:(num_snp+1),2]^2)
    
    
    # plot(df.a$a1, df.a$a1.alt)
    # plot(df.a$a2, df.a$a2.alt)
    # plot(df.a$a3, df.a$a3.alt)
    # 
    # plot(df.a$a1.p, df.a$a1.p.alt)
    # plot(df.a$a2.p, df.a$a2.p.alt)
    # plot(df.a$a3.p, df.a$a3.p.alt)
    
    # mediators -> outcome
    # df.b <- data.frame("b1" = -9, "b2" = -9, "c" = -9,
    #                    "b1.var" = -9, "b2.var" = -9, "c.var" = -9)
    # for (i in 1:num_snp) {
    #   eq.y <- paste("y~m1+m2+", snp_col[i], sep = "")
    #   row.b<- data.frame("b1" = summary( lm(eq.y), data = df.snp)$coefficient[2,1], 
    #                      "b2" = summary( lm(eq.y), data = df.snp)$coefficient[3,1],
    #                      "c" = summary( lm(eq.y), data = df.snp)$coefficient[4,1],
    #                      "b1.var" = summary( lm(eq.y), data = df.snp)$coefficient[2,2]^2, 
    #                      "b2.var" = summary( lm(eq.y), data = df.snp)$coefficient[3,2]^2,
    #                      "c.var" = summary( lm(eq.y), data = df.snp)$coefficient[4,2]^2)
    #   df.b <- rbind(df.b, row.b)
    # }
    # df.b <- df.b[df.b$b1 != -9,]
    
    # alternative way
    #eq.y.all <- "y~m1+m2 + snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+snp51+snp52+snp53+snp54+snp55+snp56+snp57+snp58+snp59+snp60+snp61+snp62+snp63+snp64+snp65+snp66+snp67+snp68+snp69+snp70"
    eq.y.all <- paste("y~m1+m2+m3+m4+m5+", snp_sum.cat, "+covariate1+covariate2", sep = "")
    
    outcome_model <- lm(eq.y.all, data = df.snp)
    y_coefs <- summary(outcome_model)$coefficient
    
    df.b <- data.frame("b1.alt" = y_coefs[2,1], "b2.alt" = y_coefs[3,1], "b3.alt" = y_coefs[4,1], 
                       "b4.alt" = y_coefs[5,1], "b5.alt" = y_coefs[6,1], 
                       "c.alt" = y_coefs[(length(known_beta)+2):(num_snp+1+length(known_beta)),1],
                       "b1.var.alt" = y_coefs[2,2]^2, "b2.var.alt" = y_coefs[3,2]^2, "b3.var.alt" = y_coefs[4,2]^2,
                       "b4.var.alt" = y_coefs[5,2]^2, "b5.var.alt" = y_coefs[6,2]^2,
                       "c.var.alt" = y_coefs[(length(known_beta)+2):(num_snp+1+length(known_beta)),2]^2   )
    
    # plot(df.b$b1, df.b$b1.alt)
    # plot(df.b$b2, df.b$b2.alt)
    # plot(df.b$c, df.b$c.alt)
    
    
    
    # try using em mixtool
    separate <- "yes"    # whether the c distribution is separable
    
    
    
    ####### testing null: no bound at all
    ####### for a
    emset.a <- em.times 
    emLR.a <- vector(mode = "logical", length = emset.a)
    emMu1.a <- vector(mode = "logical", length = emset.a)
    emMu2.a <- vector(mode = "logical", length = emset.a)
    emSd1.a <- vector(mode = "logical", length = emset.a)
    emSd2.a <- vector(mode = "logical", length = emset.a)
    emLambda1.a <- vector(mode = "logical", length = emset.a)
    emLambda2.a <- vector(mode = "logical", length = emset.a)
    emModel.a <- vector(mode = "list", length = emset.a)
    
    for (em.a in 1:emset.a) {
      snp_em.status.a <- -999
      while (snp_em.status.a == -999) {
        tryCatch({
          snp_em.a <- quiet(normalmixEM(c(df.a$a1.alt,df.a$a2.alt, df.a$a3.alt, df.a$a4.alt, df.a$a5.alt), 
                                        lambda = .5, arbvar = TRUE), messages = TRUE)
          snp_em.status.a <- 1
        }, error=function(e){
          print("caught")
          snp_em.status.a <- -999})
      }
      # plot(snp_em, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
      #      main2="SNP", xlab2="Beta")
      emLR.a[em.a] <- snp_em.a$loglik
      emMu1.a[em.a] <- snp_em.a$mu[1]
      emSd1.a[em.a] <- snp_em.a$sigma[1]
      emLambda1.a[em.a] <- snp_em.a$lambda[1]
      emMu2.a[em.a] <- snp_em.a$mu[2]
      emSd2.a[em.a] <- snp_em.a$sigma[2]
      emLambda2.a[em.a] <- snp_em.a$lambda[2]
      emModel.a[[em.a]] <- snp_em.a
    }
    
    
    
    #maxI.a <- which(emLR.a == max(emLR.a))[1]
    maxI.a <- which(round(emMu2.a, digits = 6) == round(median(emMu2.a), digits = 6))[1]
    snp_em_summary <- data.frame("mean" = c(emMu1.a[maxI.a], emMu2.a[maxI.a]), 
                                 "sd" = c(emSd1.a[maxI.a], emSd2.a[maxI.a]), 
                                 "lambda" = c(emLambda1.a[maxI.a], emLambda2.a[maxI.a]) )
    medianModel.a <- emModel.a[[maxI.a]]
    
    ####### for c
    emset.c <- em.times 
    emLR.c <- vector(mode = "logical", length = emset.c)
    emMu1.c <- vector(mode = "logical", length = emset.c)
    emMu2.c <- vector(mode = "logical", length = emset.c)
    emSd1.c <- vector(mode = "logical", length = emset.c)
    emSd2.c <- vector(mode = "logical", length = emset.c)
    emLambda1.c <- vector(mode = "logical", length = emset.c)
    emLambda2.c <- vector(mode = "logical", length = emset.c)
    emModel.c <- vector(mode = "list", length = emset.c)
    
    
    for (em.c in 1:emset.c) { 
      snp_em.status.c <- -999
      while (snp_em.status.c == -999) {
        tryCatch({
          snp_em.c <- quiet(normalmixEM(df.b$c.alt, lambda = .5, arbvar = TRUE), messages = TRUE)
          snp_em.status.c <- 1
        }, error=function(e){
          print("caught")
          snp_em.status.c <- -999})
      }
      # plot(snp_em.c, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
      #      main2="SNP", xlab2="Beta")
      emLR.c[em.c] <- snp_em.c$loglik
      emMu1.c[em.c] <- snp_em.c$mu[1]
      emSd1.c[em.c] <- snp_em.c$sigma[1]
      emLambda1.c[em.c] <- snp_em.c$lambda[1]
      emMu2.c[em.c] <- snp_em.c$mu[2]
      emSd2.c[em.c] <- snp_em.c$sigma[2]
      emLambda2.c[em.c] <- snp_em.c$lambda[2]
      emModel.c[[em.c]] <- snp_em.c
    }
    #maxI.c <- which(emLR.c == max(emLR.c))[1]
    maxI.c <- which(round(emMu2.c, digits = 6) == round(median(emMu2.c), digits = 6))[1]
    snp_em_summary.c <- data.frame("mean" = c(emMu1.c[maxI.c], emMu2.c[maxI.c]), 
                                   "sd" = c(emSd1.c[maxI.c], emSd2.c[maxI.c]), 
                                   "lambda" = c(emLambda1.c[maxI.c], emLambda2.c[maxI.c]) )
    medianModel.c <- emModel.c[[maxI.c]]
    
    # for testing null
    b_hat.null <- mean(df.b$c.alt)/ snp_em_summary$mean[2]
    
    # for estimation 
    b_hat.alt <- snp_em_summary.c$mean[2] / snp_em_summary$mean[2]

    
    # ###### bootstrap from all_a and all_c - bootstrap SNPs
    # boot_result.alt <- vector(mode = "logical", length = 1000)
    # boot_result.null <- vector(mode = "logical", length = 1000)
    # for (v in 1:1000) {
    #   if(v %% 100 == 0) {
    #     print(paste("##### bt round:", v, " #####", sep = ""))
    #   }
    #   bt_index <- sample(1:length(df.b$c.alt), size = length(df.b$c.alt), replace = TRUE)
    #   boot_a_sample <- c(df.a$a1.alt[bt_index], df.a$a2.alt[bt_index], df.a$a3.alt[bt_index],
    #                      df.a$a4.alt[bt_index], df.a$a5.alt[bt_index])
    #   boot_c_sample <- df.b$c.alt[bt_index]
    #   ####### for a
    #   emset.a <- em.times
    #   emLR.a <- vector(mode = "logical", length = emset.a)
    #   emMu1.a <- vector(mode = "logical", length = emset.a)
    #   emMu2.a <- vector(mode = "logical", length = emset.a)
    #   emSd1.a <- vector(mode = "logical", length = emset.a)
    #   emSd2.a <- vector(mode = "logical", length = emset.a)
    #   emLambda1.a <- vector(mode = "logical", length = emset.a)
    #   emLambda2.a <- vector(mode = "logical", length = emset.a)
    #   boot.emModel.a <- vector(mode = "list", length = emset.a)
    # 
    #   for (em.a in 1:emset.a) {
    #     snp_em.status.a <- -999
    #     while (snp_em.status.a == -999) {
    #       tryCatch({
    #         snp_em.a <- quiet(normalmixEM(boot_a_sample,
    #                                       lambda = .5, arbvar = TRUE), messages = TRUE)
    #         snp_em.status.a <- 1
    #       }, error=function(e){
    #         print("caught")
    #         snp_em.status.a <- -999})
    #     }
    #     # plot(snp_em, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
    #     #      main2="SNP", xlab2="Beta")
    #     emLR.a[em.a] <- snp_em.a$loglik
    #     emMu1.a[em.a] <- snp_em.a$mu[1]
    #     emSd1.a[em.a] <- snp_em.a$sigma[1]
    #     emLambda1.a[em.a] <- snp_em.a$lambda[1]
    #     emMu2.a[em.a] <- snp_em.a$mu[2]
    #     emSd2.a[em.a] <- snp_em.a$sigma[2]
    #     emLambda2.a[em.a] <- snp_em.a$lambda[2]
    #     boot.emModel.a[[em.a]] <- snp_em.a
    #   }
    #   #maxI.a <- which(emLR.a == max(emLR.a))[1]
    #   maxI.a <- which(round(emMu2.a, digits = 6) == round(median(emMu2.a), digits = 6))[1]
    #   boot.medianModel.a <- boot.emModel.a[[maxI.a]]
    # 
    #   ####### for c
    #   emset.c <- em.times
    #   emLR.c <- vector(mode = "logical", length = emset.c)
    #   emMu1.c <- vector(mode = "logical", length = emset.c)
    #   emMu2.c <- vector(mode = "logical", length = emset.c)
    #   emSd1.c <- vector(mode = "logical", length = emset.c)
    #   emSd2.c <- vector(mode = "logical", length = emset.c)
    #   emLambda1.c <- vector(mode = "logical", length = emset.c)
    #   emLambda2.c <- vector(mode = "logical", length = emset.c)
    #   boot.emModel.c <- vector(mode = "list", length = emset.c)
    # 
    # 
    #   for (em.c in 1:emset.c) {
    #     snp_em.status.c <- -999
    #     while (snp_em.status.c == -999) {
    #       tryCatch({
    #         snp_em.c <- quiet(normalmixEM(boot_c_sample, lambda = .5, arbvar = TRUE), messages = TRUE)
    #         snp_em.status.c <- 1
    #       }, error=function(e){
    #         print("caught")
    #         snp_em.status.c <- -999})
    #     }
    #     # plot(snp_em.c, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
    #     #      main2="SNP", xlab2="Beta")
    #     emLR.c[em.c] <- snp_em.c$loglik
    #     emMu1.c[em.c] <- snp_em.c$mu[1]
    #     emSd1.c[em.c] <- snp_em.c$sigma[1]
    #     emLambda1.c[em.c] <- snp_em.c$lambda[1]
    #     emMu2.c[em.c] <- snp_em.c$mu[2]
    #     emSd2.c[em.c] <- snp_em.c$sigma[2]
    #     emLambda2.c[em.c] <- snp_em.c$lambda[2]
    #     boot.emModel.c[[em.c]] <- snp_em.c
    #   }
    #   #maxI.c <- which(emLR.c == max(emLR.c))[1]
    #   maxI.c <- which(round(emMu2.c, digits = 6) == round(median(emMu2.c), digits = 6))[1]
    #   boot.medianModel.c <- boot.emModel.c[[maxI.c]]
    # 
    #   # different results based on unimodal or bimodal
    #   d.a  <- abs(boot.medianModel.a$mu[1] - boot.medianModel.a$mu[2]) /
    #     (2 * boot.medianModel.a$sigma[1] * boot.medianModel.a$sigma[2])
    #   d.c  <- abs(boot.medianModel.c$mu[1] - boot.medianModel.c$mu[2]) /
    #     (2 * boot.medianModel.c$sigma[1] * boot.medianModel.c$sigma[2])
    # 
    #   if (d.a <= 1  ) {
    #     bt.est.a <- mean(boot_a_sample)
    #   } else {
    #     bt.est.a <- boot.medianModel.a$mu[2]
    #   }
    #   if (d.c <= 1  ) {
    #     bt.est.c <- mean(boot_c_sample)
    #   } else {
    #     bt.est.c <- boot.medianModel.c$mu[2]
    #   }
    # 
    #   boot_result.alt[v] <- bt.est.c / bt.est.a
    #   
    #   ####### for c: null case
    #   bt.est.c <- mean(boot_c_sample)
    # 
    #   # different results based on unimodal or bimodal
    #   d.a  <- abs(boot.medianModel.a$mu[1] - boot.medianModel.a$mu[2]) /
    #     (2 * boot.medianModel.a$sigma[1] * boot.medianModel.a$sigma[2])
    # 
    #   if (d.a <= 1  ) {
    #     bt.est.a <- mean(boot_a_sample)
    #   } else {
    #     bt.est.a <- boot.medianModel.a$mu[2]
    #   }
    # 
    #   boot_result.null[v] <- bt.est.c / bt.est.a
    # }
    # bt.alt.low <- sort(boot_result.alt)[25]
    # bt.alt.high <- sort(boot_result.alt)[976]
    # 
    # bt.null.low <- sort(boot_result.null)[25]
    # bt.null.high <- sort(boot_result.null)[976]
    
    
    ###### conditional variance
    ###### alt case
    ### a
    p.a <- medianModel.a$posterior[,2]
    m1.vcov <- vcov(m1_model)[2:(num_snp+1), 2:(num_snp+1)]
    m2.vcov <- vcov(m2_model)[2:(num_snp+1), 2:(num_snp+1)]
    m3.vcov <- vcov(m3_model)[2:(num_snp+1), 2:(num_snp+1)]
    m4.vcov <- vcov(m4_model)[2:(num_snp+1), 2:(num_snp+1)]
    m5.vcov <- vcov(m5_model)[2:(num_snp+1), 2:(num_snp+1)]
    var_1.a <- t(p.a) %*% (adiag(m1.vcov, m2.vcov, m3.vcov, m4.vcov, m5.vcov) +
                             diag(rep(snp_em_summary$sd[2]^2, num_snp*length(known_beta)))) %*% p.a
    var_1.a <- var_1.a / ( sum(p.a)^2   )

    ### c
    p.c <- medianModel.c$posterior[,2]
    outcome.vcov <- vcov(outcome_model)[(length(known_beta)+2):(length(known_beta)+num_snp+2-1),
                                        (length(known_beta)+2):(length(known_beta)+num_snp+2-1)]

    var_1.c <- t(p.c) %*% ( outcome.vcov  +  diag(rep(snp_em_summary.c$sd[2]^2, num_snp))  ) %*% p.c
    var_1.c <- var_1.c / ( sum(p.c)^2   )

    delta.index2 <- which(p.a > 0.5)
    gamma.a <- sum(p.a[delta.index2]) / sum(p.a)
    
    #var_2.a <- var(gamma.a[!(gamma.a %in% boxplot(gamma.a)$out)])
    var_2.a <- gamma.a * (1-gamma.a) * (snp_em_summary$mean[2]^2 + snp_em_summary$mean[1]^2 -
                           2*snp_em_summary$mean[2]*snp_em_summary$mean[1])
    #var_2 <- se_gamma^2 * (snp_em_summary$mean[2]^2 )

    total_var.a <- var_1.a + var_2.a
    
    
    delta.index2 <- which(p.c > 0.5)
    gamma.c <- sum(p.c[delta.index2]) / sum(p.c)
    
    # var_2.c <- var(gamma.c[!(gamma.c %in% boxplot(gamma.c)$out)])
    var_2.c <- gamma.c * (1-gamma.c) * (snp_em_summary.c$mean[2]^2 + snp_em_summary.c$mean[1]^2 -
                               2*snp_em_summary.c$mean[2]*snp_em_summary.c$mean[1])
    #var_2.c <- se_gamma.c^2 * (snp_em_summary.c$mean[2]^2 )
    total_var.c <- var_1.c + var_2.c
    
    se.delta <- deltamethod(~ x1 / x2,
                            mean = c(snp_em_summary.c$mean[2], snp_em_summary$mean[2]),
                            cov = matrix(c(total_var.c, 0,0, total_var.a), nrow = 2, ncol = 2))
    
    b.low2 <- b_hat.alt - 1.96 * se.delta 
    b.high2 <- b_hat.alt + 1.96 * se.delta 
    
    ###### null case
    #var.c.null <- sum( var(df.b$c.alt) + df.b$c.var.alt )  / ( length(df.b$c.alt)^2 )
    outcome.vcov <- vcov(outcome_model)[(length(known_beta)+2):(length(known_beta)+num_snp+2-1), 
                                        (length(known_beta)+2):(length(known_beta)+num_snp+2-1)]
    var.c.null <- t(rep(1, length(df.b$c.alt)))  %*% 
      (diag(rep(var(df.b$c.alt), length(df.b$c.alt)))   +  outcome.vcov) %*%
      rep(1, length(df.b$c.alt))/ 
      ( length(df.b$c.alt)^2 )
    se.delta.null <- deltamethod(~ x1 / x2,
                                 mean = c(mean(df.b$c.alt), snp_em_summary$mean[2]),
                                 cov = matrix(c(var.c.null, 0,0, total_var.a), nrow = 2, ncol = 2))
    b.low1 <- b_hat.null - 1.96 * se.delta.null
    b.high1 <- b_hat.null + 1.96 * se.delta.null

    #b_hat.ml <- ifelse(c_dist$loglik > snp_em$loglik, yes = b_hat, no = b_hat.alt)
    #method <- ifelse(c_dist$loglik > snp_em$loglik, yes = "single", no = "bimodals")
    
    #### record result
    temp.record <- data.frame("freq" = snp_em_summary.c$lambda[2],
                              "b_hat" = -1,
                              "est_a" = snp_em_summary$mean[2],
                              "true_freq" = mean(seq3),
                              "b_hat.null" = b_hat.null,
                              "b_hat.alt" = b_hat.alt,
                              "b_hat.ml" = -9,
                              "b.low" = b.low1, "b.high" = b.high1,
                              "b.low2" = b.low2, "b.high2" = b.high2,
                              "bt.alt.low" = -9, "bt.alt.high" = -9,
                              "bt.null.low" = -9, "bt.null.high" = -9,
                              "separate" = separate,
                              "method" = -9,
                              "true_hidden" = true_hidden[k])
    record <- rbind(record,  temp.record)
    print(temp.record)
    
    # c_overall <- c(c_overall, snp_em_summary.c$mean[2])
    # a_overall <- c(a_overall, snp_em_summary$mean[2])  
  }, error=function(e){print("Outter caught")})
}
record <- record[record$freq != -9,]


## beta plot
record$capture1 <- ifelse((record$b.low <= record$true_hidden & record$b.high >= record$true_hidden), TRUE, FALSE)
record$capture2 <- ifelse((record$b.low2 <= record$true_hidden & record$b.high2 >= record$true_hidden), TRUE, FALSE)
record$separate <- ifelse((record$freq <= record$true_freq-0.15 | record$freq >= record$true_freq+0.15), "wrong", record$separate)
record$capture1.bt <- ifelse((record$bt.null.low <= record$true_hidden & record$bt.null.high >= record$true_hidden), TRUE, FALSE)
record$capture2.bt <- ifelse((record$bt.alt.low <= record$true_hidden & record$bt.alt.high >= record$true_hidden), TRUE, FALSE)
#record$reject <- ifelse((record$b.low2 <= 0 & record$b.high2 >= 0), FALSE, TRUE)



if (sum(true_hidden) == 0) {
  print("Null case")
  ######### For null
  #### CI1
  
  record$position <- runif(nrow(record), min = -0.5, max = 0.5)
  
  p1 <- ggplot(record, aes(y=b_hat.null, x=position))  +
    geom_errorbar(aes(ymin=b.low, ymax=b.high),color = "grey") +
    geom_point(aes(colour = factor(capture1))) +
    #geom_point() +
    #geom_jitter(aes(colour = factor(capture))) +
    ylim(-0.05, 0.05) +
    geom_hline(yintercept=0, color = "red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p1
  
  ##### CI2
  p3 <- ggplot(record, aes(y=b_hat.alt, x=position))  +
    geom_errorbar(aes(ymin=b.low2, ymax=b.high2),color = "grey") +
    geom_point(aes(colour = factor(capture2))) +
    #geom_point() +
    #geom_jitter(aes(colour = factor(capture))) +
    ylim(-0.05, 0.05) +
    geom_hline(yintercept=0, color = "red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p3
  
  ##### CI bt null
  p4 <- ggplot(record, aes(y=b_hat.null, x=position))  +
    geom_errorbar(aes(ymin=bt.null.low, ymax=bt.null.high),color = "grey") +
    geom_point(aes(colour = factor(capture1.bt))) +
    #geom_point() +
    #geom_jitter(aes(colour = factor(capture))) +
    ylim(-0.05, 0.05) +
    geom_hline(yintercept=0, color = "red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p4
  ##### CI bt alt
  p5 <- ggplot(record, aes(y=b_hat.alt, x=position))  +
    geom_errorbar(aes(ymin=bt.alt.low, ymax=bt.alt.high),color = "grey") +
    geom_point(aes(colour = factor(capture2.bt))) +
    #geom_point() +
    #geom_jitter(aes(colour = factor(capture))) +
    ylim(-0.05, 0.05) +
    geom_hline(yintercept=0, color = "red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p5
  
  # ######
  ## frequency plot
  p2 <- ggplot(record, aes(y=freq, x=true_freq))  +
    geom_point(aes(colour = factor(separate))) +
    #geom_point() +
    #xlim(-0.25, 0.8) + ylim(-0.25, 0.8) +
    geom_abline(intercept = 0, slope = 1, color="red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p2
  
  
  
  
  
  p1.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/null.emMedian.ci1.Freq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".pdf", sep = "")
  p3.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/null.emMedian.ci2.Freq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".pdf", sep = "")
  p2.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/null.emMedian.estimateFreq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".pdf", sep = "")
  p4.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/bt.null.emMedian.ci1.Freq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".pdf", sep = "")
  p5.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/bt.null.emMedian.ci2.Freq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".pdf", sep = "")
  
  csv.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/null.emMedian.Freq-", known_freq[1], "-", known_freq[2],
                        "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".csv", sep = "")
} else {
  print("Alt case")
  ############ For beta != 0
  ###### CI 1
  p1 <- ggplot(record, aes(y=b_hat.null, x=true_hidden))  +
    geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") +
    geom_point(aes(colour = factor(capture1)) ) +
    #geom_point() +
    #geom_jitter(aes(colour = factor(capture))) +
    xlim(-0.1, 0.8) + ylim(-0.1, 0.8) +
    geom_abline(intercept = 0, slope = 1, color="red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p1
  
  ###### CI 2
  p3 <- ggplot(record, aes(y=b_hat.alt, x=true_hidden))  +
    geom_errorbar(aes(ymin=b.low2, ymax=b.high2), color = "grey") +
    #geom_point(aes(colour = factor(reject)) ) +
    geom_point(aes(colour = factor(capture2)) ) +
    #geom_point() +
    #geom_jitter(aes(colour = factor(capture2))) +
    xlim(-0.1, 0.8) + ylim(-0.1, 0.8) +
    geom_abline(intercept = 0, slope = 1, color="red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p3
  # ######
  ## frequency plot
  p2 <- ggplot(record, aes(y=freq, x=true_freq))  +
    geom_point(aes(colour = factor(separate))) +
    #geom_point() +
    #xlim(-0.25, 0.8) + ylim(-0.25, 0.8) +
    geom_abline(intercept = 0, slope = 1, color="red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p2
  
  ##### CI bt null
  p4 <- ggplot(record, aes(y=b_hat.null, x=true_hidden))  +
    geom_errorbar(aes(ymin=bt.null.low, ymax=bt.null.high),color = "grey") +
    geom_point(aes(colour = factor(capture1.bt))) +
    #geom_point() +
    #geom_jitter(aes(colour = factor(capture))) +
    xlim(-0.1, 0.8) + ylim(-0.1, 0.8) +
    geom_abline(intercept = 0, slope = 1, color="red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p4
  ##### CI bt alt
  p5 <- ggplot(record, aes(y=b_hat.alt, x=true_hidden))  +
    geom_errorbar(aes(ymin=bt.alt.low, ymax=bt.alt.high),color = "grey") +
    geom_point(aes(colour = factor(capture2.bt))) +
    #geom_point() +
    #geom_jitter(aes(colour = factor(capture))) +
    xlim(-0.1, 0.8) + ylim(-0.1, 0.8) +
    geom_abline(intercept = 0, slope = 1, color="red") +
    ggtitle(paste("Freq=", known_freq[1], ",", known_freq[2], ",", known_freq[3], ", b1=", known_beta[1], ", b2=", known_beta[2], sep = ''))
  p5
  
  p1.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/emMedian.ci1.Freq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".pdf", sep = "")
  p3.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/emMedian.ci2.Freq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".pdf", sep = "")
  p2.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/emMedian.estimateFreq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".pdf", sep = "")
  p4.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/bt.emMedian.ci1.Freq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".pdf", sep = "")
  p5.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/bt.emMedian.ci2.Freq-", known_freq[1], "-", known_freq[2],
                       "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".pdf", sep = "")
  csv.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/emMedian.Freq-", known_freq[1], "-", known_freq[2],
                        "-", known_freq[3], "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-",num_snp,".numSet-",set, ".csv", sep = "")
}



# # output to files
#### CI 1
# pdf(p1.filename)
# p1
# dev.off()
# ###### CI 2
# pdf(p3.filename)
# p3
# dev.off()
# ###### CI bt 1
# pdf(p4.filename)
# p4
# dev.off()
# ###### CI bt 2
# pdf(p5.filename)
# p5
# dev.off()
# # frequency plots
# pdf(p2.filename)
# p2
# dev.off()
# ###### write the dataframe
# write.csv(record, csv.filename, row.names = FALSE)



