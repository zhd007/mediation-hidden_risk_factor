library("dplyr")
library("ggplot2")
library("reshape2")
library(ggpubr)
library(fitdistrplus)
library(SimDesign)
library(mixtools)
library(fitdistrplus)
library(msm)
library(Matrix)
library(magic)

folder <- "/Users/dr/Desktop/mediation_summary/data_sets/binary/five_known/b006.null"


num_snp <- 70
max_num_mediators <- 6
num_known_mediators <- 5
num_sets <-100
chances <- c(0.5, 0.6, 0.8, 0.2, 0.5)   # probs of the knowns
h_chance <- 1          # prob of the hidden
known_beta <- c(0.4, 0.2, 0.3, 0.2, 0.4)           # beta1 and beta2

em.times <- 51

# true beta values
true_hidden <- read.csv(paste(folder, "/true_hidden.csv", sep = ""), stringsAsFactors = FALSE)

snp_col <- c()
for (i in 1:num_snp) {
  snp_col[i] <- paste("snp", i, sep="")
}

mediator_col <- c()
for (i in 1:max_num_mediators) {
  mediator_col[i] <- paste("m", i, sep="")
}


record <- data.frame("freq" = -9, "b_hat" = -9, "est_a" = -9, "true_freq" = -9, 
                     "b_hat.sub" = -9, 
                     # "b_hat.alt1" = -9,
                     # "b_hat.alt2" = -9,
                     # "b_hat.alt3" = -9,
                     # "b_hat.alt4" = -9,
                     # "b_hat.alt5" = -9,
                     # "b_hat.alt6" = -9,
                     # "b_hat.alt7" = -9,
                     # "b_hat.alt8" = -9,
                     "b_hat.alt9" = -9,
                     "b_hat.alt10" = -9,
                     "b.low" = -9, "b.high" = -9,
                     "b.low2" = -9, "b.high2" = -9,
                     "b.low.uni" = -9, "b.high.uni" = -9,
                     "true_hidden" = -9)

# mr_predictor loop
for (f in 1:num_sets) {
#for (f in 1:50) {
  tryCatch({
    
    print(paste("--------------------This is loop:", f, "----------------------------"))
    setwd(paste(folder, "/set", f,sep = ""))
    
    #### import from mr_predictor ####
    if (num_snp == 70) {
      # 70
      pheno <- read.delim("test1_1.pheno", sep = " ", dec =".", header = TRUE)
      cov <- read.delim("test1_1.cov", sep = " ", dec = ".", header = TRUE)
      ped <- read.delim("test1_1.ped", sep = " ", dec = ".", header = TRUE)
    } else {
      # 500
      pheno <- read.delim("test1_1.pheno", sep = " ", dec =".", header = TRUE)
      cov <- read.delim("test1_1.cov", sep = " ", dec = ".", header = TRUE)
      ped <- read.delim("test1_1.ped", sep = " ", dec = ".", header = TRUE)
    }

    ped.name <- c("FID", "IID", "PID", "MID", "SEX", "T2D", "None")
    
    for (i in 1:(num_snp)) {
      var1 <- paste("g", i, "a", sep = "")
      var2 <- paste("g", i, "b", sep = "")
      ped.name <- c(ped.name, var1, var2)
    }
    colnames(ped) <- ped.name
    
    merged <- merge(pheno, cov[, c("FID", "SEX", "income")], 
                    by, by.x = "FID", by.y = "FID")
    merged <- merge(merged, ped[, c(1,8:ncol(ped))], 
                    by, by.x = "FID", by.y = "FID")
    
    
    #change values
    merged$T2D[merged$T2D == 1] <- 0
    merged$T2D[merged$T2D == 2] <- 1
    merged$SEX[merged$SEX == 1] <- 0
    merged$SEX[merged$SEX == 2] <- 1
    
    for (i in 1:num_snp) {
      var = paste("snp", i, sep = "")
      varGene1 = paste("g", i, "a", sep = "")
      varGene2 = paste("g", i, "b", sep = "")
      merged[[var]] <- 0
      merged[[var]][merged[[varGene1]] == "A" | merged[[varGene2]] == "A"] <- 1 
    }
    
    write.csv(merged, "merged.csv", row.names = FALSE)
    
    
    ####### import merged.csv #############
    merged <- read.csv("merged.csv", stringsAsFactors = FALSE)
    
    ######## old model: mediation like model ########
    method <- "mediation-like"
    original_merged <- merged

    result <- list(type = any)
    # this is for p-value
    result_var <- list(type = any)
    result_r <- list(type = any)

    for (i in 1:num_snp) {

      #snp~mediators
      for (j in 1:max_num_mediators) {
        md_eq <- paste(mediator_col[j], "~", snp_col[i], sep = "")
        assign(mediator_col[j], lm(md_eq, data = original_merged))
      }

      #mediators~T2D
      t_eq1 <- paste("T2D~", snp_col[i], "+", sep = "")
      for (k in 1:num_known_mediators) {
        t_eq1 <- paste(t_eq1, mediator_col[k], "+", sep = "")
      }
      #t_eq1 <- paste(t_eq1, paste(snp_col[snp_col != snp_col[i]], collapse = "+"), sep = "")
      t_eq1 <- paste(t_eq1, "SEX+income", sep = "")
      #t_eq1 <- paste(t_eq1, "+0", sep = "")            # no intersept

      t_eq2 <- paste("T2D~", snp_col[i], "+", sep = "")
      for (k in 1:max_num_mediators) {
        t_eq2 <- paste(t_eq2, mediator_col[k], "+", sep = "")
      }
      t_eq2 <- paste(t_eq2, "SEX+income", sep = "")
      #t_eq2 <- paste(t_eq2, "+0", sep = "")            # no intersept


      total1 <- glm(t_eq1, data = original_merged, family = binomial("logit"))
      #total1 <- lm(t_eq1, data = original_merged)
      total2 <- glm(t_eq2, data = original_merged, family = binomial("logit"))

      #get coefficients
      local_coe <- c()
      local_var <- c()
      for (k in 1:num_known_mediators) {

        #get a
        local_coe[2*k-1] <- get(paste('m', k, sep=""))$coefficient[2]
        local_var[2*k-1] <- summary(get(paste('m', k, sep="")))$coefficient[2,2]^2

        #get b
        local_coe[2*k] <- total1$coefficient[k+2]
        local_var[2*k] <- summary(total1)$coefficient[k+2,2]^2
      }

      #get c'
      local_coe[2*num_known_mediators+1] <- total1$coefficient[2]
      local_var[2*num_known_mediators+1] <- summary(total1)$coefficient[2,2]^2

      result[[i]] <- local_coe
      result_var[[i]] <- local_var

      #get real coefficients
      local_coe_r <- c()
      for (k in 1:max_num_mediators) {

        #get real a
        #if (summary(get(paste('m', k, sep="")))$coefficient[2, 4] < 0.05){
        #  local_coe_r[2*k-1] <- get(paste('m', k, sep=""))$coefficient[2]
        #} else {
        #  local_coe_r[2*k-1] <- 0
        #}
        local_coe_r[2*k-1] <- get(paste('m', k, sep=""))$coefficient[2]

        #get real b
        #if (summary(total2)$coefficient[k+2, 4] < 0.05) {
        #  local_coe_r[2*k] <- total2$coefficient[k+2]
        #} else {
        #  local_coe_r[2*k] <- 0
        #}
        local_coe_r[2*k] <- total2$coefficient[k+2]
      }

      #get real c'
      #if (summary(total2)$coefficient[2, 4] < 0.05) {
      #  local_coe_r[2*max_num_mediators+1] <- total2$coefficient[2]
      #} else {
      #  local_coe_r[2*max_num_mediators+1] <- 0
      #}
      local_coe_r[2*max_num_mediators+1] <- total2$coefficient[2]

      result_r[[i]] <- local_coe_r
    }
    result <- data.frame(t(sapply(result,c)))
    result_var <- data.frame(t(sapply(result_var,c)))

    result_r <- data.frame(t(sapply(result_r,c)))

    # add header for result
    header <- c()
    for (i in 1:num_known_mediators) {
      header[length(header) + 1] <- paste("a", i, sep = "")
      header[length(header) + 1] <- paste("b", i, sep = "")
    }
    header[length(header) + 1] <- "c_prime"

    colnames(result) <- header
    rownames(result) <- c()
    #hist(result$c_prime)

    # add header for result_var
    header <- c()
    for (i in 1:num_known_mediators) {
      header[length(header) + 1] <- paste("a", i, ".var", sep = "")
      header[length(header) + 1] <- paste("b", i, ".var", sep = "")
    }
    header[length(header) + 1] <- "c_prime.var"

    colnames(result_var) <- header
    rownames(result_var) <- c()


    # add header for result_r
    header_r <- c()
    for (i in 1:max_num_mediators) {
      header_r[length(header_r) + 1] <- paste("a", i, sep = "")
      header_r[length(header_r) + 1] <- paste("b", i, sep = "")
    }
    header_r[length(header_r) + 1] <- "c_prime"

    colnames(result_r) <- header_r
    rownames(result_r) <- c()

    result$ab1 <- result$a1 * result$b1
    result$ab2 <- result$a2 * result$b2

    # attach variance
    result <- merge(result, result_var, by = 0)

    result_r$ab1 <- result_r$a1 * result_r$b1
    result_r$ab2 <- result_r$a2 * result_r$b2
    result_r$ab3 <- result_r$a3 * result_r$b3
    
    #result$c_prime_present <- c_prime_present
    
    ######## new model: mediation model ########
    # method <- "mediation"
    # df.snp <- merged
    # df.snp$y <- df.snp$T2D
    # #### mediation point estimate
    # # snp -> mediators
    # # alternative way
    # # 50 snps
    # # eq.m1.all <- "m1~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50"
    # # eq.m2.all <- "m2~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50"
    # # eq.m3.all <- "m3~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50"
    # 
    # 
    # # 70 snps
    # # eq.m1.all <- "m1~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+snp51+snp52+snp53+snp54+snp55+snp56+snp57+snp58+snp59+snp60+snp61+snp62+snp63+snp64+snp65+snp66+snp67+snp68+snp69+snp70"
    # # eq.m2.all <- "m2~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+snp51+snp52+snp53+snp54+snp55+snp56+snp57+snp58+snp59+snp60+snp61+snp62+snp63+snp64+snp65+snp66+snp67+snp68+snp69+snp70"
    # # eq.m3.all <- "m3~snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+snp51+snp52+snp53+snp54+snp55+snp56+snp57+snp58+snp59+snp60+snp61+snp62+snp63+snp64+snp65+snp66+snp67+snp68+snp69+snp70"
    # 
    # # 500 snps
    # snp_sum <- paste("snp", 1:num_snp, sep = "", collpase = "")
    # snp_sum.cat <- do.call(paste, c(as.list(snp_sum), sep = "+"))
    # eq.m1.all <- paste("m1~", snp_sum.cat, sep = "")
    # eq.m2.all <- paste("m2~", snp_sum.cat, sep = "")
    # eq.m3.all <- paste("m3~", snp_sum.cat, sep = "")
    # eq.m4.all <- paste("m4~", snp_sum.cat, sep = "")
    # eq.m5.all <- paste("m5~", snp_sum.cat, sep = "")
    # eq.m6.all <- paste("m6~", snp_sum.cat, sep = "")
    # 
    # 
    # m1_coefs <- summary(lm(eq.m1.all, data = df.snp))$coefficient
    # m2_coefs <- summary(lm(eq.m2.all, data = df.snp))$coefficient
    # m3_coefs <- summary(lm(eq.m3.all, data = df.snp))$coefficient
    # m4_coefs <- summary(lm(eq.m4.all, data = df.snp))$coefficient
    # m5_coefs <- summary(lm(eq.m5.all, data = df.snp))$coefficient
    # m6_coefs <- summary(lm(eq.m6.all, data = df.snp))$coefficient
    # 
    # a1.alt <- m1_coefs[2:(num_snp + 1),1]
    # a2.alt <- m2_coefs[2:(num_snp + 1),1]
    # a3.alt <- m3_coefs[2:(num_snp + 1),1]
    # a4.alt <- m4_coefs[2:(num_snp + 1),1]
    # a5.alt <- m5_coefs[2:(num_snp + 1),1]
    # a6.alt <- m6_coefs[2:(num_snp + 1),1]
    # 
    # a1.var.alt <- m1_coefs[2:(num_snp + 1),2]^2
    # a2.var.alt <- m2_coefs[2:(num_snp + 1),2]^2
    # a3.var.alt <- m3_coefs[2:(num_snp + 1),2]^2
    # a4.var.alt <- m4_coefs[2:(num_snp + 1),2]^2
    # a5.var.alt <- m5_coefs[2:(num_snp + 1),2]^2
    # a6.var.alt <- m6_coefs[2:(num_snp + 1),2]^2
    # 
    # # mediators -> outcome
    # # alternative way
    # # 50 snp
    # #eq.y.all <- "y~m1+m2 + snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+SEX + income"
    # 
    # # 70 snp
    # #eq.y.all <- "y~m1+m2 + snp1+snp2+snp3+snp4+snp5+snp6+snp7+snp8+snp9+snp10+snp11+snp12+snp13+snp14+snp15+snp16+snp17+snp18+snp19+snp20+snp21+snp22+snp23+snp24+snp25+snp26+snp27+snp28+snp29+snp30+snp31+snp32+snp33+snp34+snp35+snp36+snp37+snp38+snp39+snp40+snp41+snp42+snp43+snp44+snp45+snp46+snp47+snp48+snp49+snp50+snp51+snp52+snp53+snp54+snp55+snp56+snp57+snp58+snp59+snp60+snp61+snp62+snp63+snp64+snp65+snp66+snp67+snp68+snp69+snp70+SEX + income"
    # 
    # # 500 snps
    # eq.y.all <- paste("y~m1+m2+m3+m4+m5+", snp_sum.cat, "+SEX + income", sep = "")
    # 
    # 
    # y_coefs <- summary(glm(eq.y.all, data = df.snp, family = binomial("logit")))$coefficient
    # #y_coefs <- summary(lm(eq.y.all, data = df.snp))$coefficient
    # b1.alt <- y_coefs[2,1]
    # b2.alt <- y_coefs[3,1]
    # b3.alt <- y_coefs[4,1]
    # b4.alt <- y_coefs[5,1]
    # b5.alt <- y_coefs[6,1]
    # 
    # b1.var.alt <- y_coefs[2,2]^2
    # b2.var.alt <- y_coefs[3,2]^2
    # b3.var.alt <- y_coefs[4,2]^2
    # b4.var.alt <- y_coefs[5,2]^2
    # b5.var.alt <- y_coefs[6,2]^2
    # 
    # c.alt <- y_coefs[7:(num_snp + 6),1]
    # c.var.alt <- y_coefs[7:(num_snp + 6),4]^2
    # 
    # result <- data.frame("a1" = a1.alt, "a2" = a2.alt, "a3" = a3.alt, "a4" = a4.alt, "a5" = a5.alt, "a6" = a6.alt,
    #                      "a1.var" = a1.var.alt, "a2.var" = a2.var.alt, "a3.var" = a3.var.alt,
    #                      "a4.var" = a4.var.alt, "a5.var" = a5.var.alt, "a6.var" = a6.var.alt,
    #                      "b1" = b1.alt, "b2" = b2.alt, "b3" = b3.alt, "b4" = b4.alt, "b5" = b5.alt,
    #                      "c_prime" = c.alt,
    #                      "b1.var" = b1.var.alt, "b2.var" = b2.var.alt, "b3.var" = b3.var.alt,
    #                      "b4.var" = b4.var.alt, "b5.var" = b5.var.alt,
    #                      "c_prime.var" = c.var.alt)


    ######## estimate ########

    # est1 <- vector(mode = "logical", length = 5000)
    # est2 <- vector(mode = "logical", length = 5000)
    # est3 <- vector(mode = "logical", length = 5000)
    # est4 <- vector(mode = "logical", length = 5000)
    # est5 <- vector(mode = "logical", length = 5000)
    # est6 <- vector(mode = "logical", length = 5000)
    # est7 <- vector(mode = "logical", length = 5000)
    # est8 <- vector(mode = "logical", length = 5000)
    # est9 <- vector(mode = "logical", length = 5000)
    # est10 <- vector(mode = "logical", length = 5000)

    # option (1)
    #fit a single normal distribution to a1 and a2

    # a_dist <- fitdist(c(result$a1, result$a2), distr = "norm", method = c("mle"),
    #                   start=NULL, fix.arg=NULL, keepdata = TRUE, keepdata.nb=100)
    # #plot(a_dist)
    # for (u in 1:5000) {
    #   # sample 50 a from the a distribution
    #   a <- rnorm(n = nrow(result), mean = a_dist$estimate[1], sd = a_dist$estimate[2])
    #   b.sort <- sort(result$c_prime) / sort(a)
    #   b <- result$c_prime / a
    #
    #   b.sub <- sort(b)[5:46]
    #   b.sort.sub <- sort(b.sort)[5:46]
    #
    #   est1[u] <- mean(b.sub)
    #   est2[u] <- median(b.sub)
    #   est3[u] <- mean(b.sort.sub)
    #   est4[u] <- median(b.sort.sub)
    # }
    # em_median1 <- mean(est1)
    # em_median2 <- mean(est2)
    # em_median3 <- mean(est3)
    # em_median4 <- mean(est4)


    ######### restricted EM
    # option (2)
    # try using em mixtool
    # iterator <- 0
    # pivot <- 10000
    # pivot2 <- 0
    # lambda1 <- 0
    # lambda2 <- 0
    # while ((pivot >= 0.05 | pivot <= -0.05) | (pivot2 < 0.05)  |   (lambda1 < 0.05 | lambda2 < 0.05)) {
    #   if (iterator > 300) {
    #     break
    #   }
    #   tryCatch({
    #     iterator <- iterator + 1
    #     print(paste("Iteration", iterator))
    #
    #     # use em to fit snp to a bimodel normal distribution
    #     snp_em <- quiet(normalmixEM(c(result$a1, result$a2), lambda = .5, mu = c(0, 0.5), arbvar = TRUE), messages = TRUE)
    #     # plot(snp_em, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
    #     #      main2="SNP", xlab2="Beta")
    #
    #     snp_em_summary <- data.frame("mean" = snp_em$mu, "sd" = snp_em$sigma, "lambda" = snp_em$lambda)
    #     pivot <- snp_em_summary$mean[1]
    #     pivot2 <- snp_em_summary$mean[2]
    #     lambda1 <- snp_em_summary$lambda[1]
    #     lambda2 <- snp_em_summary$lambda[2]
    #
    #   }, error=function(e){print("caught")})
    # }
    
    ########### no restriction: EM Median
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
      snp_em.status <- -999
      while (snp_em.status == -999) {
        tryCatch({
          snp_em.a <- quiet(normalmixEM(c(result$a1, result$a2, result$a3, result$a4, result$a5), lambda = .5, mu = c(0, 0.5), arbvar = TRUE), messages = TRUE)
          snp_em.status <- 1
        }, error=function(e){
          print("caught")
          snp_em.status <- -999})
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
    
    
    # ####### for c
    # emset.c <- 21
    # emLR.c <- vector(mode = "logical", length = emset.c)
    # emMu1.c <- vector(mode = "logical", length = emset.c)
    # emMu2.c <- vector(mode = "logical", length = emset.c)
    # emSd1.c <- vector(mode = "logical", length = emset.c)
    # emSd2.c <- vector(mode = "logical", length = emset.c)
    # emLambda1.c <- vector(mode = "logical", length = emset.c)
    # emLambda2.c <- vector(mode = "logical", length = emset.c)
    # emModel.c <- vector(mode = "list", length = emset.c)
    # 
    # 
    # for (em.c in 1:emset.c) { 
    #   snp_em.status.c <- -999
    #   while (snp_em.status.c == -999) {
    #     tryCatch({
    #       snp_em.c <- quiet(normalmixEM(result$c_prime, lambda = .5, arbvar = TRUE), messages = TRUE)
    #       snp_em.status.c <- 1
    #     }, error=function(e){
    #       print("caught")
    #       snp_em.status.c <- -999})
    #   }
    #   # plot(snp_em.c, density=TRUE, cex.axis=1.4, cex.lab=1.4, cex.main=1.8,
    #   #      main2="SNP", xlab2="Beta")
    #   emLR.c[em.c] <- snp_em.c$loglik
    #   emMu1.c[em.c] <- snp_em.c$mu[1]
    #   emSd1.c[em.c] <- snp_em.c$sigma[1]
    #   emLambda1.c[em.c] <- snp_em.c$lambda[1]
    #   emMu2.c[em.c] <- snp_em.c$mu[2]
    #   emSd2.c[em.c] <- snp_em.c$sigma[2]
    #   emLambda2.c[em.c] <- snp_em.c$lambda[2]
    #   emModel.c[[em.c]] <- snp_em.c
    # }
    # #maxI.c <- which(emLR.c == max(emLR.c))[1]
    # maxI.c <- which(round(emMu2.c, digits = 6) == round(median(emMu2.c), digits = 6))[1]
    # snp_em_summary.c <- data.frame("mean" = c(emMu1.c[maxI.c], emMu2.c[maxI.c]), 
    #                                "sd" = c(emSd1.c[maxI.c], emSd2.c[maxI.c]), 
    #                                "lambda" = c(emLambda1.c[maxI.c], emLambda2.c[maxI.c]) )
    # medianModel.c <- emModel.c[[maxI.c]]
    
    ########### fit a single normal
    # snp_single_norm <- fitdist(c(result$a1, result$a2, result$a3, result$a4, result$a5), distr = "norm", method = c("mle"),
    #                            start=NULL, fix.arg=NULL, keepdata = TRUE, keepdata.nb=100)
    
    
    # for (u in 1:5000) {
    #   # sample 50 a from the right a distribution
    #   a <- rnorm(n = nrow(result), mean = snp_em_summary$mean[2], sd = snp_em_summary$sd[2])
    #   b.sort <- sort(result$c_prime) / sort(a)
    #   b <- result$c_prime / a
    #   
    #   b.sub <- sort(b)[3:48]
    #   b.sort.sub <- sort(b.sort)[3:48]
    #   
    #   est1[u] <- mean(b.sub)
    #   est2[u] <- median(b.sub)
    #   est3[u] <- mean(b.sort.sub)
    #   est4[u] <- median(b.sort.sub)
    #   
    #   
    #   # sample 50 a from the single normal distribution
    #   a <- rnorm(n = nrow(result), mean = snp_single_norm$estimate[1], sd = snp_single_norm$estimate[2])
    #   b.sort <- sort(result$c_prime) / sort(a)
    #   b <- result$c_prime / a
    #   
    #   b.sub <- sort(b)[3:48]
    #   b.sort.sub <- sort(b.sort)[3:48]
    #   
    #   est5[u] <- mean(b.sub)
    #   est6[u] <- median(b.sub)
    #   est7[u] <- mean(b.sort.sub)
    #   est8[u] <- median(b.sort.sub)
    # }
    # em_median1 <- mean(est1)
    # em_median2 <- mean(est2)
    # em_median3 <- mean(est3)
    # em_median4 <- mean(est4)
    # em_median5 <- mean(est5)
    # em_median6 <- mean(est6)
    # em_median7 <- mean(est7)
    # em_median8 <- mean(est8)
    em_median9 <- mean(result$c_prime) / snp_em_summary$mean[2]       # bimodal
    #em_median9 <- snp_em_summary.c$mean[2] / snp_em_summary$mean[2]       # bimodal
    #em_median10 <- mean(result$c_prime) / snp_single_norm$estimate[1] # unimodal
    
    ###### delta for bimodal 
    # var_c <- var(result$c_prime)
    # var_a1 <- var(result$a1)
    # var_a2 <- var(result$a2)
    # cov_c_a1 <- cov(result$c_prime, result$a1)
    # cov_c_a2 <- cov(result$c_prime, result$a2)
    # cov_a1_a2 <- cov(result$a1, result$a2)
    # 
    # se.delta <- deltamethod(~ x1 / x2,
    #                         mean = c(mean(result$c_prime), snp_em_summary$mean[2]),
    #                         cov = matrix(c(  sd(result$c_prime)/num_snp, 0,
    #                                          0, snp_em_summary$sd[2]/ (num_snp*length(known_beta))  ), nrow = 2, ncol = 2))
    # 
    # b.low <- em_median9 - 1.96 * (se.delta )
    # b.high <- em_median9 + 1.96 * (se.delta )
    
    
    
    ###### delta for unimodal
    # se.delta <- deltamethod(~ x1 / x2,
    #                         mean = c(mean(result$c_prime), snp_single_norm$estimate[1]),
    #                         cov = matrix(c(  sd(result$c_prime)/num_snp, 0,
    #                                          0, snp_single_norm$estimate[2]/ (num_snp*length(known_beta))  ), nrow = 2, ncol = 2))
    # 
    # b.low.uni <- em_median10 - 1.96 * (se.delta )
    # b.high.uni <- em_median10 + 1.96 * (se.delta )
    
    ###### bootstrap just using df.b$c.alt and df.a$a1.alt, for bimodal
    # boot_result <- vector(mode = "logical", length = 5000)
    # for (v in 1:5000) {
    #   boot_a_sample <- rnorm(n = 50, mean = snp_em_summary$mean[2], sd = snp_em_summary$sd[2])
    #   boot_c_sample <- sample(result$c_prime, size  = 50)
    #   boot_result[v] <- mean(boot_c_sample) / mean(boot_a_sample)
    # }
    # b.low2 <- sort(boot_result)[125]
    # b.high2 <- sort(boot_result)[4876]
    
    ###### parametric bootstrap function in mixtool
    #### var(mu.a|a_hat)
    # pboot.a <- quiet(boot.se(medianModel.a, B=1000), messages = TRUE)
    # bp.se.a <- pboot.a$mu.se[2]
    # 
    # 
    # 
    # #### var(a_hat)
    # var_sum.a <- sum(df.a$a1.var.alt, df.a$a2.var.alt) / (   (num_snp*length(known_beta))^2   )
    # var_sum.c <- sum(df.b$c.var.alt) / (   (num_snp)^2   )
    
    #### law of total variance
    # indicator.a <- ifelse(medianModel.a$posterior[,2] > 0.5, yes = 1, no = 0)
    # n2.a <- sum(indicator.a)
    # var_total.a <- snp_em_summary$sd[2]/n2.a
    # var_total.c <- var(result$c_prime) /num_snp
    # 
    # pb.se.delta <- deltamethod(~ x1 / x2,
    #                            mean = c(mean(result$c_prime), snp_em_summary$mean[2]),
    #                            cov = matrix(c(  var_total.c, 0,
    #                                             0, var_total.a  ), nrow = 2, ncol = 2))
    # b.low2 <- em_median9 - 1.96 * (pb.se.delta )
    # b.high2 <- em_median9 + 1.96 * (pb.se.delta )
    
    ##### conditional variance
    p.a <- medianModel.a$posterior[,2] 
    var_1.a <- p.a^2 %*% ( snp_em_summary$sd[2]^2 + 
                             c(result$a1.var, result$a2.var, result$a3.var,
                               result$a4.var, result$a5.var) )  / (sum(p.a)^2) 
    ### total var a
    # l1 <- (num_snp * length(known_beta))
    # x.1 <- paste("x", 1:l1, sep = "", collpase = "")
    # x.1 <- do.call(paste, c(as.list(x.1), sep = "+"))
    # delta.index <- which(p.a > 0.5)
    # x.2 <- paste("x", delta.index, sep = "", collpase = "")
    # x.2 <- do.call(paste, c(as.list(x.2), sep = "+"))
    # delta.eq <- paste("~(", x.2, ")/(", x.1, ")")
    # var_gamma <- deltamethod(as.formula(delta.eq), mean = p.a, cov = diag(p.a * (1-p.a)))
    
    # delta.index1 <- which(p.a <= 0.5)
    # delta.index2 <- which(p.a > 0.5)
    # 
    # m1.a <- mean(p.a[delta.index1])
    # m2.a <- mean(p.a[delta.index2])
    # m1.a <- ifelse(is.na(m1.a), 0, m1.a)
    # m2.a <- ifelse(is.na(m2.a), 0, m2.a)
    # 
    # v1.a <- var(p.a[delta.index1]) /length(delta.index1)
    # v2.a <- var(p.a[delta.index2]) /length(delta.index2)
    # v1.a = ifelse(is.na(v1.a), 0, v1.a)
    # v2.a = ifelse(is.na(v2.a), 0, v2.a)
    # cov.a <- ifelse((is.na(v2.a) |  is.na(v1.a))  , 0, -mean(v1.a, v2.a))
    # se_gamma <- deltamethod(~x2/ (x1 + x2), mean = c(m1.a, m2.a), 
    #                         cov = matrix(c(v1.a, cov.a, cov.a, v2.a), nrow = 2, ncol = 2))
    # 
    # var_2.a <- se_gamma^2 * (snp_em_summary$mean[2]^2 + snp_em_summary$mean[1]^2 -
    #                          2*snp_em_summary$mean[2]*snp_em_summary$mean[1])
    
    # estimate  variance from all EM posterier
    # gamma.a <- vector(mode = "logical", length = em.times )
    # for (bt in 1:em.times ) {
    #   w_model <- emModel.a[[bt]]
    #   w_p <- w_model$posterior[,2] 
    #   delta.index2 <- which(w_p > 0.5)
    #   #w_gamma <- sum(w_p[delta.index2]) / sum(w_p)
    #   #gamma.a[bt] <- w_model$mu[1] * (1-w_gamma) + w_model$mu[2] * w_gamma
    #   gamma.a[bt] <- sum(w_p[delta.index2]) / sum(w_p)
    # }
    # var_2.a <- var(gamma.a) * (snp_em_summary$mean[2]^2 + snp_em_summary$mean[1]^2 -
    #                              2*snp_em_summary$mean[2]*snp_em_summary$mean[1])
    
    delta.index2 <- which(p.a > 0.5)
    gamma.a <- sum(p.a[delta.index2]) / sum(p.a)
    var_2.a <- gamma.a * (1-gamma.a) * (snp_em_summary$mean[2]^2 + snp_em_summary$mean[1]^2 -
                                          2*snp_em_summary$mean[2]*snp_em_summary$mean[1])
    
    total_var.a <- var_1.a + var_2.a
    
    # var of c
    var.c <- sum( var(result$c_prime) + result$c_prime.var )  / (length(result$c_prime)^2)
    se.delta <- deltamethod(~ x1 / x2,
                                 mean = c(mean(result$c_prime), snp_em_summary$mean[2]),
                                 cov = matrix(c(var.c, 0,0, total_var.a), nrow = 2, ncol = 2))
    b.low2 <- em_median9 - 1.96 * se.delta
    b.high2 <- em_median9 + 1.96 * se.delta
    
    
    ###### bootstrap from all_a and all_c
    boot_result <- vector(mode = "logical", length = 1000)
    for (v in 1:1000) {
      if(v %% 100 == 0) {
        print(paste("##### bt round:", v, " #####", sep = ""))
      }
      bt_index <- sample(1:length(result$c_prime), size = length(result$c_prime), replace = TRUE)
      boot_a_sample <- c(result$a1[bt_index], result$a2[bt_index], result$a3[bt_index],
                         result$a4[bt_index], result$a4[bt_index])
      boot_c_sample <- result$c_prime[bt_index]
      ####### for a
      emset.a <- 15
      emLR.a <- vector(mode = "logical", length = emset.a)
      emMu1.a <- vector(mode = "logical", length = emset.a)
      emMu2.a <- vector(mode = "logical", length = emset.a)
      emSd1.a <- vector(mode = "logical", length = emset.a)
      emSd2.a <- vector(mode = "logical", length = emset.a)
      emLambda1.a <- vector(mode = "logical", length = emset.a)
      emLambda2.a <- vector(mode = "logical", length = emset.a)
      boot.emModel.a <- vector(mode = "list", length = emset.a)
      
      for (em.a in 1:emset.a) {
        snp_em.status.a <- -999
        while (snp_em.status.a == -999) {
          tryCatch({
            snp_em.a <- quiet(normalmixEM(boot_a_sample,
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
        boot.emModel.a[[em.a]] <- snp_em.a
      }
      #maxI.a <- which(emLR.a == max(emLR.a))[1]
      maxI.a <- which(round(emMu2.a, digits = 6) == round(median(emMu2.a), digits = 6))[1]
      boot.medianModel.a <- boot.emModel.a[[maxI.a]]
      
      ####### for c
      bt.est.c <- mean(boot_c_sample)
   
      # different results based on unimodal or bimodal
      d.a  <- abs(boot.medianModel.a$mu[1] - boot.medianModel.a$mu[2]) / 
        (2 * boot.medianModel.a$sigma[1] * boot.medianModel.a$sigma[2])
      
      if (d.a <= 1  ) {
        bt.est.a <- mean(boot_a_sample)
      } else {
        bt.est.a <- boot.medianModel.a$mu[2]
      }
      
      boot_result[v] <- bt.est.c / bt.est.a
    }
    bt.low <- sort(boot_result)[25]
    bt.high <- sort(boot_result)[976]
    
    
    ##### default
    b.low <- bt.low
    b.high <- bt.high
    # b.low2 <- -1
    # b.high2 <- -1
    
    
    
    
    
    #### record result
    temp.record <- data.frame("freq" = -1,
                              "b_hat" = -1,
                              "est_a" = -1,
                              "true_freq" = -1,
                              "b_hat.sub" = -1,
                              # "b_hat.alt1" = -1,
                              # "b_hat.alt2" = -1,
                              # "b_hat.alt3" = -1,
                              # "b_hat.alt4" = -1,
                              # "b_hat.alt5" = -1,
                              # "b_hat.alt6" = -1,
                              # "b_hat.alt7" = -1,
                              # "b_hat.alt8" = -1,
                              "b_hat.alt9" = em_median9,
                              "b_hat.alt10" = -1,
                              "b.low" = b.low, "b.high" = b.high,
                              "b.low2" = b.low2, "b.high2" = b.high2,
                              "b.low.uni" = -1, "b.high.uni" = -1,
                              "true_hidden" = true_hidden$true_hidden[f])
    record <- rbind(record,  temp.record)
    print(temp.record)
    
    
  }, error=function(e){print("Outter caught")})
}
record <- record[record$freq != -9,]



#record$true_hidden <- true_hidden
#record$capture <- ifelse((record$b.low <= record$true_hidden & record$b.high >= record$true_hidden), TRUE, FALSE)
# p1 <- ggplot(record, aes(y=b_hat.alt1, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") + 
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) + 
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# p2 <- ggplot(record, aes(y=b_hat.alt2, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") + 
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) + 
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# p3 <- ggplot(record, aes(y=b_hat.alt3, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") + 
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) + 
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# p4 <- ggplot(record, aes(y=b_hat.alt4, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") + 
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) + 
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# ggarrange(p1, p2, p3, p4, 
#           labels = c("1.mean", "1.median", "1.sorted_mean", "1.sorted_median"),
#           ncol = 2, nrow = 2)
# 
# p5 <- ggplot(record, aes(y=b_hat.alt5, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") + 
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) + 
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# p6 <- ggplot(record, aes(y=b_hat.alt6, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") + 
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) + 
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# p7 <- ggplot(record, aes(y=b_hat.alt7, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") + 
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) + 
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# p8 <- ggplot(record, aes(y=b_hat.alt8, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") + 
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) + 
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# 
# ggarrange(p5, p6, p7, p8, 
#           labels = c("2.mean", "2.median", "2.sorted_mean", "2.sorted_median"),
#           ncol = 2, nrow = 2)

record$capture <- ifelse((record$b.low <= record$true_hidden & record$b.high >= record$true_hidden), TRUE, FALSE)
record$capture2 <- ifelse((record$b.low2 <= record$true_hidden & record$b.high2 >= record$true_hidden), TRUE, FALSE)




if (sum(true_hidden) == 0) {
  ##### null case
  # #### for null
  record$position <- runif(nrow(record), min = -0.5, max = 0.5)
  # ##### bimodal ci1
  p9 <- ggplot(record, aes(y=b_hat.alt9, x=position))  +
    geom_errorbar(aes(ymin=b.low2, ymax=b.high2), color = "grey") +
    geom_point(aes(colour = factor(capture2))) +
    #geom_point() +
    ylim(-0.8, 0.8) +
    #ggtitle("test") +
    geom_hline(yintercept=0, color = "red")
  p9

  p9.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/null.", method, ".Freq-", chances[1], "-", chances[2], "-", h_chance,
                       "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-", num_snp, ".numSet-", num_sets,".3_no_mc.pdf", sep = "")
  csv.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/null.", method, ".Freq-", chances[1], "-", chances[2], "-", h_chance,
                        "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-", num_snp, ".numSet-", num_sets,".csv", sep = "")
  
} else {
  ##### alt case
  #### for c != 0
  ##### bimodal ci1
  p9 <- ggplot(record, aes(y=b_hat.alt9, x=true_hidden))  +
    geom_errorbar(aes(ymin=b.low2, ymax=b.high2), color = "grey") +
    geom_point(aes(colour = factor(capture2))) +
    #geom_point() +
    #xlim(-0.4, 0.8) + ylim(-0.4, 0.8) +
    xlim(-0.25, 0.8) + ylim(-0.25, 0.8) +
    #ggtitle("test") +
    geom_abline(intercept = 0, slope = 1, color="red")
  #geom_hline(yintercept=0, color = "red")
  p9
  
  p9.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/", method, ".Freq-", chances[1], "-", chances[2], "-", h_chance,
        "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-", num_snp, ".numSet-", num_sets,".3_no_mc.pdf", sep = "")
  csv.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/", method, ".Freq-", chances[1], "-", chances[2], "-", h_chance,
                        "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-", num_snp, ".numSet-", num_sets,".csv", sep = "")
}



##### bimodal ci2
# p10 <- ggplot(record, aes(y=b_hat.alt9, x=true_hidden))  +
#   geom_errorbar(aes(ymin=b.low2, ymax=b.high2), color = "grey") +
#   geom_point(aes(colour = factor(capture2))) +
#   #geom_point() +
#   #xlim(-0.4, 0.8) + ylim(-0.4, 0.8) +
#   xlim(-0.8, 0.8) + ylim(-0.8, 0.8) +
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# #geom_hline(yintercept=0, color = "red")
# #p10
# ##### unimodal
# p11 <- ggplot(record, aes(y=b_hat.alt10, x=true_hidden))  +
#   geom_errorbar(aes(ymin=b.low.uni, ymax=b.high.uni), color = "grey") +
#   geom_point(aes(colour = factor(capture))) +
#   #geom_point() +
#   #xlim(-0.4, 0.8) + ylim(-0.4, 0.8) +
#   xlim(-0.8, 0.8) + ylim(-0.8, 0.8) +
#   #ggtitle("test") +
#   geom_abline(intercept = 0, slope = 1, color="red")
# #geom_hline(yintercept=0, color = "red")
# #p10



# #### for null
# record$position <- runif(nrow(record), min = -0.5, max = 0.5)
# # ##### bimodal ci1
# p9 <- ggplot(record, aes(y=b_hat.alt9, x=position))  +
#   geom_errorbar(aes(ymin=b.low2, ymax=b.high2), color = "grey") +
#   geom_point(aes(colour = factor(capture2))) +
#   #geom_point() +
#   ylim(-0.4, 0.4) +
#   #ggtitle("test") +
#   geom_hline(yintercept=0, color = "red")
# p9




# ##### bimodal ci2
# p10 <- ggplot(record, aes(y=b_hat.alt9, x=position))  +
#   geom_errorbar(aes(ymin=b.low2, ymax=b.high2), color = "grey") +
#   geom_point(aes(colour = factor(capture2))) +
#   #geom_point() +
#   ylim(-0.4, 0.4) +
#   #ggtitle("test") +
#   geom_hline(yintercept=0, color = "red")
# #p10
# ##### unimodal
# p11 <- ggplot(record, aes(y=b_hat.alt10, x=position))  +
#   geom_errorbar(aes(ymin=b.low2, ymax=b.high2), color = "grey") +
#   geom_point(aes(colour = factor(capture))) +
#   #geom_point() +
#   ylim(-0.4, 0.4) +
#   #ggtitle("test") +
#   geom_hline(yintercept=0, color = "red")
# #p11

# ggarrange(p9, p10, p11, 
#           labels = c("bimodal.ci1","bimodal.ci2", "3.unimodal"),
#           ncol = 2, nrow = 2)
#p9

# pdf(paste("/Users/dr/Desktop/mediation_summary/binary_plots/Freq-", chances[1], "-", chances[2], "-", h_chance,
#           "_b1-0.3_b2-0.6.pdf", sep = ""))
# ggarrange(p1, p2, p3, p4,
#           labels = c("mean", "median", "sorted_mean", "sorted_median"),
#           ncol = 2, nrow = 2)
# dev.off()


# pdf(paste("/Users/dr/Desktop/mediation_summary/binary_plots/Freq-", chances[1], "-", chances[2], "-", h_chance,
#           "_b1-0.3_b2-0.6.pdf", sep = ""))
# p1 <- ggplot(record, aes(y=b_hat.alt1, x=true_hidden))  +
#   #geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") +
#   #geom_point(aes(colour = factor(capture))) +
#   geom_point() +
#   xlim(-0.25, 0.8) + ylim(-0.25, 0.8) +
#   ggtitle(paste("Freq-", chances[1], "-", chances[2], "-", h_chance,
#                 "_b1-0.3_b2-0.6.pdf", sep = "")) +
#   geom_abline(intercept = 0, slope = 1, color="red")
# p1
# dev.off()
# 
# 
# ### plot 1
# pdf(paste("/Users/dr/Desktop/mediation_summary/binary_plots/large_n/m_like.Freq-", chances[1], "-", chances[2], "-", h_chance,
#           "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-", num_snp, ".numSet-", num_sets,".1_bimodal.pdf", sep = ""))
# ggarrange(p1, p2, p3, p4,
#           labels = c("1.mean", "1.median", "1.sorted_mean", "1.sorted_median"),
#           ncol = 2, nrow = 2)
# dev.off()
# 
# #### plot 2
# pdf(paste("/Users/dr/Desktop/mediation_summary/binary_plots/large_n/m_like.Freq-", chances[1], "-", chances[2], "-", h_chance,
#           "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-", num_snp, ".numSet-", num_sets,".2_unimodal.pdf", sep = ""))
# ggarrange(p5, p6, p7, p8,
#           labels = c("2.mean", "2.median", "2.sorted_mean", "2.sorted_median"),
#           ncol = 2, nrow = 2)
# dev.off()




record$capture10 <- ifelse((record$b.low <= record$true_hidden & record$b.high >= record$true_hidden), TRUE, FALSE)
record$position <- runif(nrow(record), min = -0.5, max = 0.5)
p10 <- ggplot(record, aes(y=b_hat.alt9, x=position))  +
  geom_errorbar(aes(ymin=b.low, ymax=b.high), color = "grey") +
  geom_point(aes(colour = factor(capture10))) +
  #geom_point() +
  #xlim(-0.4, 0.8) + ylim(-0.4, 0.8) +
  xlim(-0.25, 0.8) + ylim(-0.25, 0.8) +
  #ggtitle("test") +
  geom_hline(yintercept=0, color = "red")
#geom_hline(yintercept=0, color = "red")
p10

p10.filename <- paste("/Users/dr/Desktop/mediation_summary/simulation_plots/bt.", method, ".Freq-", chances[1], "-", chances[2], "-", h_chance,
                     "_b1-", known_beta[1], "_b2-", known_beta[2], ".numSNP-", num_snp, ".numSet-", num_sets,".3_no_mc.pdf", sep = "")


# pdf("/Users/dr/Desktop/mediation_summary/simulation_plots/test.pdf")
# hist(c(result$a1, result$a2, result$a3, result$a4, result$a5))
# dev.off()

#### plot 3
# pdf(p9.filename )
# # ggarrange(p9, p10, p11,
# #           labels = c("bimodal.ci1","bimodal.ci2", "3.unimodal"),
# #           ncol = 2, nrow = 2)
# p9
# dev.off()
# 
# write.csv(record, csv.filename, row.names = FALSE)


# pdf(p10.filename )
# p10
# dev.off()

