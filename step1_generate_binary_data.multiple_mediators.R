library("dplyr")
library("ggplot2")
library("reshape2")
library(ggpubr)
library(fitdistrplus)
library(SimDesign)
library(mixtools)
library(fitdistrplus)




setwd("/Users/dr/Applications/MR_predictor/mr_predictor/mytest11-muti_snp_dis_test")
#setwd("/home/dingzh/t2d_mediation/mr_predictor_operator")
folder <- "/Users/dr/Desktop/mediation_summary/data_sets/binary/five_known/hidden_freq_0.8"


num_snp <- 70
max_num_mediators <- 6
num_known_mediators <- 5
num_sets <- 50
chances <- c(0.5, 0.6, 0.8, 0.2, 0.5)  # probs of the knowns
h_chance <- 0.7          # prob of the hidden
known_beta <- c(0.4, 0.2, 0.3, 0.2, 0.4)           # beta1 and beta2

# true beta values
true_hidden <- runif(num_sets, min = 0.02, max = 0.5)

# true beta = 0
#true_hidden <- rep(0, num_sets)

snp_col <- c()
for (i in 1:num_snp) {
  snp_col[i] <- paste("snp", i, sep="")
}

mediator_col <- c()
for (i in 1:max_num_mediators) {
  mediator_col[i] <- paste("m", i, sep="")
}


####### create sub-directories
create_dir_command <- paste("mkdir ", folder, "/set{", do.call(paste, c(as.list(paste(1:num_sets, sep = "", collpase = "")), sep = ",")), "}", sep = "")
system(command = create_dir_command)

# mr_predictor loop
for (f in 1:num_sets) {
  tryCatch({
    
    print(paste("--------------------This is loop:", f, "----------------------------"))
    
    setwd(paste(folder, "/set", f,sep = ""))
    
    sign <- sample(c(-1,1), size = num_snp, replace = TRUE)
    
    #### generate iifile.txt ####
    ### no correlation for now
    iifile <- data.frame()
    write.table(iifile, "iifile.txt", sep = " ", col.names = FALSE,
                row.names = FALSE, quote = FALSE)
    
    ### have correlation
    # iifile <- data.frame("trait1" = c("m1", "m3"),
    #                      "trait2" = c("m2", "m4"),
    #                      "cov" = c(0.2, 0.3))
    # write.table(iifile, "iifile.txt", sep = " ", col.names = FALSE,
    #             row.names = FALSE, quote = FALSE)
    
    #### generate_infosheet ####
    #norm <- rnorm(num_snp, mean = 0.30447, sd = 0.1018345)
    #norm <- rnorm(num_snp, mean = 0.150447, sd = 0.03018345)
    norm <- runif(num_snp, min = 0.1, max = 0.4)
    #norm <- runif(num_snp, min = 0.01, max = 0.05)
    #norm <- rnorm(num_snp, mean = 0.0100, sd = 0.0030)
    
    
    norm[norm<0] <- 0.01
    norm[norm>1] <- 0.99
    
    #### A is positive, T is negative
    alt_nucleotide <- ifelse(sign == 1, "A", "T")
    ref_nucleotide <- ifelse(sign == 1, "T", "A")
    
    
    infosheet <- data.frame("snp" = "snp", "alt" = "A", "ref" = "T", 
                            "freq" = norm, "index" = c(1:num_snp))
    infosheet$snp <- paste(infosheet$snp, infosheet$index, sep = "")
    infosheet <- subset(infosheet, select = c("snp", "alt", "ref", "freq"))
    
    write.table(infosheet, "infosheet.txt", sep = " ", col.names = FALSE, 
                row.names = FALSE, quote = FALSE)
    
    #### generate_scorefile: SNP -> M ####
    c_prime_present <- vector(mode = "logical", length = 50)
    
    
    effects <- rep(0.2, num_snp)
    sd <- rep(0.001, num_snp)
    # effects <- rnorm(num_snp, mean = 0.011, sd = 0.0015)
    # sd <- rnorm(num_snp, mean = 0.0011, sd = 0.00015)
    
    scoresheets <- c()
    
    ## %%% edit here to change minimun number of mediators %%%
    #num_mediator <- rep(3, num_snp)
    #num_mediator <- sample(2:num_known_mediators, num_snp, replace=T)
    
    
    
    
    combines_scoresheets <- data.frame("snp" = c(), "med" = c(), "effect" = c(), 
                                       "sd" = c(),
                                       "fifth" = c(), "sixth" = c())
    
    num_mediator <- vector(mode = "logical", length = num_snp)
    for (i in 1:num_snp) {
      
      # each mediation now have the same change (50%) to be included: need to go further later^^^
      
      selected_mediators <- c()
      for (l in 1:num_known_mediators) {
        if (sample(0:1, size = 1, replace = TRUE, prob = c(1-chances[l], chances[l])) == 1) {
          selected_mediators <- c(selected_mediators, l)
        }
      }
      
      num_mediator[i] <- length(selected_mediators)
      
      
      #selected_mediators <- sample(1:num_known_mediators, num_mediator[i], replace=F)
      
      # %%% the rate of having m4
      if (sample(0:1, size = 1, replace = TRUE, prob = c(1-h_chance, h_chance)) == 1) {
        #if (1 == 1) {
        selected_mediators <- c(selected_mediators, max_num_mediators)
        num_mediator[i] <- num_mediator[i] + 1
      }
      
      # if current snp is not associated with any mediators, skip current loop
      if (num_mediator[i] == 0) {
        next
      }
      
      
      c_prime_present[i] <- (max_num_mediators %in% selected_mediators)
      
      snp_local <- snp_col[i]
      
      local_mediators <- c()
      for (k in 1:num_mediator[i]) {
        local_mediators[k] <- paste("m", selected_mediators[k], sep="")
      }
      
      
      effect_local <- rnorm(num_mediator[i], mean = effects[i], sd = 0.05)
      
      # change the sign
      effect_local <- effect_local * sign[i]
      sd_local <- rep(sd[i], num_mediator[i])
      scoresheet <- data.frame("snp" = snp_local, "med" = local_mediators, 
                               "effect" = effect_local, "sd" = sd_local,
                               "fifth" = 0.00, "sixth" = -9)
      combines_scoresheets <- rbind(combines_scoresheets, scoresheet)
    }
    
    # assume a all from the same distritbiion 
    combines_scoresheets$effect <- rnorm(nrow(combines_scoresheets), mean = 0.11, sd = 0.015)
    
    write.table(combines_scoresheets, "scorefile.txt", sep = " ", 
                col.names = FALSE, row.names = FALSE, quote = FALSE)
    
    
    #### generate_phenofile ####
    phenofile1 <- data.frame("name" = c("T2D"), "type" = c("d"), 
                             "kind" = c("b"), "freq" = c("0.094"))
    
    phenofile2 <- data.frame("name" = mediator_col, "type" = "i", 
                             "kind" = "-9", "freq" = "-9")
    phenofile3 <- data.frame("name" = c("SEX", "income"), "type" = c("c", "c"), 
                             "kind" = c("b", "b"), "freq" = c("0.5", "0.2"))
    
    phenofile <- rbind(phenofile1, phenofile2, phenofile3)
    write.table(phenofile, "phenofile.txt", sep = " ", col.names = FALSE, 
                row.names = FALSE, quote = FALSE)
    
    #### generate_idfile: M -> T2D ####
    
    #effects_id <- runif(max_num_mediators-1, min=0.3, max=0.7)
    effects_id <- known_beta
    sd_id <- rep(0.001, max_num_mediators-1)
    
    
    
    idfile1 <- data.frame("T2D" = "T2D", "med" = mediator_col[1:(length(mediator_col)-1)], "eff" = effects_id,
                          "sd" = sd_id)
    
    
    ## %%% Edit here to set the size of unknown effect %%%
    idfile2 <- data.frame("T2D" = c("T2D", "T2D", "T2D"), "med" = c("m6", "SEX", "income"),
                          "eff" = c(true_hidden[f], "-0.01", "0.01"), "sd" = c("0.001", "0.001", "0.001"))
    # idfile2 <- data.frame("T2D" = c("T2D", "T2D", "T2D"), "med" = c("m2", "SEX", "income"), 
    #                       "eff" = c(true_hidden[f], "-0.01", "0.01"), "sd" = c("0.001", "0.001", "0.001"))
    
    idfile <- rbind(idfile1, idfile2)
    write.table(idfile, "idfile.txt", sep = " ", col.names = FALSE, 
                row.names = FALSE, quote = FALSE)
    
    
    
    #### mr_predictor ####
    
    out_dir <- paste(folder, "/set", f, sep = "")
    #system("perl /appl/mrpredictor-0.028/mr_predictor.pl --score scorefile.txt --info infosheet.txt --pheno phenofile.txt --ii-rel iifile.txt --id-rel idfile.txt --scores --pheno-verbose --nsims 1 --cc_nsamp T2D 10000 10000 --out test11 --skip-plink")
    my_predictor_command <- paste("/Users/dr/Applications/MR_predictor/mr_predictor/mr_predictor.pl --score scorefile.txt --info infosheet.txt --pheno phenofile.txt --ii-rel iifile.txt --id-rel idfile.txt --scores --pheno-verbose --nsims 1 --cc_nsamp T2D 5000 10000 --out ",
                                  out_dir, "/data", " --skip-plink", sep = "")
    
    system(command = my_predictor_command)
    
    #### import from mr_predictor ####
    pheno <- read.delim("data_1.pheno", sep = " ", dec =".", header = TRUE)
    cov <- read.delim("data_1.cov", sep = " ", dec = ".", header = TRUE)
    ped <- read.delim("data_1.ped", sep = " ", dec = ".", header = TRUE)
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
    
  }, error=function(e){print("Outter caught")})
}


setwd(folder)
df.true <- data.frame("true_hidden" = true_hidden)
write.csv(df.true, "true_hidden.csv", row.names = FALSE)



