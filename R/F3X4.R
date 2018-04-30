# environment variables
CodonAlphabet <- c("aaa", "aac", "aag", "aat", "aca", "acc", "acg", "act", 
                  "aga", "agc", "agg", "agt", "ata", "atc", "atg", "att", 
                  "caa", "cac", "cag", "cat", "cca", "ccc", "ccg", "cct", 
                  "cga", "cgc", "cgg", "cgt", "cta", "ctc", "ctg", "ctt", 
                  "gaa", "gac", "gag", "gat", "gca", "gcc", "gcg", "gct", 
                  "gga", "ggc", "ggg", "ggt", "gta", "gtc", "gtg", "gtt", 
                  "tac", "tat", "tca", "tcc", "tcg", "tct", "tgc", "tgg", 
                  "tgt", "tta", "ttc", "ttg", "ttt") 
# order is defined as per phangorn -> phyDat.R -> line 324
nucleotideAlphabet <- c("a", "c", "g", "t")

# estimate codon frequencies based on the F3X4 method
getF3X4codon_frequencies <- function(codon_data, codon_alphabet = CodonAlphabet,
                                     nucleotide_alphabet = nucleotideAlphabet) {
    epirical_codon_frequenices <- baseFreq(codon_data)
    # get the position-wise nucleotide frequencies
    nuc_freq_by_pos <- matrix(0, nrow=length(nucleotide_alphabet), ncol=3)
    for (i in seq_along(codon_alphabet)) {
        codon <- codon_alphabet[i]
        codon_freq <- epirical_codon_frequenices[i]
        for (j in 1:3) {
            nuc <- substr(codon,j,j)
            nuc_index <- which(nucleotide_alphabet == nuc)
            nuc_freq_by_pos[nuc_index,j] = nuc_freq_by_pos[nuc_index,j] + codon_freq
        }
    }
    # construct the codon frequencies according to the F3X4 method
    codon_frequencies <- rep(0,length(codon_alphabet))
    for (i in seq_along(codon_alphabet)) {
        codon <- codon_alphabet[i]
        codon_components <- unlist(strsplit(codon, ""))
        nuc1_freq <- nuc_freq_by_pos[which(nucleotide_alphabet==codon_components[1]),1] 
        nuc2_freq <- nuc_freq_by_pos[which(nucleotide_alphabet==codon_components[2]),2]
        nuc3_freq <- nuc_freq_by_pos[which(nucleotide_alphabet==codon_components[3]),3]
        codon_frequencies[i] <- nuc1_freq*nuc2_freq*nuc3_freq
    }
    return(codon_frequencies)
}


# return M
bf_by_codon <- function(x){
    y <- codon2dna(x)
    ny <- sum(attr(y, "weight"))
    M <- matrix(0, nrow=4, ncol=3)
    M[,1] <- subset(y, select=seq(1, ny, by=3), site.patter=FALSE) %>% baseFreq
    M[,2] <- subset(y, select=seq(2, ny, by=3), site.patter=FALSE) %>% baseFreq 
    M[,3] <- subset(y, select=seq(3, ny, by=3), site.patter=FALSE) %>% baseFreq
    M
}

# return
F3x4_freq <- function(M, CodonAlphabet, nucleotideAlphabet){
    pos <- CodonAlphabet %>% strsplit("") %>% unlist %>% 
        match(nucleotideAlphabet) %>% matrix(ncol=3, byrow = TRUE) 
    M[pos[,1],1] * M[pos[,2],2] * M[pos[,3],3]
}

# test
#SCRIPT_DIR = # please fill in the script dir
#codon_data_path = SCRIPT_DIR + "/codon_aln.fas"
#codon_data = read.phyDat(codon_data_path, format = "fasta", type="CODON")
#codon_frequencies = getF3X4codon_frequencies(codon_data)
