nepal2Df = read.csv("nepal_gps.csv")

data_neg_pos <- nepal2Df



data_neg_pos$ermB2 <- replace(nepal2Df$ermB, nepal2Df$ermB == "NEG", "ermbNEG") 
data_neg_pos$ermB2 <- replace(nepal2Df$ermB, nepal2Df$ermB == "POS", "ermbPOS")
data_neg_pos$mefA2 <- replace(nepal2Df$mefA, nepal2Df$mefA == "NEG", "mefANEG") 
data_neg_pos$mefA2 <- replace(nepal2Df$mefA, nepal2Df$mefA == "POS", "mefAPOS")
data_neg_pos$folA_I100L2 <- replace(nepal2Df$folA_I100L, nepal2Df$folA_I100L == "NEG", "folA_I100LNEG")
data_neg_pos$folA_I100L2 <- replace(nepal2Df$folA_I100L, nepal2Df$folA_I100L == "POS", "folA_I100LPOS")
data_neg_pos$cat2 <- replace(nepal2Df$cat, nepal2Df$cat == "NEG", "catNEG")
data_neg_pos$cat2 <- replace(nepal2Df$cat, nepal2Df$cat == "POS", "catPOS")
data_neg_pos$PCV7B <- replace(nepal2Df$PCV7, nepal2Df$PCV7 == "Y", "PCV7Y")
data_neg_pos$PCV7B <- replace(nepal2Df$PCV7, nepal2Df$PCV7 == "N", "PCV7N")
data_neg_pos$PCV10B <- replace(nepal2Df$PCV10, nepal2Df$PCV10 == "Y", "PCV10Y")
data_neg_pos$PCV10B <- replace(nepal2Df$PCV10, nepal2Df$PCV10 == "N", "PCV10N")
data_neg_pos$PCV13B <- replace(nepal2Df$PCV13, nepal2Df$PCV13 == "Y", "PCV13Y")
data_neg_pos$PCV13B <- replace(nepal2Df$PCV13, nepal2Df$PCV13 == "N", "PCV13N")
data_neg_pos$PCV20B <- replace(nepal2Df$PCV20, nepal2Df$PCV20 == "Y", "PCV20Y")
data_neg_pos$PCV20B <- replace(nepal2Df$PCV20, nepal2Df$PCV20 == "N", "PCV20N")

data_neg_pos
