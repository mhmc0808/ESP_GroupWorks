# GROUP WORK 1

# Jackson Cramer, s2274544
# Maximillian McCourt, s2145762
# Natalia Montalvo Cabornero, s1969053

# brief description of what each team member contributed to the project




# Exercise 3

# setwd("GW1") ## comment out of submitted
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")


# Exercise 4

# drop stage directions
stage_start = grep("[",a,fixed=TRUE)

stage_directions <- c()

# loop for each opening bracket
for (left_bracket in stage_start){
  # find closing brackets indexes in next 100 words after opening bracket
  stage_end = grep("]", a[(left_bracket):(left_bracket+100)], fixed=TRUE)
  # find index of first closing bracket after opening bracket
  if (length(stage_end) > 0) {
  right_bracket <- left_bracket-1+stage_end[1]
  # append indexes corresponding to stage directions
  stage_directions <- c(stage_directions, left_bracket:right_bracket)
}}
# dropped stage directions
a <- a[-stage_directions]


# Drop names and numbers

a <- a[which(a == "I" | a == "A" | a == "O" | a != toupper(a))]

# get rid of _ and -
a <- gsub("[_-]", "", a)


# Make punctuation marks into their own words
split_punct <- function(w, p) {
  for (punct in p) {
    w <- gsub(punct, paste(" ", punct, sep = ""), w, fixed=TRUE)
  }
  w <- paste(w, collapse = " ")
  out <- strsplit(w, split = " ")[[1]]
  return(out)
}

p_to_use <- c(",", ".", ";", "!", ":", "?")

a <- split_punct(a, p_to_use)


# make all lower case for simplicity
a <- tolower(a)


# Exercise 5

# all unique words
unique_words <- unique(a)
# index of each word's occurence in the text
index_vector <- match(a, unique_words)
# word frequencies of unique words
word_frequencies <- tabulate(index_vector)
# rank the frequencies
frequency_ranks <- rank(-word_frequencies, ties.method = "average") # check average treats same no. occurences equally
# Get indices of top 1000 words
top_indices <- which(frequency_ranks <= 1000)
# Extract the top 1000 words
b <- unique_words[top_indices]


# Exercise 6

# index of top words positions in text
index_b <- match(a, b)

# choose mlag 
mlag <- 4
n <- length(a)

# initialise M matrix with dimension (n-mlag) x (mlag+1)
M <- matrix(nrow=n-mlag, ncol=mlag+1)
# Put top words index into 1st colm of matrix, then three lags in subsequent colms
for (m in 1:(mlag+1)){
  M[,m] <- index_b[m:(n-(mlag+1)+m)]
}


# Exercise 7

# Function that randomly generates the next word of a given key using the text
next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){
  # key - word sequence for which the next word is generated
  # M is matrix M with top word indexed and lags
  # M1 is the indexing vector of the entire text (index_vector)
  # w - vector of mixture weights
  # initialise next_word and match_rows as NA
  next_word <- NA
  match_rows <- c()
  
  # finds the rows of M where the text matches the key
  while (length(match_rows) == 0){
    # select the number of column entries in M we need to look at based on key length
    mc <- mlag - length(key) +1
    
    # ii outputs sums of true(=0) and false(=1) based on if word matches key
    ii <- colSums(!(t(M[,mc:mlag,drop=FALSE])==key))
    # if a row has sum 0 (ie all true) and finite (no NA), the row matches the key
    match_rows <- which(ii == 0 & is.finite(ii))
    # now delete all match rows when last entry is NA
    match_rows <- match_rows[!is.na(M[match_rows,(mlag+1)])]
    # if there are no matching rows, we initialise key to be one word shorter
    # than previously until we have matching rows
    key <- key[2:length(key)]
  }

  # chooses random row from match rows (by indexing to avoid sample function issue)
  random_row <- match_rows[sample(length(match_rows), 1)]
  # take next word from last entry of random_row
  next_word <- M[random_row, mlag+1]
  return(next_word)
}


# Exercise 8 

# Function that iteratively produces sentence using next.word function
sim_shakespeare <- function(key_length, M, M1, p, b){
  
  # Choose first word
  # first, take b where punctuation is not included
  b_no_p <- b[!b %in% p]
  
  # initialise random starter token
  random_starter_token <- sample(b_no_p, 1)
  key <- which(b == random_starter_token)
  token_sentence <- key
  
  # find index of period for ending our while loop
  period_index <- (which(b=="."))
  while (token_sentence[length(token_sentence)] != period_index){
    
    # Take next word using next.word function
    next_word <- next.word(key, M, M1)
    # append next_word to sentence
    token_sentence <- append(token_sentence, next_word)
    
    # if the sentence is shorter than the desired key_length, just make sentence the key
    if (length(token_sentence) < key_length){
      key <- token_sentence
    }else{
      # take last key length no. of words of current token sentence
      key <- token_sentence[(length(token_sentence)-key_length+1):length(token_sentence)]
    }
  }
  # collapse sentence
  sentence <- paste(b[token_sentence], collapse=" ")
  return(sentence)
}

# Exercise 9

# Let's simulate Shakespeare
sim_shakespeare(2, M, M1, p_to_use, b)




