# GROUP WORK 1

# Jackson Cramer, s2274544
# Maximillian McCourt, s2145762
# Natalia Montalvo Cabornero, s1969053

# List of tasks completed by each group member:
# Jackson Cramer - primary focus on exercises 3 and 4.
# Max McCourt - primary focus on exercises 7 and 8.
# Natalia Montalvo Cabornero - primary focus on exercises 5 and 6.
# We collaborative contributed on the development and optimization of each exercise. While we divided forces, 
# we ensured  a roughly equal distribution of work trough code reviews, debugging sessions and meetings for problem solving. 



# Exercise 3

# set working directory and read shakespeare text into R environment
# setwd("GW1") ## comment out of submitted
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")


# Exercise 4

# to remove stage directions (words located between square brackets), first locate all indices of opening brackets
stage_start = grep("[",a,fixed=TRUE)
stage_dir <- c()
# loop through each opening bracket
for (left_bracket in stage_start){
  # find closing bracket indices in next 100 words after opening bracket
  stage_end = grep("]", a[(left_bracket):(left_bracket+100)], fixed=TRUE)
  # find index of first closing bracket after opening bracket, but first ensure there is a closing bracket
  if (length(stage_end) > 0) {
    right_bracket <- left_bracket-1+stage_end[1]
    # append indices corresponding to stage directions
    stage_dir <- c(stage_dir, left_bracket:right_bracket)
  }}
# remove all stage directions
a <- a[-stage_dir]


# remove character names (fully upper case) and arabic numerals, but do not remove I, A, or O
a <- a[which(a == "I" | a == "A" | a == "O" | a != toupper(a))]

# remove _ and - from words by replacing them by empty strings
a <- gsub("[_-]", "", a)

# to seperate punctuation marks from words, create a split_punct function that takes in word vector w and vector of punctuation marks p
split_punct <- function(w, p) {
# for each punctuation p, place a space between the punctuation and connected word
  for (punct in p) {
    w <- gsub(punct, paste(" ", punct, sep = ""), w, fixed=TRUE)
  }
# put word vector into one string, with each item seperated by a space
  w <- paste(w, collapse = " ")
# split punctuation and words by splitting on empty space, returning a vector of single words and single punctuation marks
  out <- strsplit(w, split = " ")[[1]]
  return(out)
}

p_to_use <- c(",", ".", ";", "!", ":", "?")
# apply split_punct function to word vector a using punctuation marks in p_to_use
a <- split_punct(a, p_to_use)

# convert word vector a to lower case
a <- tolower(a)


# Exercise 5

# find all unique words
unique_words <- unique(a)
# index mapping each unique word to it's position in the text
word_idx <- match(a, unique_words)
# tabulate function helps us to find the frequency of each unique word in the text
word_freq <- tabulate(word_idx)
# rank the words by frequency (the highest frequency is rank 1)
freq_ranks <- rank(-word_freq, ties.method = "average")
# get indices of top 1000 words
top_idx <- which(freq_ranks <= 1000)
# extract the actual top 1000 words
b <- unique_words[top_idx]


# Exercise 6

# index of the top words in the text, returning NA if not found
b_idx <- match(a, b)

# choose mlag, i.e. how many preceding words we are considering (here current word + 3 preceding)
mlag <- 4 
n <- length(a)

# initialise M matrix with dimension (n-mlag) x (mlag+1)
M <- matrix(nrow=n-mlag, ncol=mlag+1)
# fill the matrix: current word index in column 1 and then three lags in subsequent columns
for (m in 1:(mlag+1)){
  # slide through the text for each lag position
  M[,m] <- b_idx[m:(n-(mlag+1)+m)]
}


# Exercise 7

# function that randomly generates the next word of a given key using the text
next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){
  # inputs:
  #   key - word sequence for which the next word is generated
  #   M is matrix M with top word indexed and lags
  #   M1 is the indexing vector of the entire text (index_vector)
  #   w - vector of mixture weights
  
  # initialise next_word and match_rows
  next_word <- NA
  match_rows <- c()
  
  # finds the rows of M where the text matches the key
  while (length(match_rows) == 0){
    # select the number of column entries in M we need to look at based on key length
    mc <- mlag - length(key) +1
    
    # finds all rows of M where text matches key, labelled match_rows
    ii <- colSums(!(t(M[,mc:mlag,drop=FALSE])==key))
    match_rows <- which(ii == 0 & is.finite(ii))
    
    # now delete all match rows when next word after key is NA
    match_rows <- match_rows[!is.na(M[match_rows,(mlag+1)])]
    
    # if no matching rows, initialise key to be one word shorter for next loop
    key <- key[2:length(key)]
  }
  
  # chooses random row from match rows using sample function (indexing to avoid function issue)
  random_row <- match_rows[sample(length(match_rows), 1)]
  # take next word from last entry of random_row
  next_word <- M[random_row, mlag+1]
  return(next_word)
}


# Exercise 8 

# function that iteratively produces sentence using next.word function
sim.shakespeare <- function(key_length, M, M1, p, b){
  # inputs:
  #   key_length - desired key length for word generation
  #   M is matrix M with top word indexed and lags
  #   M1 is the indexing vector of the entire text (index_vector)#
  #   p - punctuation marks
  #   b - list of ~1000 most frequent words in text
  
  # first, set b where punctuation is not included
  b_no_p <- b[!b %in% p]
  
  # use b_no_p to initialise random starter token for sentence
  random_starter_token <- sample(b_no_p, 1)
  key <- which(b == random_starter_token)
  # initialise sentence as the random word
  token_sentence <- key
  
  # find index of period for ending our while loop
  period_index <- (which(b=="."))
  # while the last word in the sentence is not a period:
  while (token_sentence[length(token_sentence)] != period_index){
    
    # Take next word using next.word function
    next_word <- next.word(key, M, M1)
    # append next word to sentence
    token_sentence <- append(token_sentence, next_word)
    
    # if the sentence is shorter than the desired key_length, just make sentence the key
    if (length(token_sentence) < key_length){
      key <- token_sentence
    }else{
      # else take last key length no. of words of current token sentence
      key <- token_sentence[(length(token_sentence)-key_length+1):length(token_sentence)]
    }
  }
  # collapse sentence
  sentence <- paste(b[token_sentence], collapse=" ")
  return(sentence)
}

# Exercise 9

# call on sim.shakespeare function to generate shakespearean sentence
key_length <- 2
sim.shakespeare(key_length, M, M1, p_to_use, b)
