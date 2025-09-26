# names and university user names of each member of the group
# brief description of what each team member contributed to the project

#setwd("GW1") ## comment out of submitted
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")


### PREPROCESSING

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

a = a[which(a == "I" | a == "A" | a == "O" | a != toupper(a))]

# get rid of _ and -
a = gsub("[_-]", "", a)


# Ex 4: Make punctuation marks into their own words
split_punct = function(w, p) {
  for (punct in p) {
    if (punct == "." | punct == "?") {
      punct_escaped = paste("\\", punct, sep = "")
    } else {
      punct_escaped = punct
    }
    w = gsub(punct_escaped, paste(" ", punct, sep = ""), w)
  }
  w = paste(w, collapse = " ")
  out = strsplit(w, split = " ")[[1]]
  return(out)
}

p_to_use = c(",", ".", ";", "!", ":", "?")

a = split_punct(a, p_to_use)


# make all lower case for simplicity
a <- tolower(a)


# Ex 5: Matching

# all unique words
unique_words <- unique(a_4)
# index of each word's occurence in the text
index_vector <- match(a_4, unique_words)
# word frequencies of unique words
word_frequencies <- tabulate(index_vector)
# rank the frequencies
frequency_ranks <- rank(-word_frequencies, ties.method = "average") # check average treats same no. occurences equally
# Get indices of top 1000 words
top_indices <- which(frequency_ranks <= 1000)
# Extract the top 1000 words
b <- unique_words[top_indices]


## TEXT SETUP

# Ex 6

# index of top words positions in text
index_b <- match(a, b)

# choose mlag 
mlag <- 4
n <- length(a_4)

# initialise M matrix with dimension (n-mlag) x (mlag+1)
M <- matrix(nrow=n-mlag, ncol=mlag+1)
# Put top words index into 1st colm of matrix, then three lags in subsequent colms
for (m in 1:(mlag+1)){
  M[,m] <- index_b[m:(n-(mlag+1)+m)]
}


## SENTENCE GENERATION

# Ex 7

# Function that randomly generates the next word of a given key using the text
next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){
  # key - word sequence for which the next word is generated
  # M is matrix M with top word indexed and lags
  # M1 is the indexing vector of the entire text (index_vector)
  # w - vector of mixture weights
  key_length <- length(key)
  mc <- mlag - key_length +1
  # initialise next_word as NA
  next_word <- NA
  
  # finds the rows of M where the text matches the key
  ii <- colSums(!(t(M[,mc:mlag,drop=FALSE])==key))
  match_rows <- which(ii == 0 & is.finite(ii))
  print(match_rows)

  # now randomly selects the next word following the key in the text
  # (if NA, chooses again)
  while (is.na(next_word)){
  random_row <- sample(match_rows, 1)
  next_word <- M[random_row,mlag+1]
  }
  return(next_word)
}

# Ex 8 + 9

# Function that iteratively produces sentence using next.word function
sim_shakespeare <- function(key_length, M, M1, b){
  
  # Choose first word
  # first, take b where punctuation is not included
  punctuation <- c(",", ".", ";", ":", "!", "?")
  tw_no_punct <- b[!b %in% punctuation]
  
  # initialise random starter token
  random_starter_token <- sample(tw_no_punct, 1)
  key <- which(b == random_starter_token)
  key_sentence <- key
  
  # find index of period for ending our while loop
  period_index <- (which(b=="."))
  while (key_sentence[length(key_sentence)] != period_index){
    
    print(key)
    
    
    # Take next word using next.word function
    next_word <- next.word(key, M, M1)
    # append next_word to sentence
    key_sentence <- append(key_sentence, next_word)
    
    # if the sentence is shorter than the desired key_length, just make sentence the key
    if (length(key_sentence) < key_length){
      key <- key_sentence
    }else{
      # take last key length no. of words of sentence
      key <- key_sentence[(length(key_sentence)-key_length+1):length(key_sentence)]
    }
    # collapse sentence
    sentence <- paste(b[key_sentence], collapse=" ")
    print(sentence)
  }
  return(sentence)
}


# Let's simulate Shakespeare
sim_shakespeare(2, M, M1, b)





# the O's (not mentioned as a special case)
# ACT I - should we remove this separately for more complete analysis?
# how do we deal with -? remove, make two words or ignore?
#
