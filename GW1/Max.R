setwd("GW1") ## comment out of submitted
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")


### PREPROCESSING

# drop stage directions
stage_start = grep("[[]", a)

stage_directions <- c()
# loop for each opening bracket
for (i in stage_start){
  # find closing brackets indexes in next 100 words after opening bracket
  stage_end = grep("[]]", a[(i):(i+100)])
  # find index of first closing bracket after opening bracket
  j <- i-1+stage_end[1]
  # append indexes corresponding to stage directions
  stage_directions <- append(stage_directions, c(range(i, j)))
}
# dropped stage directions
a <- a[-unique(stage_directions)]


# Drop names and numbers

uppers_nums <- c()
for (i in 1:length(a_1)){
  # temporarily remove grammar to check character length
  i_temp <- gsub("[.,!?;:_']", "", a_1[i])
  # If uppercase and more than one character
  if (a_1[i] == toupper(a_1[i]) && nchar(i_temp)>1){
    uppers_nums <- append(uppers_nums, i)
  }
  # If numeric
  if (grepl("[0-9]", a[i])){
    uppers_nums <- append(uppers_nums, i)
  }
}
# remove names and numerals
a <- a[-uppers_nums]

# get rid of _ and -
a <- gsub("[_]", "", a)
a <- gsub("[-]", "", a)


# Ex 4: Make punctuation marks into their own words
punctuation_marks <- c("!", "\\.", ":", ";", ",", "\\?")

# Function that splits the punctuation into its own entry in our word vector
split_punct <- function(word_vector, punctuation){
  
  # Separate all strings with punctuation by introducing a space between the word and the punctuation
  for (p in punctuation){
    p_space <- paste0(" ", paste(p))
    #print(p_space)
    word_vector <- gsub(p, p_space, word_vector)
    #print(word_vector)
  }
  #print(word_vector)
  # We then remove 
  word_vector <- unlist(strsplit(word_vector, split = " "))
  
  return(word_vector)
  
}


system.time(a <- split_punct(a, punctuation_marks))


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
top_words <- unique_words[top_indices]


## TEXT SETUP

# Ex 6

# index of top words positions in text
index_top_words <- match(a, top_words)

# choose mlag 
mlag <- 4
n <- length(a_4)

# initialise M matrix with dimension (n-mlag) x (mlag+1)
M <- matrix(nrow=n-mlag, ncol=mlag+1)
# Put top words index into 1st colm of matrix, then three lags in subsequent colms
for (m in 1:(mlag+1)){
  M[,m] <- index_top_words[m:(n-(mlag+1)+m)]
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
sim_shakespeare <- function(key_length, M, M1, top_words){
  
  # Choose first word
  # first, take top_words where punctuation is not included
  punctuation <- c(",", ".", ";", ":", "!", "?")
  tw_no_punct <- top_words[!top_words %in% punctuation]
  
  # initialise random starter token
  random_starter_token <- sample(tw_no_punct, 1)
  key <- which(top_words == random_starter_token)
  key_sentence <- key
  
  # find index of period for ending our while loop
  period_index <- (which(top_words=="."))
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
    sentence <- paste(top_words[key_sentence], collapse=" ")
    print(sentence)
  }
  return(sentence)
}


# Let's simulate Shakespeare
sim_shakespeare(2, M, M1, top_words)





# the O's (not mentioned as a special case)
# ACT I - should we remove this separately for more complete analysis?
# how do we deal with -? remove, make two words or ignore?
#

