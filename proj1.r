# Jackson Cramer (s2274544), Maximillian McCourt (s2145762), Natalia Montalvo Cabornero (s1969053)

# As a three, we collaborated on the development and optimisation of each exercise. While we divided forces, 
# we ensured  a roughly equal distribution of work through code reviews, debugging sessions and meetings for problem solving. 
# Jackson focused on question 4, and code optimisation/efficiency. Max focused on the functions in exercises 7,8 and 9, and Natalia
# focused on exercises 5 and 6, and concise, clear code commenting.

# This code builds a model that uses the complete works of Shakespeare to generate sentences
# that mimic the structure of Shakespearean writing.

###### fix comments, remove all "exercises" --> explain conceptual meaning of code and purpose


# Exercise 3

# set working directory and read shakespeare text into R environment
# setwd("GW1")
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")


# Exercise 4

# to remove stage directions (words located between square brackets), first locate all indices of opening brackets
stage_start = grep("[",a,fixed=TRUE)
# initialise empty vector to store indices of stage directions
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

# remove _ and - from words by replacing them with empty strings
a <- gsub("[_-]", "", a)

# function that separates punctuation marks from words they are attached to
split_punct <- function(w, p) {
  # inputs:
  #   w - word vector
  #   p - punctuation marks vector
  
  # for each punctuation p, place a space between the punctuation and connected word
  for (punct in p) {
    w <- gsub(punct, paste(" ", punct, sep = ""), w, fixed=TRUE)
  }
  # put word vector into one string, with each item separated by a space
  w <- paste(w, collapse = " ")
  # split punctuation and words by splitting on empty space, returning a vector of single words and single punctuation marks
  out <- strsplit(w, split = " ")[[1]]
  return(out)
}

# apply split_punct function to word vector a, using punctuation marks in p_to_use
p_to_use <- c(",", ".", ";", "!", ":", "?")
a <- split_punct(a, p_to_use)

# convert word vector a to lower case
a <- tolower(a)


# Exercise 5

# find all unique words
unique_words <- unique(a)
# index mapping each unique word to its position in the text
word_idx <- match(a, unique_words)
# tabulate function helps us to find the frequency of each unique word in the text
word_freq <- tabulate(word_idx)
# rank the words by frequency, the highest frequency being rank 1 (assigns words with same frequency the same average rank)
freq_ranks <- rank(-word_freq, ties.method = "average")
# get indices of highest-ranking ~1000 words
top_idx <- which(freq_ranks <= 1000)
# extract ~1000 most-frequent unique words
b <- unique_words[top_idx]


# Exercise 6

# index of the most frequent words in the text (assigns NA if not in b)
frequent_idx <- match(a, b)

# choose mlag, i.e. how many preceding words we are considering (here current word + 3 preceding)
mlag <- 4 
n <- length(a)
# To construct M matrix of common word token sequences, initialise empty matrix with dimension (n-mlag) x (mlag+1)
M <- matrix(nrow=n-mlag, ncol=mlag+1)
# fill the matrix with current word index in column 1, and then four lags in subsequent columns
for (m in 1:(mlag+1)){
  # shift through the text for each lag position
  M[,m] <- frequent_idx[m:(n-(mlag+1)+m)]
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
  match_words <- list()
  word_weights <- c()
  
  # for each length of the key given
  for (i in 1:length(key)){
    # i-th key length
    key_i <- key[i:length(key)]
    # select the number of column entries in M we need to look at based on key length
    mc <- mlag - length(key_i) +1
    # finds all rows of M where text matches key, labelled match_rows
    ii <- colSums(!(t(M[,mc:mlag,drop=FALSE])==key_i))
    match_rows <- which(ii == 0 & is.finite(ii))
    # now delete all match rows when next word after key is NA
    match_rows <- match_rows[!is.na(M[match_rows,(mlag+1)])]
    # all words that follow matching rows
    match_words[[i]] <- M[match_rows, mlag+1]
    # find weights of each word for probabilities 
    # weight indexing adjusted to ensure correct weight is assigned dependent on key length (first weight assigned to longest key length, etc)
    word_weights <- append(word_weights, w[length(w)-length(key_i)+1]*rep(1, length(match_words[[i]])))
  }
  # randomly select one of the matching words, with probabilities weighted by word_weights
  next_word <- sample(unlist(match_words), 1, prob = word_weights)
  return(next_word)
}


# Exercise 8 

# function that iteratively produces sentence using next.word function
sim.shakespeare <- function(max_key_len, M, M1, p, b, w=rep(1,ncol(M)-1)){
  # inputs:
  #   key_length - desired key length for word generation
  #   M is matrix M with top word indexed and lags
  #   M1 is the indexing vector of the entire text (index_vector)#
  #   p - punctuation marks
  #   b - list of ~1000 most frequent words in text
  #   w - weights for key lengths
  
  # set b such that punctuation is not included
  b_no_p <- b[!b %in% p]
  # use b_no_p to initialise random starter token for sentence
  random_starter_token <- sample(b_no_p, 1)
  key <- which(b == random_starter_token)
  # initialise token sentence as the random word token
  token_sentence <- key
  
  # find index of full stop for ending our while loop
  period_index <- (which(b=="."))
  # while the last token in the token sentence is not the token of full stop:
  while (token_sentence[length(token_sentence)] != period_index){
    # take next word using next.word function
    next_word <- next.word(key, M, M1, w)
    # append next word token to token sentence
    token_sentence <- append(token_sentence, next_word)
    
    # if token sentence is shorter than max key length, just make updated key the current token sentence
    if (length(token_sentence) < max_key_len){
      key <- token_sentence
    }else{
      # else assign key as last max key length no. of words of current token sentence
      key <- token_sentence[(length(token_sentence)-max_key_len+1):length(token_sentence)]
    }
  }
  # convert token sentence into words, and collapse using space
  sentence <- paste(b[token_sentence], collapse=" ")
  return(sentence)
}


# Exercise 9

# call on sim.shakespeare function to generate shakespearean sentence
# choose max key length to be 3, let weights be default (equally weighted)
max_key_length <- 3
sim.shakespeare(max_key_length, M, M1, p_to_use, b)
