setwd("GW1") ## comment out of submitted
a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")
a_full <- a

#a <- a[1:50000]
a

# PREPROCESSING

## drop stage directions
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
a_1 <- a[-unique(stage_directions)]

# make each punctuation its own word




# drop names and numerals

uppers_nums <- c()
for (i in 1:length(a_1)){
  # temporarily remove grammar to check character length
  i_temp <- gsub("[.,!?;:_']", "", a_1[i])
  # If uppercase and more than one character
  if (a_1[i] == toupper(a_1[i]) && nchar(i_temp)>1){
    uppers_nums <- append(uppers_nums, i)
  }
  # If numeric
  if (grepl("[0-9]", a_1[i])){
    uppers_nums <- append(uppers_nums, i)
  }
}
# remove names and numerals
a_2 <- a_1[-uppers_nums]

# Also maybe remove all I's for when it says ACT I and SCENE I?
#ACT_words <- which(a_1 == "ACT")
#a_1[ACT_words+1]

# get rid of _ (chose not to get rid of - for now)
a_2 <- gsub("[_]", "", a_2)

# Ex 4: Make punctuation marks into their own words
punctuation_marks <- c("!", ".", ":", ";", ",", "?")

split_punct <- function(word_vector, punctuation) {
  # input: word vector, punctuation mark list
  # separates punctuation from words and returns a new vector where
  # each punctuation mark is its own separate entry
  
  # we need to get punctuation list into correct formatting
  punct <- paste0("[", paste(ifelse(punctuation == ".", "\\.", punctuation), collapse=""), "]")
  
  # identify which words in the vector contain punctuation
  punct_indexes <- grep(punct, word_vector)
  
  # Initiate empty vector with maximum possible size (2x length in the case every word has punctuation)
  max_length <- 2 * length(word_vector)
  vec <- character(max_length)
  
  # initialise an index to track our position in the result vector
  pos <- 1
  
  # looping through each word...
  for (i in 1:length(word_vector)) {
    
    # Check if this word has punctuation at the end
    if (i %in% punct_indexes) {
      
      # Step 5a: Extract the word part (without punctuation)
      # Using regex backreference \\1 to get the first capture group
      clean_k <- gsub(punct, "\\1", word_vector[i])
      
      # find the punctuation mark
      k_string <- strsplit(word_vector[i], "")[[1]]
      is_punct <- which(k_string %in% punctuation)
      punct_k <- k_string[is_punct]
    
      # add the clean word and punctuation mark to the result vector
      vec[pos] <- clean_k
      vec[pos+1] <- punct_k
      pos <- pos + 2
      
    } else {
      # if no punctuation, add as is
      vec[pos] <- word_vector[i]
      pos <- pos + 1
    }
  }
  
  # final result is vector minus empty spaces
  result <- vec[1:(pos - 1)]
  return(result)
}

split_punct(c("a?", "b,", "c.", "d;", "e:", "f!"), punctuation_marks)
system.time(a_3 <- split_punct(a_2, punctuation_marks))

# make all lower case for simplicity
a_4 <- tolower(a_3)

a_4

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
top_words


# Ex 6

# index of top words positions in text
index_top_words <- match(a_4, top_words)

# choose mlag 
mlag <- 4
n <- length(a_4)

# initialise M matrix with dimension (n-mlag) x (mlag+1)
M <- matrix(nrow=n-mlag, ncol=mlag+1)
# Put top words index into 1st colm of matrix, then three lags in subsequent colms
for (m in 1:(mlag+1)){
  M[,m] <- index_top_words[m:(n-5+m)]
}

M[,1]

# Ex 7


next.word <- function(key, M, M1, w=rep(1,ncol(M)-1)){
  # key - word sequence for which the next word is generated
  # M is matrix M with top word indexed and lags
  # M1 is the indexing vector of the entire text (index_vector)
  # w - vector of mixture weights
  key_length <- length(key)
  mc <- mlag - length(key) + 1
  
  ii <- colSums(!(t(M[,mc:mlag,drop=FALSE])==key))
  print(ii)
  
  # matching rows, where key matches rows of M
  matching_rows <- integer(0)
  for (r in 1:nrow(M)) {
    if (all(M[r, mc:mlag] == key, na.rm=TRUE)) {  # Compare the appropriate columns
      matching_rows <- c(matching_rows, r)
    }
  }
  print(matching_rows)
  
  #for (r in nrow(M)){
  #  if (key == M[r,(1:length(key))]){
  #    
  #  }
}




ncol(M)

key = c("bare")

next.word(key,M,M1,w=rep(1,ncol(M)-1))

