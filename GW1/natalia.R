# names + student numbers
# brief description of work


# setwd("/Users/natmo/OneDrive/Escritorio/MSc/SEM1/Extended_Stat/GW1")

a <- scan("shakespeare.txt",what="character",skip=83,nlines=196043-83,
          fileEncoding="UTF-8")    

###### 4. CLEANING THE TEXT

###### 4.a. Remove the stage directions 

# Find the index of the start of the stage direction 
index_stage_dir <- grep("[",a,fixed=TRUE)

stage_dir_delete <- integer(0)   # we'll store the indexes

# Find the end of that stage direction
for (i in index_stage_dir) {
  # we need to check if it's only one word, "[]" together
  
  range = a[(i+1):(i+100)]    # the hundred words after the "["
  end_index <- grep("]",range,fixed=TRUE)  #index in range of words in that range with "]"

  if (length(end_index) > 0) {    #if we find an end
  
    first_end <- i + end_index[1] # index in a of the first end
    
    stage_dir_delete <- c(stage_dir_delete, i:first_end)   
  }
}

a <- a[-stage_dir_delete]

###### 4.b Remove words fully Upper case and numbers

# find the upper words 

upper_words <- a == toupper(a) # return logic, if si son mayusculas en new_a
words_to_keep <- a %in% c("I", "A") #returns a logic, if a word in new_A is == to a word in c("I", "A") --> true
remove_upper <- which(upper_words & !words_to_keep) # upper words and false
a <- a[-remove_upper]

numbers <- grep("^[0-9]+$", a)
a <- a[-numbers]


###### 4.c Remove “-” and “_” from words

a <- gsub("-_", "", a)


###### 4.d,e Split punctuation Function


x <- c("hello,", "world!", "this", "is.", "great")

split_punct <- function(){
  parts <- sub("^(\\w+)([[:punct:]]+)$", "\\1 \\2", xi)
  strsplit(parts, " ", fixed = TRUE)[[1]]
}
<<<<<<< HEAD

## need to check this
=======
### need to check this
>>>>>>> 5a2631c3e301e0ae5e5cea31122d84217ff94570

###### 4.f Vector in lower case

a <- tolower(a)


###### 5. find most frequent words

# a 
# is the vector of the unique words of our new clean set of a

b <- unique(a)

# b - Use match to find the vector of indices indicating which element in the unique word vector each element in the text corresponds to
unique_idx <-  match(a, b)     

# c - how many time each unique word occurs in the text.

freq <- tabulate(unique_idx)

# d - vector b containing the ≈ 1000 most common words

n = 1000

freq_rank <- rank(-freq, ties.method = "average")    # rank of the frequencies, with ties.method = "first" the last 1000 should be the most frequent


top_idx <- which(freq_rank <= n)      # b <= n creates a logical logical vector that’s the same length as b, marking which words are in the top n_common
top_words <- b[top_idx]


## 6. matrices of common word token sentences.

# index of top words in text

top_words_idx <- match(a, top_words)

# how forward we look
mlag = 4    # we can change and see how it works
n = length(a)

M <- matrix(n-mlag, mlag+1)

# add words in first column and 
for (m in 1:(mlag+1)){
  M[,m] <- index_top_words[m:(n-5+m)]
}

## 7. matrices of common word token sentences.

# next.word <- function(key,M,M1,w=rep(1,ncol(M)-1))



  
  
  
## 8. Start key

p <- c(",", ".", ";", "!", ":", "?")

# select a random word 

chosen_word <- "romeo"
chosen_word <- sample(a, 1)
if (chosen_word %in% words| !(chosen_word %in% p)){
  start_key <- chosen_word
}











