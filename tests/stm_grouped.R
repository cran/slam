

library("slam")
set.seed(201008)

## test
x <- matrix(sample(c(0,1), 100L, TRUE, prob = c(.9,.1)), 5L)
dim(x)
INDEX <- sample(1:4, 20L, TRUE)

s <- as.simple_triplet_matrix(x)
z <- as.matrix(s)

identical(row_tsums(z, INDEX), as.matrix(row_tsums(s, INDEX)))
identical(col_tsums(t(z), INDEX), as.matrix(col_tsums(t(s), INDEX)))

## NA indexes
k <- INDEX
is.na(k) <- k == 1L
any(is.na(k))
identical(as.matrix(row_tsums(s, k)), row_tsums(z, k))

## other data types
s$v <- as.integer(s$v)

identical(row_tsums(z, INDEX), as.matrix(row_tsums(s, INDEX)))

## NA values
is.na(s$v) <- 1:2
z   <- as.matrix(s)
z[] <- as.double(z) # coerce

identical(row_tsums(z, INDEX), as.matrix(row_tsums(s, INDEX)))
identical(row_tsums(z, INDEX, na.rm = TRUE), 
	  as.matrix(row_tsums(s, INDEX, na.rm = TRUE)))

##
s$v <- as.double(s$v)

identical(row_tsums(z, INDEX, na.rm = TRUE), 
	  as.matrix(row_tsums(s, INDEX, na.rm = TRUE)))

###
