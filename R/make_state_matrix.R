#' Create a State Matrix
#'
#' This function generates state matrices based on the given dataset, with options for
#' defining ID columns, LTC columns, and matrix dimensions. It splits the data, processes
#' it by group, and computes matrices for each subset of data.
#'
#' @param data A data frame or tibble containing the dataset long format
#' @param id A string representing the column name of patients ID
#' @param ltc A string representing the column name of the event (long term conditions)
#' @param aos A string representing the column name of age at onset of considered long term conditions
#' @param k An integer representing the group size when splitting the dataset (usefull for large datasets)
#' @param l An integer representing the number of rows of the state matrices. It usualle represent the maximum age of patients
#' @param fail_code A string representing the code used for failure (death)
#' @param cens_code A string representing the code used for censoring
#'
#' @return A list of state matrices, one per patient
#'
#' @details The function processes the data by splitting it based on the `id` and `ltc` columns.
#' It computes matrices for each subset, with rows representing different states and columns representing unique IDs.
#'
#' @importFrom dplyr pull filter mutate %>%
#' @importFrom Matrix Matrix
#' @importFrom tibble add_column
#' @importFrom rlang sym
#' @importFrom utils tail
#' @export


# library(dplyr)
# library(Matrix)
# library(tibble)
# library(rlang)

# data = data6 ; id = 'link_id' ; ltc = 'reg' ; aos = 'aos' ; k = 100 ; l = 111 ; fail_code = 'death'
# cens_code = 'cens'


make_state_matrix <- function( data , id = 'link_id' , ltc = 'reg' , aos = 'aos' ,
                               k = 100 , l = 111,
                               fail_code = 'death', cens_code = 'cens'){

  ltc3 <- setdiff( unique( data %>% pull( !!ltc ) ) , c( NA ) )
  lid3 <- unique( data %>% pull( id ) )
  n <- length( ltc3 )
  s <- length( lid3 )

  #######################
  s1.1 <- as.numeric( as.factor( (data %>% pull( !!id )  ) ) )
  index <-  ( s1.1 ) %/% k

  df1 <- data.frame( s1.1 , index )

  data <- data %>%
    add_column( index)

  sd6 <- split( data , index )
  k2 <- length( sd6 )

  cs2 <- cumsum(sapply( sd6 , function( x ) length( unique( x %>% pull( !!id )) )) )
  df2 <- data.frame( cs1 = c( 0 , cs2[-length(cs2)] )+1 , cs2 )

  res <- as.list( seq_along(1:k2)  )
  j<-0
  #ii <- 1
  for( ii in 1:k2  ){

    data2 <- data %>%
      filter( get(id) %in% lid3[df2[ii,1]:df2[ii,2]])
    lid3.2 <- unique( data2 %>% pull( !!id ) )
    Mat <- Matrix( 0 , ncol = length( lid3.2 ) ,
                   nrow = n * l ,
                   dimnames = list(  as.character(rep( 0:( l - 1 ) , n ) ) , as.character( lid3.2 ) ) )

    # Initialize matrices
    mat2 <- Matrix(  0:( l - 1 )   , ncol =  n  , nrow = l  ,  dimnames = list(   as.character(0:(l-1)) , ltc3 ) )
    mat0 <- Matrix(  0   , ncol =  n  , nrow = l  ,  dimnames = list(  as.character(0:(l-1)) , ltc3 ) )
    ####
    i<- lid3.2[1]
    #i <- "K_47274"
    for( i in lid3.2 ) {

      j <- j + 1
      #cat('---\n')
      #print(j)
      #print(i)
      data3 <- data2 %>%
        filter( !!sym(id)  == i)
      ltc4 <- data3 %>% pull( !!ltc )
      m1 <- mat2[, ltc4 ] >=  matrix( data3 %>% pull(  !!aos  ), ncol = nrow( data3 ) ,
                                      nrow = l , byrow = T )
      # Get index of censoring column
      i2 <- colnames(m1) %in%  c(fail_code , cens_code)
      ## Get multiplying vector to censor m1 with
      # delay censoring of one year if it happens the same year of the onset of one or more other ltc
      if( sum(Matrix::colSums(m1[,!i2,drop = F] == m1[,i2]) == l) )
        m1[,i2] <-  c(FALSE , m1[-l,i2])

      cens_2 <- ifelse( m1[,i2] , NA , 1)


      # Censoring m1
      m1 <- m1 * cens_2
      # Get the whole matrix
      mat0[, ltc4 ] <- m1

      Mat[,i] <- as.vector(mat0)
      #plot( Mat[,i] <- as.vector(mat0) )

      # re initialize the matrix
      mat0 <- Matrix(  0   , ncol =  n  , nrow = l  ,  dimnames = list(  as.character(0:(l-1)) , ltc3 ) )

      if( i == tail(lid3.2,1)){
        #save(Mat , file = paste0('./RES_INT/MAT/Mat_',(j-(k-1)),'_',j,'.Rdata' ))
        res[[ii]] <- Mat
      }

    }
  }
  mat <- as.matrix(do.call(cbind , res))
  return( mat )

}
