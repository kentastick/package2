



isValidWalk()
isValidWalk <- function(walk){
  all(sapply(X = c("n","w"), FUN = function(x)sum(grepl(x,walk))) ==sapply(X = c("s","e"), FUN = function(x)sum(grepl(x,walk))))

}


digit_cal <- function(number) {
  original_number = number
  k= 0
  value = number
  while(value >10){
    k =k+1
    value <- number%/%10**k
  }
  digit_number <- vector(length = k+1)
  i = 1
  while(k>=0){
    digit_number[i] <- number%/%10**k
    number <- number - digit_number[i]*10**k
    i = i+1
    k = k-1
  }
  if(original_number == sum(factorial(digit_number))){
    print("STRONG!!!!")
  }else{
    print("Not Strong !!")

  }
}


value = 0
count = 1
for(k in 1:10){
  for(i in 1:10){
    current_value =  1/i*(k+2)**2*k
    value = value + 1/i*(k+2)**2*k
    print(current_value)
    print(count)
    count = count+1
  }


}


#whether prime number or not
is_prime <- function(n){
  if(n <=1){
    return(FALSE)
  }else if(n ==2){
    return(TRUE)
  }else if (n > 2){
    if(any(n %%2:(n-1) ==0)){
        return(FALSE)
      }
  }
  return(TRUE)
  }

is_prime <- function(num) {
  flag = 0
  # prime numbers are greater than 1
  if(num > 1) {
    # check for factors
    flag = 1
    for(i in 2:(num-1)) {
      if ((num %% i) == 0) {
        flag = 0
        break
      }
    }
  }
  if(num == 2)    flag = 1
  if(flag == 1) {
    return(TRUE)
  } else {
    return(FALSE)
  }
}

find_missing(c(1,3,5,9,11))




find_missing <- function(seq_number) {
  l <- length(seq_number)
  max_n <- seq_number[l]
  min_n <- seq_number[1]
  ditch  <- (max_n-min_n)/l
  complete_seq <- seq(min_n,max_n,ditch)
  setdiff(complete_seq,seq_number)
}
