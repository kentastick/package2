
#
prime_n <- function(){
  x <- as.numeric(readline("Natural number? "))
  is_integer <- function(values){
    all(values %% 1 == 0)
  }
  for(i in (2:x)){
    if(is_integer(x / i) == T) break
  }
  if(i == x) cat(paste0(x, " is a prime number!"))
  else cat(paste0(x, " is NOT a prime number!"))
}

is.prime <- function(num){
  # if num is integer
  if (num != as.integer(num)) return("Input is not integer!")
  # if num is less than 2
  if (num < 2){
    return("Input is wrong!")
  }
  # if num is 2 or 3
  else if (num == 2 || num == 3) {
    return(TRUE)
  } else {

    for (i in 2:ceiling(num/2)){

      if ((num %% i) == 0) {

        return(FALSE)

      }

    }

    return(TRUE)
  }
}




#stocklist
a <- c("ABART 20", "CDXEF 50", "BKWRK 25", "BTSQZ 89", "DRTYM 60")

sapply(c("A", "B"), function(x)a[grepl(pattern = paste0("^",x),a)])



#bouce ball
bouncingBall <- function(h, bounce, window) {
  if(h<0|bounce==1|window>h)return(-1)

  count =0
  while(TRUE){
    h = h*bounce
    if(h<window){
      count = count + 1; break
    }else{
      count = count + 2
    }
  }
  print(count)
}

bouce(3, 0.66, 1.5)


#fizz puzz
seq(15,20, by = 15)

solution <- function(n) {
  c((n-1)%/%3 -((n-1)%/%5)%/%3, (n-1)%/%5-(n-1)%/%15, (n-1)%/%15)
}

solution <- function(n){
  s<-1
  i<-0
  j<-0
  k<-0
  while (min(c(s*3, s*5))<n){
    if( 3*s<n){
      i<-i+1
    }
    if (5*s<n){
      j<-j+1
    }
    if (15*s<n){
      i<-i-1; j<-j-1; k<-k+1;
    }
    s<-s+1
  }
  return(c(i,j,k))
}

solution <- function(n) {
  c(length(setdiff(seq(3,n-1, by = 3),seq(5,n-1, by = 5))), length(setdiff(seq(5,n-1, by = 5),seq(3,n-1, by = 3))), (n-1)%/%15)
}

solution(20)

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
