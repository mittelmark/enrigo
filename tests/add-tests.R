library(enrigo)

x = 2 + 2

if (x != 4) {
    stop("Error: Something strange happened!!")
}

add <- function (x,y) { return(x+y) }
x = add(2,2)

if (x != 4) {
    stop("Error: Something strange happened in add(2,2)!!")
}
