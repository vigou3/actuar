## building the object
library(actuar)
data(hachemeister)
hach <- hachemeister$claims
age <- unlist(lapply(rbinom(5, 50, 0.5), seq, length = 12))
unit <- rep(rpois(5, 5), each = 12)
hach <- data.frame(claim = as.vector(hach),
                   year = rep(1:12, 5),
                   contract = rep(1:5, each = 12),
                   age = age,
                   unit = unit)
#### Bühlmann-Straub

## Xit = m + Yi + Yit

claim ~ contract + contract:year ## expanded form
claim ~ contract + year %in% contract ## preferred
f <- claim ~ contract/year ## equivalent (recommended to be parsinomious)

attributes(terms(lm(f, data = hach)))$factors ## Simply to see the representation of the formula

#### Crossed Classification

## Xijt = m + Yi + Yj + Yij + Yijt

claim ~ contract + age + contract:age + contract:age:year ## expanded form
claim ~ contract/age/year + age ## equivalent
f <- claim ~ contract*age + contract:age:year ## equivalent

attributes(terms(lm(f, data = hach)))$factors

## more generally, n factors, factor1 = contract, factorn = year
## X_1...n = m
##           + Y1 + ... + Y(n-1)
##           + Y12 + ... + Y1(n-1) + ... + + Y23 + ... + Y2(n-1)
##           .
##           .
##           .
##           + Y_1..(n-1)
##           + Y_1..n
##########
########## claim ~ (. - year) ^ (n-1) + contract:factor2:factor3:...:year



#### Hierarchical

## Xpijt = m + Yp + Ypi + Ypi + Ypij + Ypijt

claim ~ sector/unit/contract/year

f <- claim ~ unit/contract/year

attributes(terms(lm(f, data = hach)))$factors


