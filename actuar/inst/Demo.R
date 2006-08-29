## Louis-Philippe Pouliot, 29/08/06
##
#######                                              #######
#######    Démontration des fonctions programmées    #######
#######    et modifiées                              #######
#######                                              ####### 



                      ### Simulation de portefeuilles

#### Fonction simpf()

## Création de l'objet
## Portfolio where both frequency and severity models are mixed.
     modelfreq <- list(dist1 = "pois",
                       par1 = list(lambda = quote(Lambda * weights)),
                       dist2 = "gamma",
                       par2 = c(shape = 2, rate = 1))
     modelsev<-list(dist1 = "lnorm",
                    par1 = list(meanlog = quote(Theta), sdlog = 1),
                    dist2 = "norm",
                    par2 = c(mean = 5, sd = 1))
data(hachemeister)
weights <- hachemeister$weights[,1:8]/mean(hachemeister$weights)
(x <- simpf(5, 8, modelfreq, modelsev, weights))


## Traitement de l'objet
summary(x)                      ## Chaque réclamation, pour chaque contrat, chaque année

aggregate(x, by = c("contract", "year"), sum)      ## Sommaires par contrat, par année ou tout en même temps

aggregate(x, by = c("contract", "year"), length)
frequency(x)                                       ## même chose



severity(x, bycol = FALSE, y.exclude = 1)



                   #####  Théorie du risque  #####

     ## Normal Power approximation
       res <- aggregateDist("np2", moments = c(200, 200, 0.5))
       res(210)

     ## Simulation method
       model.freq <- list(dist = "pois", par = list(lambda = 4))
       model.sev <- list(dist = "gamma", par = list(shape = 100, rate = 2))
       res <- aggregateDist("simulation", model.sev = model.sev,
                            model.freq = model.freq, n = 10000)
       mean(res)
       res(250)

     ## Approximation by the recursive method
       # with a parameterized distribution severity model

       model.sev <- list(dist = "gamma", par = list(shape = 100, rate = 2))
       (res <- aggregateDist("recursive", model.sev, model.freq,
                            x.scale = 1))

       # with a vectorized severity model

       model.sev <- c(0.1, 0.1, 0.2, 0.3, 0.2, 0.1)
       (res <- aggregateDist("recursive", model.sev, model.freq,
                            x.scale = 100))
       plot(res)

     ## Exact calculation

       model.freq <- c(0.1, 0.1, 0.2, 0.3, 0.2, 0.1)
       (res <- aggregateDist("exact", model.sev = model.sev, model.freq = model.freq, x.scale = 100))
       plot(res)
       res(1000)

                   #####  Crédibilité  #####

data(hachemeister)
(x <- bstraub(hachemeister$claims, hachemeister$weights))  ## affichage
summary(x)

## Build a 'data' object
data(hachemeister)
contract <- c("C1", "C2" ,"C3", "C4", "C5")
sector <- LETTERS[c(1,2,1,2,2)]
(x <- data.frame(sector, contract, hachemeister$claims,
                 hachemeister$weights))

ghcm(~sector + sector:contract, x, X1:X12, X1.1:X12.1)

ghcm(~contract, x, X1:X12, X1.1:X12.1, subset = sector == "B")

ghcm(~contract, x, X1:X12, X1.1:X12.1, subset = sector %in% c("A","B"))

(z <- ghcm(~sector + sector:contract, x, X1:X12, X1.1:X12.1))

predict(z)

