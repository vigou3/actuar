form1 <- ~ x1 + x1:x2 | a1/a2/a3
form <- form1[[length(form1)]] ## class 'call', devient de longueur '3'

asOneSidedFormula(form[[3]]) ## côté droit

asOneSidedFormula(form[[2]]) ## côté gauche

sep <- "/"
val <- splitFormula(asOneSidedFormula(form[[2]]), sep = sep)
names(val) <- unlist(lapply(val, function(el) deparse(el[[2]])))
as.formula(eval(parse(text = paste("~", paste(names(val), collapse = "/")))))

getGroupsFormula(form1)
getGroupsFormula(form1, asList = TRUE)
getCovariateFormula(form1)

library(actuar)
data(hachemeister) ; hach <- hachemeister$claims

dframepf <- function(x, contracts, subs = NULL)
{
    d <- dim(x)
    nrows <- d[1]
    ir <- seq(length = nrows)
    ncols <- d[2] + 1 + length(subs)
    ic <- seq(length = d[2])
    collabs <- c(names(subs), "contract", paste("Y", ic, sep = ""))
    value <- vector("list", ncols)
    for (i in (ic)) value[[i+1+length(subs)]] <- as.vector(x[, i])
    if (!is.null(subs))
        value[1:length(subs)] <- subs
    value[[length(subs)+1]] <- as.factor(contracts)
    row.names <- as.character(ir)
    names(value) <- collabs
    attr(value, "row.names") <- row.names
    class(value) <- "data.frame"
    value
}
x <- dframepf(hach, contract = c("4B12", "66FZ", "76H5", "F44R", "FGG5"), list(unit = c(1,1,2,3,2)))
cm(formula = ~ contract | Y1 + Y2 + Y3, data = x)
cm(formula = ~ contract | Y1 + Y2 + Y3)


cm(formula = ~ contract | ., data = x)
cm(formula = ~ contract | Y1+Y2+Y3+Y4+Y5+Y6+Y7+Y8+Y9+Y10+Y11+Y12, data = x)
<<<<<<< .mine
cm(formula = ~ contract | ., data = x, subset = (unit == 1))



=======
cm(formula = ~ contract | ., data = x, subset = (x$unit == 1))


######## VG #########

### On peut extraire les éléments du data frame selon la
### classification spécifiée dans la partie gauche de la formule avec
### tapply() ou, plus pratique, avec by().

## Un data frame d'exemple. Classification par secteur (s), unité (u)
## et contrat (c).
dd <- data.frame(s = as.factor(rep(1:2, c(9, 6))),
                 u = as.factor(rep(1:5, each = 3)),
                 c = as.factor(1:15),
                 matrix(sample(1:100, 75, rep = TRUE), 15))
names(dd)[4:8] <- paste("Y", 1:5, sep = "")

## Supposons une classification par secteur et unité seulement.
form <- ~ s + u | Y1 + Y2 + Y3 + Y4 + Y5
( struct <- all.vars(form[[2]][[2]]) )
( years <- all.vars(form[[2]][[3]]) )

## Extraction des données.
( X <- by(dd, dd[struct], subset, select = years) )

## Données du secteur 1, unité 2.
X[[1, 2]]

## Classification par secteur, unité et contrat (hiérarchique).
form <- ~ s + s:u + s:u:c| Y1 + Y2 + Y3 + Y4 + Y5
( struct <- all.vars(form[[2]][[2]]) )
( years <- all.vars(form[[2]][[3]]) )

## Extraction des données.
( X <- by(dd, dd[struct], subset, select = years) )

## Données du contrat 4, unité 2, secteur 1.
X[[1, 2, 4]]

## Classification par contrat seulement (Bühlmann-Straub).
form <- ~ c| Y1 + Y2 + Y3 + Y4 + Y5
( struct <- all.vars(form[[2]][[2]]) )
( years <- all.vars(form[[2]][[3]]) )

## Extraction des données.
( X <- by(dd, dd[struct], subset, select = years) )

## Données du contrat 4.
X[[4]]


### Bon, le tout n'est pas complet:
###
### 1. dans le dernier exemple, nous voulons avoir les données sous
###    forme de matrice, pas sous forme de liste. Mais est-ce
###    absolument nécessaire? On pourrait adapter bstraub() pour
###    supposer que les données de chaque contrat se trouvent dans un
###    élément d'une liste.
###
### 2. la façon de définir 'struct', ci-dessus, ne tient pas vraiment
###    compte des interactions entre les variables. Toutefois, je crois
###    que l'on peut extraire cette information de 'terms(form[[2]][2]]).
###
### Néanmoins, je crois être sur la bonne piste.
###
### ...
###
### Ça, c'était à 17h30. À 2h30, j'ai presque une fonction pour créer
### une liste contenant les données pour n'importe quel type de
### modèle. Je l'ai appelée 'model.list', c'est bien, hein? La
### fonction ne marche pas comme c'est là parce qu'il y a des
### problèmes d'environnment, je crois. Mais on se rapproche drôlement.
###
### Là, je vais me coucher...

model.list <- function(formu, data)
{
    form <- formu[[2]]

    ## colonnes des données
    rhs <- form[[3]]
    years <- all.vars(rhs)

    ## termes de la structure
    lhs <- terms(asOneSidedFormula(form[[2]]))
    vars <- attr(lhs, "variables")
    varnames <- sapply(vars, deparse, width.cutoff = 500)[-1]
    labs <- attr(lhs, "term.labels")
    ord <- attr(lhs, "order")

    ## il faut envelopper les données dans une liste avant d'entrer
    ## dans la boucle
    ml <- list(data)

    fun <- function(x, cols, select)
    {
        INDICES <- x[, cols]
        if (is.list(INDICES))
            IND <- lapply(INDICES, factor)
        else
            IND <- INDICES
        tapply(1:nrow(x), IND, function(y) subset(x[y, ], select = select))
    }

    for (i in seq(to = max(ord)))
    {
        todrop <- all.vars(parse(text = labs[ord == i]))
        tokeep <- setdiff(varnames, todrop)
        ml <- lapply(ml, fun, cols = todrop, select = c(tokeep, years))
    }
    ml
}

debug(model.list)
model.list(~ s + s:u + s:u:c|  Y1 + Y2 + Y3 + Y4 + Y5, dd)



### test
f <- function(x, IND, sel)
    by(x, IND, subset, select = substitute(sel))

f(dd, dd[, 1], c("u", "Y1"))
>>>>>>> .r459
