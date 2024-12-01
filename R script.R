library(factoextra)
library(FactoMineR)
library(raster)
library(landscapemetrics)
library(rgdal)
library(foreach)
library(reshape2)
library(factoextra)
library(FactoMineR)
library(terra)
library(MASS)
library(car)
library(interplot)
library(emmeans)
library(ggeffects)
library(nFactors)
library(MuMIn)
library(modEvA)
library(ggplot2)
library(ggthemes)
library(partR2)
library(r2glmm)
library(effectsize)
library(vegan)
library(reshape2)

mon_theme <- theme(panel.grid.major = element_blank(),
                   panel.grid.minor = element_blank(),
                   panel.background = element_rect(fill = "white"),
                   panel.border = element_rect(colour = "black", 
                                               fill = NA),
                   plot.title = element_text(size = 18, 
                                             face = "bold", 
                                             hjust = 0.5),
                   plot.caption = element_text(size = 10,
                                               face = "italic", 
                                               color = "grey"),
                   legend.title = element_text(size = 12,
                                               face = "bold"),
                   legend.text = element_text(size = 12),
                   axis.title = element_text(size = 14, 
                                             face = "bold"),
                   axis.text = element_text(size = 12))

#Calling up the data set
A <- read.csv("df.csv", sep=";") 

##Calculating abundance and diversity variables
#Total Abundance
abondance_tot <- vector(length = NROW(A))
for (i in c(1:NROW(A))){
  abondance_tot[i] <- sum(A[i,3:16])
}

#Pielou index for total arthropods, carabids, and spiders
pielou.arth <- vector(length = NROW(A))
pielou.car <- vector(length = NROW(A))
pielou.ara <- vector(length = NROW(A))

for (i in c(1:NROW(A))){
  pielou.arth[i] <- A$shannon[i]/log(length(which(A[i,3:16]>0)))
  pielou.car[i] <- diversity(A[i,110:142], index = "shannon")/log(length(which(A[i,110:142]>0)))
  pielou.ara[i] <- diversity(A[i,73:107], index = "shannon")/log(length(which(A[i,73:107]>0)))
}

df <- data.frame(A,abondance_tot, pielou.arth, pielou.car, pielou.ara)


####Models for abundance####
##Spring
#First step: Local variables
ab_log <- log(df2$abondance_tot+1)
df2 <- data.frame(ab_log, df2)

mod <- lm(ab_log ~ intensity + ndvi + microclimat, 
           data = df2[which(df2$saison == "printemps"),])
vif(mod)
stepAIC(mod, direction = "backward")
mod <- lm(ab_log ~ microclimat,data = df2[which(df2$saison == 
                                                 "printemps"),])
summary(mod)

#Second step: land cover variables (composition)
mod.occ.sol <- lm(ab_log ~ microclimat + X._prairies_500 + diversite_500 + grain_bocager_500m, 
            data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                        scope = list(upper = mod.occ.sol,
                                     lower = stepAIC(mod))))


#Third step: microlimate and land cover variables (configuration)
fullModel2 <- lm(ab_log ~ microclimat + ndvi_500 + Taille_parcelle_500, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2, direction = "backward",
        scope = list(upper = fullModel2,
                     lower = minModel2)))


#Last step : farming practices in landscape variables
fullModel3 <- lm(ab_log ~ microclimat + ai_l + moy.intensity.tot + enn_c_1, 
                  data = df2[which(df2$saison == "printemps"),]) 
fullModel3b <- lm(ab_log ~ (microclimat + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "printemps"),])
minModel3 <- stepAIC(fullModel2, direction = "backward",
                     scope = list(upper = fullModel2,
                                  lower = minModel2))


summary(modF <- stepAIC(fullModel3, direction = "backward",
        scope = list(upper = fullModel3,
                     lower = minModel2)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = ab_log ~ microclimat, data = df2[which(df2$saison == "printemps"),])

B <- lm(formula = ab_log ~ 1 , data = df2[which(df2$saison == "printemps"),])

C <- lm(formula = ab_log ~ 1, data = df2[which(df2$saison == "printemps"),])

AB <- lm(formula = ab_log ~ microclimat , data = df2[which(df2$saison == "printemps"),])

AC <- lm(formula = ab_log ~ microclimat, data = df2[which(df2$saison == "printemps"),])

BC <- lm(formula = ab_log ~ 1, data = df2[which(df2$saison == "printemps"),])

ABC <- lm(formula = ab_log ~ microclimat, data = df2[which(df2$saison == "printemps"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

##Autumn
#variables locales 
ab_log <- log(df2$abondance_tot+1)
df2 <- data.frame(ab_log, df2)

mod <- lm(ab_log ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
stepAIC(mod, direction = "backward")
mod <- lm(ab_log ~ microclimat,data = df2[which(df2$saison == 
                                                  "automne"),])
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(ab_log ~ microclimat + X._prairies_500 + diversite_500 + grain_bocager_500m, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(ab_log ~ microclimat + ndvi_500 + Taille_parcelle_500, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(modF <- stepAIC(fullModel2, direction = "backward",
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(ab_log ~ microclimat + ndvi_500 + Taille_parcelle_500 + 
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
fullModel3b <- lm(ab_log ~ (microclimat + ndvi_500 + Taille_parcelle_500 + 
                   ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel3 <- stepAIC(fullModel2, direction = "backward",
                     scope = list(upper = fullModel2,
                                  lower = minModel2))


summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel2)))


r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = ab_log ~ microclimat, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = ab_log ~ ndvi_500 + Taille_parcelle_500 , data = df2[which(df2$saison == "automne"),])

C <- lm(formula = ab_log ~ moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = ab_log ~ microclimat + ndvi_500 + Taille_parcelle_500 , data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = ab_log ~ microclimat + moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = ab_log ~ ndvi_500 + Taille_parcelle_500 + moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = ab_log ~ microclimat+ndvi_500 + Taille_parcelle_500 + moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

####Richness####
##Spring
#variables locales
mod <- lm(richesse ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "printemps"),])#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(richesse ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(richesse ~ ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(richesse ~ ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "printemps"),]) 
fullModel3b <- lm(richesse ~ (ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))

r2beta(modF, method = "lm")


##Autumn
#variables locales
mod <- lm(richesse ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(richesse ~ intensity + ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(richesse ~ intensity + ndvi + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(richesse ~ intensity + ndvi + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
fullModel3b <- lm(richesse ~ (intensity + ndvi + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))

r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = richesse ~ intensity + ndvi, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = richesse ~ 1 , data = df2[which(df2$saison == "automne"),])

C <- lm(formula = richesse ~ moy.intensity.tot + ai_l, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = richesse ~ intensity + ndvi, data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = richesse ~ intensity + ndvi+ moy.intensity.tot + ai_l, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = richesse ~ moy.intensity.tot + ai_l, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = richesse ~ intensity + ndvi+ moy.intensity.tot + ai_l, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


####Shannon index####
##Spring
#variables locales
mod <- lm(pielou.arth ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "printemps"),])#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(pielou.arth ~ intensity + microclimat + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(pielou.arth ~ intensity + microclimat + Taille_parcelle_500, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(pielou.arth ~ intensity + microclimat + Taille_parcelle_500+
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = pielou.arth ~ intensity + microclimat, data = df2[which(df2$saison == "printemps"),])

B <- lm(formula = pielou.arth ~ Taille_parcelle_500 , data = df2[which(df2$saison == "printemps"),])

C <- lm(formula = pielou.arth ~ 1, data = df2[which(df2$saison == "printemps"),])

AB <- lm(formula = pielou.arth ~ intensity + microclimat + Taille_parcelle_500 , data = df2[which(df2$saison == "printemps"),])

AC <- lm(formula = pielou.arth ~ intensity + microclimat, data = df2[which(df2$saison == "printemps"),])

BC <- lm(formula = pielou.arth ~ Taille_parcelle_500, data = df2[which(df2$saison == "printemps"),])

ABC <- lm(formula = pielou.arth ~ intensity + microclimat + Taille_parcelle_500, data = df2[which(df2$saison == "printemps"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)



##Autumn
#variables locales
mod <- lm(pielou.arth ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(pielou.arth ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(pielou.arth ~ ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(pielou.arth ~ ndvi_500 +
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = pielou.arth ~ 1, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = pielou.arth ~ ndvi_500, data = df2[which(df2$saison == "automne"),])

C <- lm(formula = pielou.arth ~enn_c_1, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = pielou.arth ~ ndvi_500, data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = pielou.arth ~ enn_c_1, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = pielou.arth ~ ndvi_500+ enn_c_1, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = pielou.arth ~ndvi_500+ enn_c_1, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


####cComposition Spring####
arth <- df[which(df2$saison == "printemps"),c(3:16)]

res.pca <- PCA(log(arth+1), scale.unit = T, graph = T)
fviz_pca_var(res.pca)+mon_theme
fviz_eig(res.pca)
nCng(res.pca$eig[,1])
dimdesc(res.pca)

pos.arth <- res.pca$ind$coord[,1:3]
colnames(pos.arth) = c("arth.1","arth.2","arth.3")

df2.arth <- data.frame(df2[which(df2$saison == "printemps"),], pos.arth)

##Dim.1
#variables locales
mod <- lm(arth.1 ~ intensity + ndvi + microclimat, 
          data = df2.arth)#transformer en log

mod <- lm(arth.1 ~ ai_l, 
          data = df2.arth)

vif(mod)
summary(stepAIC(mod))
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(arth.1 ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.arth)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(arth.1 ~ ndvi + X._prairies_500 +  Taille_parcelle_500+ 
                   ndvi_500 + grain_bocager_500m, 
                 data = df2.arth) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))
vif(fullModel2)

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))
vif(stepAIC(fullModel2,
            scope = list(upper = fullModel2,
                         lower = minModel2)))

#4 <- pratiques
fullModel3b <- lm(arth.1 ~ (ndvi + X._prairies_500 +  Taille_parcelle_500+ 
                               ndvi_500 + ai_l + enn_c_1) * moy.intensity.tot, 
                  data = df2.arth) 
fullModel3 <- lm(arth.1 ~ ndvi + X._prairies_500 +  Taille_parcelle_500+ 
                              ndvi_500 + ai_l + enn_c_1 + moy.intensity.tot, 
                  data = df2.arth) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))


summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = arth.1 ~ ndvi + X._prairies_500 + Taille_parcelle_500 + ndvi_500, data = df2.arth)

B <- lm(formula = arth.1 ~ ai_l ,data = df2.arth)

AB <- lm(formula = arth.1 ~ ndvi + X._prairies_500 + Taille_parcelle_500 + ndvi_500 + ai_l, data = df2.arth)

varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1],
        AB =  r.squaredGLMM(AB)[1,1],
        A.name = "Autres", B.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

A <- lm(formula = arth.1 ~ ndvi, data = df2.arth)

B <- lm(formula = arth.1 ~ X._prairies_500 + Taille_parcelle_500 + ndvi_500 , data = df2.arth)

C <- lm(formula = arth.1 ~ ai_l, data = df2.arth)

AB <- lm(formula = arth.1 ~ ndvi + X._prairies_500 + Taille_parcelle_500 + ndvi_500 , data = df2.arth)

AC <- lm(formula = arth.1 ~ ndvi+ai_l, data = df2.arth)

BC <- lm(formula = arth.1 ~ ai_l+X._prairies_500 + Taille_parcelle_500 + ndvi_500, data = df2.arth)

ABC <- lm(formula = arth.1 ~ ndvi + X._prairies_500 + Taille_parcelle_500 + ndvi_500+ai_l, data = df2.arth)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


interplot(modF, var1 = "moy.intensity.tot", var2 = "enn_c_1") + 
  ylab("") +
  xlab("")+
  geom_hline(yintercept = 0, linetype = "dashed") + geom_line(size=1) +
  theme_classic()

r2beta(modF, method = "lm")
omega_squared(modF)


##Dim.2
#variables locales
mod <- lm(arth.2 ~ intensity + ndvi + microclimat, 
          data = df2.arth)#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(arth.2 ~ intensity + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.arth)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(arth.2 ~ intensity + 
                   ndvi_500 + grain_bocager_500m, 
                 data = df2.arth) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))

#4 <- pratiques
fullModel3b <- lm(arth.2 ~ (intensity + ndvi_500 + ai_l + enn_c_1) * moy.intensity.tot, 
                  data = df2.arth) 
fullModel3 <- lm(arth.2 ~ intensity + ndvi_500 + ai_l + enn_c_1 + moy.intensity.tot, 
                  data = df2.arth) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = arth.2 ~ intensity, data = df2.arth)

B <- lm(formula = arth.2 ~ ndvi_500 , data = df2.arth)

C <- lm(formula = arth.2 ~ 1, data = df2.arth)

AB <- lm(formula = arth.2 ~ intensity + ndvi_500 , data = df2.arth)

AC <- lm(formula = arth.2 ~ intensity, data = df2.arth)

BC <- lm(formula = arth.2 ~ intensity, data = df2.arth)

ABC <- lm(formula = arth.2 ~ intensity + ndvi_500, data = df2.arth)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

##Dim.3
#variables locales
mod <- lm(arth.3 ~ intensity + ndvi + microclimat, 
          data = df2.arth)#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(arth.3 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.arth)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(arth.3 ~ ndvi_500 + grain_bocager_500m, 
                 data = df2.arth) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))

#4 <- pratiques
fullModel3b <- lm(arth.3 ~ (ai_l + enn_c_1) * moy.intensity.tot, 
                  data = df2.arth) 
fullModel3 <- lm(arth.3 ~ ai_l + enn_c_1 + moy.intensity.tot, 
                  data = df2.arth)
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = arth.3 ~ 1, data = df2.arth)

B <- lm(formula = arth.3 ~ 1 , data = df2.arth)

C <- lm(formula = arth.3 ~ enn_c_1 + moy.intensity.tot, data = df2.arth)

AB <- lm(formula = arth.3 ~ 1 , data = df2.arth)

AC <- lm(formula = arth.3 ~ enn_c_1 + moy.intensity.tot, data = df2.arth)

BC <- lm(formula = arth.3 ~ enn_c_1 + moy.intensity.tot, data = df2.arth)

ABC <- lm(formula = arth.3 ~ enn_c_1 + moy.intensity.tot, data = df2.arth)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


####Composition autumn####
arth <- df[which(df2$saison == "automne"),c(3:16)]

somme <- vector(length = NCOL(arth))
for (i in c(1:NCOL(arth))){
  somme[i] <- sum(arth[,i])  
}
arth <- arth[,which(somme>0)]

res.pca <- PCA(log(arth+1), scale.unit = T, graph = T)
fviz_pca_var(res.pca, repel = T)+mon_theme
fviz_eig(res.pca)
nCng(res.pca$eig[,1])
dimdesc(res.pca)

pos.arth <- res.pca$ind$coord[,1:3]
colnames(pos.arth) = c("arth.1","arth.2","arth.3")

df2.arth <- data.frame(df2[which(df2$saison == "automne"),], pos.arth)

##Dim.1
#variables locales
mod <- lm(arth.1 ~ intensity + ndvi + microclimat, 
          data = df2.arth)#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(arth.1 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.arth)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(arth.1 ~ ndvi_500 + grain_bocager_500m, 
                 data = df2.arth) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3b <- lm(arth.1 ~ (ndvi_500 + grain_bocager_500m + ai_l + enn_c_1) * moy.intensity.tot, 
                  data = df2.arth)
fullModel3 <- lm(arth.1 ~ ndvi_500 + grain_bocager_500m + ai_l + enn_c_1 + moy.intensity.tot, 
                  data = df2.arth) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = arth.1 ~ 1, data = df2.arth)

B <- lm(formula = arth.1 ~ ndvi_500 + grain_bocager_500m, data = df2.arth)

C <- lm(formula = arth.1 ~ enn_c_1 + moy.intensity.tot, data = df2.arth)

AB <- lm(formula = arth.1 ~ ndvi_500 + grain_bocager_500m, data = df2.arth)

AC <- lm(formula = arth.1 ~ ndvi + enn_c_1 + moy.intensity.tot, data = df2.arth)

BC <- lm(formula = arth.1 ~ ndvi_500 + grain_bocager_500m+enn_c_1 + moy.intensity.tot, data = df2.arth)

ABC <- lm(formula = arth.1 ~ndvi_500 + grain_bocager_500m+enn_c_1 + moy.intensity.tot, data = df2.arth)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

##Dim.2
#variables locales
mod <- lm(arth.2 ~ intensity + ndvi + microclimat, 
          data = df2.arth)#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(arth.2 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.arth)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(arth.2 ~ X._prairies_500 + diversite_500 + 
                   ndvi_500 + grain_bocager_500m, 
                 data = df2.arth) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))

#4 <- pratiques
fullModel3b <- lm(arth.2 ~ (X._prairies_500 + diversite_500 + ai_l + enn_c_1) * moy.intensity.tot, 
                  data = df2.arth) 
fullModel3 <- lm(arth.2 ~ X._prairies_500 + diversite_500 + ai_l + enn_c_1 + moy.intensity.tot, 
                  data = df2.arth) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)
interplot(modF, var1 = "moy.intensity.tot", var2 = "X._prairies_500") + 
  ylab("") +
  xlab("")+
  geom_hline(yintercept = 0, linetype = "dashed") + geom_line(size=1) +
  theme_classic()

A <- lm(formula = arth.2 ~ 1, data = df2.arth)

B <- lm(formula = arth.2 ~X._prairies_500 + diversite_500, data = df2.arth)

C <- lm(formula = arth.2 ~ 1, data = df2.arth)

AB <- lm(formula = arth.2 ~ X._prairies_500 + diversite_500, data = df2.arth)

AC <- lm(formula = arth.2 ~ 1, data = df2.arth)

BC <- lm(formula = arth.2 ~X._prairies_500 + diversite_500, data = df2.arth)

ABC <- lm(formula = arth.2 ~ X._prairies_500 + diversite_500, data = df2.arth)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)
##Dim.3
#variables locales
mod <- lm(arth.3 ~ intensity + ndvi + microclimat, 
          data = df2.arth)#transformer en log
vif(mod)
summary(stepAIC(mod))

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(arth.3 ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.arth)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(arth.3 ~ ndvi + ndvi_500 + grain_bocager_500m, 
                 data = df2.arth) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))

#4 <- pratiques
fullModel3b <- lm(arth.3 ~ (ndvi + ai_l + enn_c_1) * moy.intensity.tot, 
                  data = df2.arth) 
fullModel3 <- lm(arth.3 ~ ndvi + ai_l + enn_c_1 + moy.intensity.tot, 
                  data = df2.arth) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = arth.3 ~ ndvi, data = df2.arth)

B <- lm(formula = arth.3 ~1, data = df2.arth)

C <- lm(formula = arth.3 ~ enn_c_1, data = df2.arth)

AB <- lm(formula = arth.3 ~ ndvi, data = df2.arth)

AC <- lm(formula = arth.3 ~ ndvi+enn_c_1, data = df2.arth)

BC <- lm(formula = arth.3 ~enn_c_1, data = df2.arth)

ABC <- lm(formula = arth.3 ~ ndvi+enn_c_1, data = df2.arth)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


####Spiders abundance spring####
#variables locales
mod <- lm(log(A.ab+1) ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "printemps"),])#transformer en log
vif(mod)
stepAIC(mod)
mod <- lm(log(A.ab+1) ~ ndvi, data = df2[which(df2$saison == 
                                                 "printemps"),])
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(log(A.ab+1) ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(log(A.ab+1) ~ ndvi + X._prairies_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(log(A.ab+1) ~ ndvi + X._prairies_500 + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "printemps"),]) 
fullModel3b <- lm(log(A.ab+1) ~ (ndvi + X._prairies_500 + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)


A <- lm(formula = log(A.ab+1) ~ ndvi, data = df2[which(df2$saison == "printemps"),])

B <- lm(formula = log(A.ab+1) ~ X._prairies_500 , data = df2[which(df2$saison == "printemps"),])

C <- lm(formula = log(A.ab+1) ~ 1, data = df2[which(df2$saison == "printemps"),])

AB <- lm(formula = log(A.ab+1) ~ ndvi + X._prairies_500 , data = df2[which(df2$saison == "printemps"),])

AC <- lm(formula = log(A.ab+1) ~ ndvi, data = df2[which(df2$saison == "printemps"),])

BC <- lm(formula = log(A.ab+1) ~ X._prairies_500, data = df2[which(df2$saison == "printemps"),])

ABC <- lm(formula = log(A.ab+1) ~ ndvi + X._prairies_500, data = df2[which(df2$saison == "printemps"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

########Spiders richness spring########
#variables locales
mod <- lm(A.rs ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "printemps"),])#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(A.rs ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
        scope = list(upper = mod.occ.sol,
                     lower = stepAIC(mod))))

plot(df$A.rs[which(df2$saison == "printemps")] ~ df$diversite_500[which(df2$saison == "printemps")])
plot(df$A.rs[which(df2$saison == "printemps")] ~ df$Taille_parcelle_500[which(df2$saison == "printemps")])

#3
fullModel2 <- lm(A.rs ~ diversite_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(A.rs ~ diversite_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m 
                 + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "printemps"),])

fullModel3b <- lm(A.rs ~ (diversite_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m 
                 + ai_l + enn_c_1 ) * moy.intensity.tot , 
                 data = df2[which(df2$saison == "printemps"),]) 

minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = A.rs ~ 1, data = df2[which(df2$saison == "printemps"),])

B <- lm(formula = A.rs ~ diversite_500 + Taille_parcelle_500 , data = df2[which(df2$saison == "printemps"),])

C <- lm(formula = A.rs ~ moy.intensity.tot, data = df2[which(df2$saison == "printemps"),])

AB <- lm(formula = A.rs ~ diversite_500 + Taille_parcelle_500 , data = df2[which(df2$saison == "printemps"),])

AC <- lm(formula = A.rs ~ moy.intensity.tot, data = df2[which(df2$saison == "printemps"),])

BC <- lm(formula = A.rs ~ moy.intensity.tot + diversite_500 + Taille_parcelle_500, data = df2[which(df2$saison == "printemps"),])

ABC <- lm(formula = A.rs ~ moy.intensity.tot + diversite_500 + Taille_parcelle_500, data = df2[which(df2$saison == "printemps"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

####Spiders pielou index spring####
#variables locales
mod <- lm(pielou.ara ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "printemps"),])#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(pielou.ara ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))

plot(df$A.rs[which(df2$saison == "printemps")] ~ df$diversite_500[which(df2$saison == "printemps")])
plot(df$A.rs[which(df2$saison == "printemps")] ~ df$Taille_parcelle_500[which(df2$saison == "printemps")])

#3
fullModel2 <- lm(pielou.ara ~ ndvi + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(pielou.ara ~ ndvi + Taille_parcelle_500 + 
                 ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "printemps"),])

fullModel3b <- lm(A.rs ~ (diversite_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m 
                          + ai_l + enn_c_1 ) * moy.intensity.tot , 
                  data = df2[which(df2$saison == "printemps"),]) 

minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                        scope = list(upper = fullModel3,
                                     lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = pielou.ara ~ ndvi, data = df2[which(df2$saison == "printemps"),])

B <- lm(formula = pielou.ara ~ Taille_parcelle_500 , data = df2[which(df2$saison == "printemps"),])

C <- lm(formula = pielou.ara ~ 1, data = df2[which(df2$saison == "printemps"),])

AB <- lm(formula = pielou.ara ~ ndvi + Taille_parcelle_500 , data = df2[which(df2$saison == "printemps"),])

AC <- lm(formula = pielou.ara ~ ndvi, data = df2[which(df2$saison == "printemps"),])

BC <- lm(formula = pielou.ara ~ Taille_parcelle_500, data = df2[which(df2$saison == "printemps"),])

ABC <- lm(formula = pielou.ara ~ ndvi + Taille_parcelle_500, data = df2[which(df2$saison == "printemps"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

####Spiders composition spring####
##acp
ara <- df[which(df2$saison == "printemps"),c(73:106)]

somme <- vector(length = NCOL(ara))
for (i in c(1:NCOL(ara))){
  somme[i] <- sum(ara[,i])  
}
ara <- ara[,which(somme>0)]

res.pca <- PCA(log(ara+1), scale.unit = T, graph = T)
fviz_pca_var(res.pca, repel = T)+mon_theme
fviz_eig(res.pca)
nCng(res.pca$eig[,1])
dimdesc(res.pca)

pos.ara <- res.pca$ind$coord[,1:3]
colnames(pos.ara) = c("ara.1","ara.2","ara.3")

df2.ara <- data.frame(df2[which(df2$saison == "printemps"),],pos.ara)

##Dim1
#variables locales
mod <- lm(ara.1 ~ intensity + ndvi + microclimat, 
          data = df2.ara)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(ara.1 ~ microclimat + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.ara)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(ara.1 ~ microclimat + ndvi_500 + grain_bocager_500m, 
                 data = df2.ara) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(ara.1 ~ microclimat + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.ara) 
fullModel3b <- lm(ara.1 ~ (microclimat + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.ara) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = ara.1 ~ microclimat, data = df2.ara)

B <- lm(formula = ara.1 ~ 1 , data = df2.ara)

C <- lm(formula = ara.1 ~ 1, data = df2.ara)

AB <- lm(formula = ara.1 ~ microclimat , data = df2.ara)

AC <- lm(formula = ara.1 ~ microclimat, data = df2.ara)

BC <- lm(formula = ara.1 ~ 1, data = df2.ara)

ABC <- lm(formula = ara.1 ~ microclimat, data = df2.ara)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


##Dim2
#variables locales
mod <- lm(ara.2 ~ intensity + ndvi + microclimat, 
          data = df2.ara)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(ara.2 ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.ara)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(ara.2 ~ ndvi + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2.ara) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(ara.2 ~ ndvi + Taille_parcelle_500 + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.ara)
fullModel3b <- lm(ara.2 ~ (ndvi + Taille_parcelle_500 + ai_l  + enn_c_1) * moy.intensity.tot, 
                 data = df2.ara) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = ara.2 ~ ndvi, data = df2.ara)

B <- lm(formula = ara.2 ~ Taille_parcelle_500 , data = df2.ara)

C <- lm(formula = ara.2 ~ moy.intensity.tot , data = df2.ara)

AB <- lm(formula = ara.2 ~ ndvi + Taille_parcelle_500  , data = df2.ara)

AC <- lm(formula = ara.2 ~ ndvi + moy.intensity.tot, data = df2.ara)

BC <- lm(formula = ara.2 ~ moy.intensity.tot +Taille_parcelle_500 , data = df2.ara)

ABC <- lm(formula = ara.2 ~ moy.intensity.tot +Taille_parcelle_500+ndvi, data = df2.ara)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

##Dim3
#variables locales
mod <- lm(ara.3 ~ intensity + ndvi + microclimat, 
          data = df2.ara)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(ara.3 ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.ara)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(ara.3 ~ ndvi + ndvi_500 + grain_bocager_500m, 
                 data = df2.ara) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(ara.3 ~ ndvi + grain_bocager_500m + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.ara) 
fullModel3b <- lm(ara.3 ~ (ndvi + grain_bocager_500m + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.ara)
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = ara.3 ~ ndvi, data = df2.ara)

B <- lm(formula = ara.3 ~ grain_bocager_500m , data = df2.ara)

C <- lm(formula = ara.3 ~ moy.intensity.tot , data = df2.ara)

AB <- lm(formula = ara.3 ~ ndvi + grain_bocager_500m  , data = df2.ara)

AC <- lm(formula = ara.3 ~ ndvi + moy.intensity.tot, data = df2.ara)

BC <- lm(formula = ara.3 ~ moy.intensity.tot +grain_bocager_500m , data = df2.ara)

ABC <- lm(formula = ara.3 ~ moy.intensity.tot +grain_bocager_500m+ndvi, data = df2.ara)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

####Carabid abundance spring####
#variables locales
mod <- lm(log(C.ab+1) ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "printemps"),])#transformer en log
vif(mod)
stepAIC(mod)
mod <- lm(log(C.ab+1) ~ intensity + ndvi, data = df2[which(df2$saison == 
                                                 "printemps"),])
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(log(C.ab+1) ~ intensity + ndvi + X._prairies_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(log(C.ab+1) ~ intensity + ndvi + X._prairies_500 +
                   ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(log(C.ab+1) ~ intensity + ndvi + X._prairies_500 + 
                   ndvi_500 + grain_bocager_500m + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "printemps"),]) 
fullModel3b <- lm(log(C.ab+1) ~ (intensity + ndvi + X._prairies_500 + 
                   ndvi_500 + grain_bocager_500m + ai_l + enn_c_1) * moy.intensity.tot , 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = log(C.ab+1) ~ intensity + ndvi, data = df2[which(df2$saison == "printemps"),])

B <- lm(formula = log(C.ab+1) ~ X._prairies_500 + ndvi_500 + grain_bocager_500m , data = df2[which(df2$saison == "printemps"),])

C <- lm(formula = log(C.ab+1) ~ ai_l+moy.intensity.tot, data = df2[which(df2$saison == "printemps"),])

AB <- lm(formula = log(C.ab+1) ~ intensity + ndvi + X._prairies_500 + ndvi_500 + grain_bocager_500m, data = df2[which(df2$saison == "printemps"),])

AC <- lm(formula = log(C.ab+1) ~ intensity + ndvi + ai_l+moy.intensity.tot, data = df2[which(df2$saison == "printemps"),])

BC <- lm(formula = log(C.ab+1) ~ ai_l+moy.intensity.tot + X._prairies_500 + ndvi_500 + grain_bocager_500m, data = df2[which(df2$saison == "printemps"),])

ABC <- lm(formula = log(C.ab+1) ~ intensity + ndvi + X._prairies_500 + ndvi_500 + grain_bocager_500m + ai_l+moy.intensity.tot, data = df2[which(df2$saison == "printemps"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

####Carabid richness spring####
#variables locales
mod <- lm(C.rs ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "printemps"),])#transformer en log
vif(mod)
stepAIC(mod)
mod <- lm(C.rs ~ ndvi, data = df2[which(df2$saison == 
                                                             "printemps"),])
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(C.rs ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(log(C.rs+1) ~ ndvi + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(log(C.rs+1) ~ ndvi + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "printemps"),]) 
fullModel3b <- lm(log(C.rs+1) ~ (ndvi + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = log(C.rs+1) ~ ndvi, data = df2[which(df2$saison == "printemps"),])

B <- lm(formula = log(C.rs+1) ~ 1 , data = df2[which(df2$saison == "printemps"),])

C <- lm(formula = log(C.rs+1) ~ 1, data = df2[which(df2$saison == "printemps"),])

AB <- lm(formula = log(C.rs+1) ~ ndvi, data = df2[which(df2$saison == "printemps"),])

AC <- lm(formula = log(C.rs+1) ~ ndvi, data = df2[which(df2$saison == "printemps"),])

BC <- lm(formula = log(C.rs+1) ~ 1, data = df2[which(df2$saison == "printemps"),])

ABC <- lm(formula = log(C.rs+1) ~ ndvi, data = df2[which(df2$saison == "printemps"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

mod <- lm(formula = log(C.rs + 1) ~ ndvi + ndvi_500 + grain_bocager_500m + 
            enn_c_1 + moy.intensity.tot + ndvi_500:moy.intensity.tot + 
            grain_bocager_500m:moy.intensity.tot + enn_c_1:moy.intensity.tot, 
          data = df2[which(df2$saison == "printemps"), ])

interplot(modF, var1 = "moy.intensity.tot", var2 = "enn_c_1") + 
  ylab("") +
  xlab("")+
  geom_hline(yintercept = 0, linetype = "dashed") + geom_line(size=1) +
  theme_classic()

pred_values <- ggemmeans(model = modF,
                         terms = c("moy.intensity.tot","enn_c_1[0,300,600]"),
                         ci.lvl = 0.95,
                         type = "fixed", #population-level predictions, with intervals considering the uncertainty in the variance parameters
                         typical = "mean",
                         back.transform = T)
plot(pred_values, add.data = T)+mon_theme

####Carabid pielou spring####
#variables locales
mod <- lm(pielou.car ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "printemps"),])#transformer en log
vif(mod)
stepAIC(mod)
mod <- lm(C.rs ~ ndvi, data = df2[which(df2$saison == 
                                          "printemps"),])
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(pielou.car ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "printemps"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(pielou.car ~ ndvi + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "printemps"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(pielou.car ~ ndvi + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "printemps"),]) 
fullModel3b <- lm(log(C.rs+1) ~ (ndvi + ai_l + enn_c_1) * moy.intensity.tot, 
                  data = df2[which(df2$saison == "printemps"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                        scope = list(upper = fullModel3,
                                     lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = pielou.car ~ ndvi, data = df2[which(df2$saison == "printemps"),])

B <- lm(formula = pielou.car ~ 1 , data = df2[which(df2$saison == "printemps"),])

C <- lm(formula = pielou.car ~ ai_l, data = df2[which(df2$saison == "printemps"),])

AB <- lm(formula = pielou.car ~ ndvi, data = df2[which(df2$saison == "printemps"),])

AC <- lm(formula = pielou.car ~ ndvi+ai_l, data = df2[which(df2$saison == "printemps"),])

BC <- lm(formula = pielou.car ~ ai_l, data = df2[which(df2$saison == "printemps"),])

ABC <- lm(formula = pielou.car ~ ndvi+ai_l, data = df2[which(df2$saison == "printemps"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

mod <- lm(formula = log(C.rs + 1) ~ ndvi + ndvi_500 + grain_bocager_500m + 
            enn_c_1 + moy.intensity.tot + ndvi_500:moy.intensity.tot + 
            grain_bocager_500m:moy.intensity.tot + enn_c_1:moy.intensity.tot, 
          data = df2[which(df2$saison == "printemps"), ])

interplot(modF, var1 = "moy.intensity.tot", var2 = "enn_c_1") + 
  ylab("") +
  xlab("")+
  geom_hline(yintercept = 0, linetype = "dashed") + geom_line(size=1) +
  theme_classic()

pred_values <- ggemmeans(model = modF,
                         terms = c("moy.intensity.tot","enn_c_1[0,300,600]"),
                         ci.lvl = 0.95,
                         type = "fixed", #population-level predictions, with intervals considering the uncertainty in the variance parameters
                         typical = "mean",
                         back.transform = T)
plot(pred_values, add.data = T)+mon_theme

####Carabid composition spring####
##acp
car <- df[which(df2$saison == "printemps"),c(110:142)]

somme <- vector(length = NCOL(car))
for (i in c(1:NCOL(car))){
  somme[i] <- sum(car[,i])  
}
car <- car[,which(somme>0)]

res.pca <- PCA(log(car+1), scale.unit = T, graph = T)
fviz_pca_var(res.pca, repel = T)+mon_theme
fviz_eig(res.pca)+mon_theme
nCng(res.pca$eig[,1])
dimdesc(res.pca)

pos.car <- res.pca$ind$coord[,1:3]
colnames(pos.car) = c("car.1","car.2","car.3")

df2.car <- data.frame(df2[which(df2$saison == "printemps"),],pos.car)

##Dim1
#variables locales
mod <- lm(car.1 ~ intensity + ndvi + microclimat, 
          data = df2.car)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(car.1 ~ intensity + microclimat + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.car)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(car.1 ~ intensity + microclimat + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2.car) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(car.1 ~ intensity + microclimat + Taille_parcelle_500 + ndvi_500 + 
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.car) 
fullModel3b <- lm(car.1 ~ (intensity + microclimat + Taille_parcelle_500 + ndvi_500 + 
                   ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.car) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = car.1 ~ intensity + microclimat, data = df2.car)

B <- lm(formula = car.1 ~ Taille_parcelle_500 + ndvi_500 , data = df2.car)

C <- lm(formula = car.1 ~ 1, data = df2.car)

AB <- lm(formula = car.1 ~ intensity + microclimat + Taille_parcelle_500 + ndvi_500, data = df2.car)

AC <- lm(formula = car.1 ~ intensity + microclimat, data = df2.car)

BC <- lm(formula = car.1 ~ Taille_parcelle_500 + ndvi_500, data = df2.car)

ABC <- lm(formula = car.1 ~ intensity + microclimat + Taille_parcelle_500 + ndvi_500, data = df2.car)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

##Dim2
#variables locales
mod <- lm(car.2 ~ intensity + ndvi + microclimat, 
          data = df2.car)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(car.2 ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.car)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(car.2 ~ ndvi + ndvi_500 + grain_bocager_500m, 
                 data = df2.car) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(car.2 ~ ndvi + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.car) 
fullModel3b <- lm(car.2 ~ (ndvi + ai_l + enn_c_1)* moy.intensity.tot, 
                 data = df2.car)
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = car.2 ~ ndvi, data = df2.car)

B <- lm(formula = car.2 ~ 1 , data = df2.car)

C <- lm(formula = car.2 ~ 1, data = df2.car)

AB <- lm(formula = car.2 ~ ndvi, data = df2.car)

AC <- lm(formula = car.2 ~ ndvi, data = df2.car)

BC <- lm(formula = car.2 ~ 1, data = df2.car)

ABC <- lm(formula = car.2 ~ ndvi, data = df2.car)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C =  r.squaredGLMM(C)[1,1], 
        AB =  r.squaredGLMM(AB)[1,1], AC =  r.squaredGLMM(AC)[1,1], 
        BC =  r.squaredGLMM(BC)[1,1], ABC =  r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Practices", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


##Dim3
#variables locales
mod <- lm(car.3 ~ intensity + ndvi + microclimat, 
          data = df2.car)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(car.3 ~ microclimat + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.car)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(car.3 ~ microclimat + X._prairies_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2.car) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(car.3 ~ microclimat + X._prairies_500 + Taille_parcelle_500 + 
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.car) 
fullModel3b <- lm(car.3 ~ (microclimat + X._prairies_500 + Taille_parcelle_500 + 
                   ai_l + enn_c_1)* moy.intensity.tot, 
                 data = df2.car)
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = car.3 ~ microclimat, data = df2.car)

B <- lm(formula = car.3 ~ X._prairies_500 + Taille_parcelle_500 , data = df2.car)

C <- lm(formula = car.3 ~ 1, data = df2.car)

AB <- lm(formula = car.3 ~ X._prairies_500 + Taille_parcelle_500 + microclimat, data = df2.car)

AC <- lm(formula = car.3 ~ microclimat, data = df2.car)

BC <- lm(formula = car.3 ~ X._prairies_500 + Taille_parcelle_500, data = df2.car)

ABC <- lm(formula = car.3 ~ X._prairies_500 + Taille_parcelle_500 + microclimat, data = df2.car)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

####Spiders abundance autumn####
#variables locales
mod <- lm(log(A.ab+1) ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
stepAIC(mod)
mod <- lm(log(A.ab+1) ~ microclimat, data = df2[which(df2$saison == 
                                                 "automne"),])
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(log(A.ab+1) ~ microclimat + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(log(A.ab+1) ~ microclimat + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(log(A.ab+1) ~ microclimat + ndvi_500 + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
fullModel3b <- lm(log(A.ab+1) ~ (microclimat + ndvi_500 + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))

r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = log(A.ab+1) ~ microclimat, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = log(A.ab+1) ~ ndvi_500 , data = df2[which(df2$saison == "automne"),])

C <- lm(formula = log(A.ab+1) ~ moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = log(A.ab+1) ~ microclimat + ndvi_500 , data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = log(A.ab+1) ~ microclimat +moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = log(A.ab+1) ~ ndvi_500 + moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = log(A.ab+1) ~ microclimat +ndvi_500 + moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

####Spiders richness autumn####
#variables locales
mod <- lm(A.rs ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(A.rs ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(A.rs ~ ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(A.rs ~ ndvi_500 + grain_bocager_500m + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
fullModel3b <- lm(A.rs ~ (ndvi_500 + grain_bocager_500m + ai_l + enn_c_1 ) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = A.rs ~ 1, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = A.rs ~ ndvi_500 + grain_bocager_500m , data = df2[which(df2$saison == "automne"),])

C <- lm(formula = A.rs ~moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = A.rs ~ ndvi_500 + grain_bocager_500m, data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = A.rs ~ moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = A.rs ~ ndvi_500 + grain_bocager_500m+moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = A.rs ~ndvi_500 + grain_bocager_500m+moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)



####Spiders pielou index autumn####
#variables locales
mod <- lm(pielou.ara ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(pielou.ara ~ ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(pielou.ara ~ ndvi + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(pielou.ara ~ ndvi + ndvi_500 + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
fullModel3b <- lm(A.rs ~ (ndvi_500 + grain_bocager_500m + ai_l + enn_c_1 ) * moy.intensity.tot, 
                  data = df2[which(df2$saison == "automne"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                        scope = list(upper = fullModel3,
                                     lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = pielou.ara ~ ndvi, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = pielou.ara ~ ndvi_500, data = df2[which(df2$saison == "automne"),])

C <- lm(formula = pielou.ara ~ enn_c_1 + ai_l, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = pielou.ara ~ ndvi + ndvi_500, data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = pielou.ara ~ ndvi + enn_c_1 + ai_l, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = pielou.ara ~ ndvi_500 + enn_c_1 + ai_l, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = pielou.ara ~enn_c_1+ ai_l+ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

####Spiders composition autumn####
##acp
ara <- df[which(df2$saison == "automne"),c(73:106)]

somme <- vector(length = NCOL(ara))
for (i in c(1:NCOL(ara))){
  somme[i] <- sum(ara[,i])  
}
ara <- ara[,which(somme>0)]

res.pca <- PCA(log(ara+1), scale.unit = T, graph = T)
fviz_pca_var(res.pca, repel = T)+mon_theme
fviz_eig(res.pca)
nCng(res.pca$eig[,1])
dimdesc(res.pca)

pos.ara <- res.pca$ind$coord[,1:3]
colnames(pos.ara) = c("ara.1","ara.2","ara.3")

df2.ara <- data.frame(df2[which(df2$saison == "automne"),],pos.ara)

##Dim1
#variables locales
mod <- lm(ara.1 ~ intensity + ndvi + microclimat, 
          data = df2.ara)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(ara.1 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.ara)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(ara.1 ~ ndvi_500 + grain_bocager_500m, 
                 data = df2.ara) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(ara.1 ~ ndvi_500 + grain_bocager_500m + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.ara) 
fullModel3b <- lm(ara.1 ~ (ndvi_500 + grain_bocager_500m + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.ara) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = ara.1 ~ 1, data = df2.ara)

B <- lm(formula = ara.1 ~ndvi_500 + grain_bocager_500m, data = df2.ara)

C <- lm(formula = ara.1 ~ moy.intensity.tot, data = df2.ara)

AB <- lm(formula = ara.1 ~ ndvi_500 + grain_bocager_500m, data = df2.ara)

AC <- lm(formula = ara.1 ~ moy.intensity.tot, data = df2.ara)

BC <- lm(formula = ara.1 ~ndvi_500 + grain_bocager_500m+moy.intensity.tot, data = df2.ara)

ABC <- lm(formula = ara.1 ~ ndvi_500 + grain_bocager_500m+moy.intensity.tot, data = df2.ara)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

##Dim2
#variables locales
mod <- lm(ara.2 ~ intensity + ndvi + microclimat, 
          data = df2.ara)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(ara.2 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.ara)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(ara.2 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2.ara) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(ara.2 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500 + 
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.ara) 
fullModel3b <- lm(ara.2 ~ (X._prairies_500 + diversite_500 + Taille_parcelle_500 + 
                   ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.ara) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(stepAIC(modF <- fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = ara.2 ~ 1, data = df2.ara)

B <- lm(formula = ara.2 ~X._prairies_500 + diversite_500 + Taille_parcelle_500, data = df2.ara)

C <- lm(formula = ara.2 ~ ai_l , data = df2.ara)

AB <- lm(formula = ara.2 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, data = df2.ara)

AC <- lm(formula = ara.2 ~ ai_l , data = df2.ara)

BC <- lm(formula = ara.2 ~X._prairies_500 + diversite_500 + Taille_parcelle_500+ai_l , data = df2.ara)

ABC <- lm(formula = ara.2 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500+ai_l , data = df2.ara)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)

##Dim3
#variables locales
mod <- lm(ara.3 ~ intensity + ndvi + microclimat, 
          data = df2.ara)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(ara.3 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.ara)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(ara.3 ~ ndvi + ndvi_500 + grain_bocager_500m, 
                 data = df2.ara) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(ara.3 ~ ndvi_500 + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.ara) 
fullModel3b <- lm(ara.3 ~ (ndvi_500 + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.ara)
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

interplot(modF, var1 = "moy.intensity.tot", var2 = "ndvi_500") + 
  ylab("") +
  xlab("")+
  geom_hline(yintercept = 0, linetype = "dashed") + geom_line(size=1) +
  theme_classic()


A <- lm(formula = ara.3 ~ 1, data = df2.ara)

B <- lm(formula = ara.3 ~ndvi_500, data = df2.ara)

C <- lm(formula = ara.3 ~ enn_c_1 , data = df2.ara)

AB <- lm(formula = ara.3 ~ ndvi_500, data = df2.ara)

AC <- lm(formula = ara.3 ~ enn_c_1 , data = df2.ara)

BC <- lm(formula = ara.3 ~ndvi_500+enn_c_1 , data = df2.ara)

ABC <- lm(formula = ara.3 ~ ndvi_500+enn_c_1 , data = df2.ara)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


####Carabids abundance autumn####
#variables locales
mod <- lm(log(C.ab+1) ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
stepAIC(mod)

summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(log(C.ab+1) ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(log(C.ab+1) ~ diversite_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(log(C.ab+1) ~ diversite_500 + Taille_parcelle_500 + ndvi_500 + 
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
fullModel3b <- lm(log(C.ab+1) ~ (diversite_500 + Taille_parcelle_500 + ndvi_500 + 
                   ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = log(C.ab+1) ~ 1, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = log(C.ab+1) ~ diversite_500 + Taille_parcelle_500 + ndvi_500 , data = df2[which(df2$saison == "automne"),])

C <- lm(formula = log(C.ab+1) ~ moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = log(C.ab+1) ~ diversite_500 + Taille_parcelle_500 + ndvi_500 , data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = log(C.ab+1) ~ moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = log(C.ab+1) ~ diversite_500 + Taille_parcelle_500 + ndvi_500 + moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = log(C.ab+1) ~ diversite_500 + Taille_parcelle_500 + ndvi_500 + moy.intensity.tot + enn_c_1, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


####Carabids richness autumn####
#variables locales
mod <- lm(C.rs ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
stepAIC(mod)
mod <- lm(C.rs ~ microclimat, data = df2[which(df2$saison == 
                                          "automne"),])
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(C.rs ~ microclimat + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(log(C.rs+1) ~ microclimat + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(log(C.rs+1) ~  microclimat  + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
fullModel3b <- lm(log(C.rs+1) ~  (microclimat  + ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = log(C.rs+1) ~ microclimat, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = log(C.rs+1) ~ 1 , data = df2[which(df2$saison == "automne"),])

C <- lm(formula = log(C.rs+1) ~1, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = log(C.rs+1) ~ microclimat, data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = log(C.rs+1) ~ microclimat, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = log(C.rs+1) ~ 1, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = log(C.rs+1) ~ microclimat, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


####Carabids pielou index autumn####
#variables locales
mod <- lm(pielou.car ~ intensity + ndvi + microclimat, 
          data = df2[which(df2$saison == "automne"),])#transformer en log
vif(mod)
stepAIC(mod)
mod <- lm(pielou.car ~ intensity + ndvi, data = df2[which(df2$saison == 
                                                 "automne"),])
summary(mod)

#variables d'occupation du sol --> structure
mod.occ.sol <- lm(pielou.car ~ intensity + ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2[which(df2$saison == "automne"),])#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(pielou.car ~ intensity + ndvi + X._prairies_500 + 
                   diversite_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2[which(df2$saison == "automne"),]) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(pielou.car ~ intensity + ndvi + X._prairies_500 + 
                   diversite_500 + Taille_parcelle_500 + ndvi_500 + ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2[which(df2$saison == "automne"),]) 
fullModel3b <- lm(log(C.rs+1) ~  (microclimat  + ai_l + enn_c_1) * moy.intensity.tot, 
                  data = df2[which(df2$saison == "automne"),]) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                        scope = list(upper = fullModel3,
                                     lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = pielou.car ~ intensity + ndvi, data = df2[which(df2$saison == "automne"),])

B <- lm(formula = pielou.car ~ X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500, data = df2[which(df2$saison == "automne"),])

C <- lm(formula = pielou.car ~ enn_c_1, data = df2[which(df2$saison == "automne"),])

AB <- lm(formula = pielou.car ~ intensity + ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500, data = df2[which(df2$saison == "automne"),])

AC <- lm(formula = pielou.car ~ intensity + ndvi + enn_c_1, data = df2[which(df2$saison == "automne"),])

BC <- lm(formula = pielou.car ~ X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500 + enn_c_1, data = df2[which(df2$saison == "automne"),])

ABC <- lm(formula = pielou.car ~enn_c_1 + intensity + ndvi + X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500, data = df2[which(df2$saison == "automne"),])


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)



####Carabids composition autumn####
##acp
car <- df[which(df2$saison == "automne"),c(110:142)]

somme <- vector(length = NCOL(car))
for (i in c(1:NCOL(car))){
  somme[i] <- sum(car[,i])  
}
car <- car[,which(somme>0)]

res.pca <- PCA(log(car+1), scale.unit = T, graph = T)
fviz_pca_var(res.pca, repel = T)+mon_theme
fviz_eig(res.pca)
nCng(res.pca$eig[,1])
dimdesc(res.pca)

pos.car <- res.pca$ind$coord[,1:3]
colnames(pos.car) = c("car.1","car.2","car.3")

df2.car <- data.frame(df2[which(df2$saison == "automne"),],pos.car)

##Dim1
#variables locales
mod <- lm(car.1 ~ intensity + ndvi + microclimat, 
          data = df2.car)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(car.1 ~ intensity + X._prairies_500 + Taille_parcelle_500, 
                  data = df2.car)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(car.1 ~ intensity + X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2.car) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(car.1 ~ intensity + X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500 + 
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.car) 
fullModel3b <- lm(car.1 ~ (intensity + X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500 + 
                   ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.car) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = car.1 ~ intensity, data = df2.car)

B <- lm(formula = car.1 ~X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500, data = df2.car)

C <- lm(formula = car.1 ~ moy.intensity.tot, data = df2.car)

AB <- lm(formula = car.1 ~ intensity+X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500, data = df2.car)

AC <- lm(formula = car.1 ~ intensity+moy.intensity.tot, data = df2.car)

BC <- lm(formula = car.1 ~X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500+moy.intensity.tot, data = df2.car)

ABC <- lm(formula = car.1 ~ intensity+X._prairies_500 + diversite_500 + Taille_parcelle_500 + ndvi_500+moy.intensity.tot, data = df2.car)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)


##Dim2
#variables locales
mod <- lm(car.2 ~ intensity + ndvi + microclimat, 
          data = df2.car)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(car.2 ~ intensity + X._prairies_500 + Taille_parcelle_500, 
                  data = df2.car)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(car.2 ~ intensity + X._prairies_500 + ndvi_500 + grain_bocager_500m, 
                 data = df2.car) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(car.2 ~ intensity + X._prairies_500 + diversite_500 + grain_bocager_500m + 
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.car) 
fullModel3b <- lm(car.2 ~ (intensity + X._prairies_500 + grain_bocager_500m + 
                   ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.car) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

interplot(modF, var1 = "moy.intensity.tot", var2 = "intensity") + 
  ylab("") +
  xlab("")+
  geom_hline(yintercept = 0, linetype = "dashed") + geom_line(size=1) +
  theme_classic()

A <- lm(formula = car.2 ~ intensity, data = df2.car)

B <- lm(formula = car.2 ~X._prairies_500 + diversite_500 + grain_bocager_500m, data = df2.car)

C <- lm(formula = car.2 ~ moy.intensity.tot+ai_l, data = df2.car)

AB <- lm(formula = car.2 ~ intensity+X._prairies_500 + diversite_500 + grain_bocager_500m, data = df2.car)

AC <- lm(formula = car.2 ~ intensity+moy.intensity.tot+ai_l, data = df2.car)

BC <- lm(formula = car.2 ~X._prairies_500 + diversite_500 + grain_bocager_500m+moy.intensity.tot+ai_l, data = df2.car)

ABC <- lm(formula = car.2 ~ intensity+X._prairies_500 + diversite_500 + grain_bocager_500m+moy.intensity.tot+ai_l, data = df2.car)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)
##Dim3
#variables locales
mod <- lm(car.3 ~ intensity + ndvi + microclimat, 
          data = df2.car)#transformer en log
vif(mod)
summary(stepAIC(mod))


#variables d'occupation du sol --> structure
mod.occ.sol <- lm(car.3 ~ X._prairies_500 + diversite_500 + Taille_parcelle_500, 
                  data = df2.car)#structure
vif(mod.occ.sol)

summary(stepAIC(mod.occ.sol, 
                scope = list(upper = mod.occ.sol,
                             lower = stepAIC(mod))))


#3
fullModel2 <- lm(car.3 ~ ndvi_500 + grain_bocager_500m, 
                 data = df2.car) 
minModel2 <-  stepAIC(mod.occ.sol, 
                      scope = list(upper = mod.occ.sol,
                                   lower = stepAIC(mod)))

summary(stepAIC(fullModel2,
                scope = list(upper = fullModel2,
                             lower = minModel2)))


#4 <- pratiques
fullModel3 <- lm(car.3 ~ ndvi_500 + 
                   ai_l + moy.intensity.tot + enn_c_1, 
                 data = df2.car) 
fullModel3b <- lm(car.3 ~ (ndvi_500 + 
                   ai_l + enn_c_1) * moy.intensity.tot, 
                 data = df2.car) 
minModel3 <-  stepAIC(fullModel2,
                      scope = list(upper = fullModel2,
                                   lower = minModel2))

summary(modF <- stepAIC(fullModel3, direction = "backward",
                scope = list(upper = fullModel3,
                             lower = minModel3)))
r2beta(modF, method = "lm")
omega_squared(modF)

A <- lm(formula = car.3 ~ 1, data = df2.car)

B <- lm(formula = car.3 ~ndvi_500, data = df2.car)

C <- lm(formula = car.3 ~ moy.intensity.tot+ai_l, data = df2.car)

AB <- lm(formula = car.3 ~ ndvi_500, data = df2.car)

AC <- lm(formula = car.3 ~ moy.intensity.tot+ai_l, data = df2.car)

BC <- lm(formula = car.3 ~ndvi_500+moy.intensity.tot+ai_l, data = df2.car)

ABC <- lm(formula = car.3 ~ ndvi_500+moy.intensity.tot+ai_l, data = df2.car)


varPart(A = r.squaredGLMM(A)[1,1], B =  r.squaredGLMM(B)[1,1], C = r.squaredGLMM(C)[1,1],
        AB =  r.squaredGLMM(AB)[1,1], BC = r.squaredGLMM(BC)[1,1],
        AC = r.squaredGLMM(AC)[1,1], ABC = r.squaredGLMM(ABC)[1,1],
        A.name = "Local", B.name = "Landscape", C.name = "Pratiques", plot = TRUE,
        plot.digits = 3, cex.names = 1.5, cex.values = 1.2, main = "", cex.main = 2,
        plot.unexpl = TRUE)