\documentclass[a4paper]{article}
\usepackage{fullpage} 
\usepackage{parskip} 
\usepackage{amsmath}
\usepackage{anysize}
\usepackage[dvipsnames]{xcolor}
\usepackage{hyperref}
\usepackage[spanish]{babel}  
\usepackage[utf8x]{inputenc}   
\usepackage[margin=2cm]{geometry}
\usepackage{multirow}
\usepackage{booktabs}
\usepackage{mathtools}
\usepackage{float}
\usepackage{adjustbox}
\usepackage[bottom]{footmisc}
\usepackage{makecell}
    \setcellgapes{5pt}
\usepackage[labelfont=bf]{caption}
\newcommand*{\torange}{\textcolor{BurntOrange}}
\newcommand*{\tblue}{\textcolor{BlueGreen}}


\begin{document}
\setcounter{page}{1}

\thispagestyle{empty}

\begin{center}
{\LARGE{\bfseries Trabajo Práctico N\textsuperscript{o}3 - Modelos de Clasificación}}\\
\hspace
{\large{\bfseries Aprendizaje Estadístico - FIUBA}}\\
{\large{\bfseries 1\textsuperscript{er} Cuatrimestre - 2020}}
\end{center}

\hspace

\begin{center}
{\large{\textfont{José F. González - 100063 - \url{jfgonzalez@fi.uba.ar}}}\\
\end{center}

\begin{center}
{\large{\textfont{Atento a: Ing. Jemina García}}} \hfill {\large{\textfont{Revisión 1.8}}}\\
\end{center}


\section{Introducción}

<<echo=FALSE>>=
## Cargando datos
levels <- read.csv("Colon_X.txt", sep=" ", header=FALSE)
## Elimina una linea vacia al final
levels <- levels[1:62]
## Datos con logaritmo
logLevels <- log10(levels)
@

Nos interesa contruir un modelo de clasificación para detectar la presencia de tumores en muestras de tejido celular. Para ello disponemos el conjunto de datos \torange{\texttt{levels}} con mediciones de los niveles de actividad de 2000 genes presentes en células de 62 tejidos distintos, donde 40 tejidos están etiquetados como tumores y 22 como normales.
\textbf{Como criterio general buscamos un modelo con pocas variables explicativas que maximize la medida \textit{accuracy} en una validación cruzada.} 

\section{Estandarización}

Los niveles de actividad tienen un rango muy amplio lo que se traduce en una alta variabilidad, esto es evidente el panel izquierdo de la Figura \ref{fig:variabilidad} donde se muestran la media y el desvío estándar para los niveles de cada una de las muestras de tejido. Aplicando un logaritmo base diez a todos los datos se comprime el rango de datos y se obtiene un desvío estándar más representativo del conjunto como se ve en el panel central de la Figura \ref{fig:variabilidad}, definimos así un nuevo conjunto \torange{\texttt{logLevels}}.

Los métodos de reducción de variables que utilizaremos, para enfrentar la desproporción entre cantidad de predictoras y muestras, funcionan mejor con predictores estandarizados en media cero y desvío uno, por ejemplo para no afectar las penalizaciones de los métodos de reducción. Luego resultará conveniente definir un nuevo conjunto \torange{\texttt{stdLogLevels}} en esta escala a cambio de perder un poco de interpretabilidad. En el panel derecho de la Figura \ref{fig:variabilidad} se ve el efecto de la nueva escala sobre la media y desvío de los genes de cada uno de los tejidos.

\begin{figure}[h]
\caption{Variabilidad de la media y el desvío estandar en las 62 muestras de tejido. \textbf{Izquierda:} Desvío estandar y media para las 62 muestras de tejido del conjunto \torange{\texttt{levels}}. \textbf{Centro:} Desvío estándar y media para las 62 muestras del conjunto \torange{\texttt{logLevels}} luego de tomar logaritmo en base diez a los niveles de expresión de los genes. \textbf{Derecha:} Desvío estándar y media para las muestras de \torange{\texttt{stdLogLevels}} luego de estandarizar los datos a media cero y desvío uno.}

\label{fig:variabilidad}
<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cache=TRUE, fig.height=3>>=
par(mfrow=c(1,3), mar=c(4, 4, 1, 3))  

####
## Justificación del logaritmo
## Media y desvio estandar para cada muestra sin aplicar logaritmo
m <- seq(1,62,1)
s <- seq(1,62,1)
for (tissue in seq(1,62,1) ){
  m[tissue] <- mean(levels[,tissue])
  s[tissue] <- sd(levels[,tissue])
}

## Grafico Media y Desvio sin aplicar logaritmo
plot(m[1],s[1],
     col="cadetblue3",
     ylab=expression("Desvío Estándar"),
     type = "p",
     xlab=expression("Media"),
     cex.lab=0.75,
     cex.axis=0.75,xaxs="i", yaxs="i")
points(m,
       s,
       col="cadetblue3")
abline(lm(s~m), lty="dashed")


## Media y desvio estandar para cada muestra aplicando logaritmo
m <- seq(1,62,1)
s <- seq(1,62,1)
for (tissue in seq(1,62,1) ){
  m[tissue] <- mean(logLevels[,tissue])
  s[tissue] <- sd(logLevels[,tissue])
}
## Grafico sin aplicar logaritmo
plot(m[1],s[1],
     col="cadetblue3",
     ylab=expression("Desvío Estándar"),
     type = "p",
     xlab=expression("Media"),
     cex.lab=0.75,
     cex.axis=0.75,xaxs="i", yaxs="i")
## Puntos
points(m,
       s,
       col="cadetblue3")
## Regresion
abline(lm(s~m), lty="dashed")


## Estandarizamos a media cero y desvio 1
logLevelsMatrix <- as.matrix(sapply(logLevels, as.numeric)) 
stdLogLevels <- scale(logLevels)
m <- seq(1,62,1)
s <- seq(1,62,1)
for (tissue in seq(1,62,1) ){
  m[tissue] <- mean(stdLogLevels[,tissue])
  s[tissue] <- sd(stdLogLevels[,tissue])
}
plot(m,s,
     col="cadetblue3",
     ylab=expression("Desvío Estándar"),
     type = "p",
     xlab=expression("Media"),
     cex.lab=0.75,
     cex.axis=0.75,xaxs="i", yaxs="i")
## Regresion
abline(lm(s~m), lty="dashed")
@
\end{figure}


\section{Reducción de variables}

El mapa de calor de de la Figura \ref{fig:heatmap} muestra los logaritmos de niveles de actividad, en los distintos genes del eje vertical, como tonos rojos para genes con mayor expresión y tonos azules para genes con menor expresión. Los genes fueron agrupados en un cluster jerárquico según los patrones de su cadena de intesidad revelando estructuras similares en las expresiones de distintos genes. Esto se interpreta como procesos de regulación entre ciertos grupos de genes, lo que sugiere efectos de multicolinealidad fuertes. Analizando la matriz de correlación de \torange{\texttt{stdLogLevels}} se cuenta un $9\%$ de las correlaciones entre genes son mayor a $0.7$ y un $2\%$ son correlación mayor a $0.8$. Los casos de alta correlación son generalmente interpretados como genes con regulación directa entre ellos. 


<<echo=FALSE>>=
## Transpone los datos
logLevelsMatrix <- as.matrix(sapply(logLevels, as.numeric)) 
stdLogLevels <- scale(t(logLevelsMatrix))
corData <- cor(stdLogLevels)

## Porcentaje de genes con correlacion mayor al 70%
count = 0
for (x in seq(1,2000,1)) {
  for (y in seq(1,2000,1)) {
    if(corData[x,y] > 0.8) {
      count = count + 1
    }
  }
}
# Porcentaje
#(count)/length(corData)*100
@

\begin{figure}[H]
<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cache=TRUE, fig.height=3>>=
## Mapa de calor con cluster
library(pheatmap)
b <- pheatmap(logLevelsMatrix[1:100,],square=TRUE, treeheight_row = 0, border_color = NA, treeheight_col = 0,cutree_rows=6,  cluster_rows = T,cluster_cols = F,  labels_col = "", legend=FALSE)
@

\caption{Actividad de los distintos genes. El eje vertical corresponde a genes y el horizontal a los tejidos. Genes con mayor nivel de actividad están codificados con tonos rojos y los de menor actividad con azules. Los genes se encuentran ordenados por un \textit{clustering jerarquíco} según similitudes en la cadena de intesidades para facilitar la distinción de grupos correlacionados.}
\label{fig:heatmap}
\end{figure}


\subsection{Selección por \textit{Lasso}}

El primer enfoque para reducir dimensión y seleccionar variables simultaneamente es penalizar por \textit{LASSO}. En la Figura \ref{fig:lasso} se muestra el $\lambda$ que minimiza el error cuadrático medio por validación cruzada para Lasso con \torange{\texttt{nfolds=10}}, para el valor obtenido sobreviven 9 genes de los 2000 iniciales. 


\begin{figure}[H]

<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cache=TRUE, fig.height=3>>=
## Lasso
library(plotmo)
library(glmnet)

## para que sea reproducible
set.seed(1)

layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))

## Etiquetas para clasificar
labels <- read.csv("colonT.txt", header=FALSE)
labels <- as.matrix(sapply(labels, as.numeric)) 

## Ajuste de lasso
ajus_lasso<-glmnet(stdLogLevels,labels,alpha=1, family = "binomial")

## Tuneado de lambda por cv con misclassification
set.seed(1)
las.crossval = cv.glmnet(stdLogLevels, labels, family = "binomial", type.measure = "class", nfolds=10, alpha=1)

## Colores para el plot del lasso
cols <- c("#5F9EA0", "#339FA3", "#D69545", "#B36505",  "#7D7956", "#FF8C00", "#EDE7B1")
palette(cols);

## Plot del lasso
plot_glmnet(ajus_lasso, label=9, col=cols, ylab=expression("Coeficientes"))
abline(lty="dashed", v=log(las.crossval$lambda.min))

# Plot de Misclasification Logistica segun lambda por CV
plot(las.crossval$cvm, x=log(las.crossval$lambda), type="l", col = "cadetblue3",
xlab=expression("Log Lambda"),
ylab=expression("Error de clasificación por CV"), main="Lasso")
abline(lty="dashed", v=log(las.crossval$lambda.min))

# Las mismas etiquetas pero con nombres "Tumor" y "Normal"
labelsNames <- read.csv("colonTnames.txt", header=FALSE)
labelsNames <- as.matrix(sapply(labelsNames, as.factor)) 
@

\end{figure}

\clearpage

\subsection{Selección por Elastic Net}

\textbf{Bajo efectos de colinealidad \textit{Lasso} elige arbitrariamente entre muchos posibles conjuntos colineales de predictoras omitiendo algunos genes que puedan estar relacionados con la salida}. \textit{Elastic Net} combina las penalizaciones de \textit{Ridge} y \textit{Lasso} obteniendo un método que es más inclusivo con las variables correlacionadas, por efector de la penalización \textit{Ridge}, y da mejores resultado en este tipo de datos a costa de aumentar la cantidad de genes y tener que ajustar un nuevo parámetro $\alpha$. En la Figura XX se muestran los resultados de \rextit{Elastic Net} para tres valores de $\alpha$.



\begin{figure}[t]
<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cahce=TRUE, fig.height=3>>=
layout(matrix(c(1,1,2), nrow = 1, ncol = 3, byrow = TRUE))

## Enet 0.5
enet05<-glmnet(stdLogLevels, labels, alpha=0.5, family = "binomial")

## Misclassification de modelo logístico con enet alpha 0.5 - Tuneado de lambda
set.seed(1)
enet.05.cv = cv.glmnet(stdLogLevels, labels, family = "binomial", type.measure = "class", nfolds=10, alpha=0.5)

## Grafico de enet 0.5
plot_glmnet(enet05, label=15, col=cols, ylab=expression("Coeficientes"))
abline(lty="dashed", v=log(enet.05.cv$lambda.min))

# Grafico de cv
plot(enet.05.cv$cvm, x=log(enet.05.cv$lambda), type="l", col = "cadetblue3",
xlab=expression("Log Lambda"),
ylab=expression("Error de clasificación por CV"), main="Alpha=0.5")
abline(lty="dashed", v=log(enet.05.cv$lambda.min))

@

\caption{Comparación de \textit{Lasso}, \textit{Ridge} y \textit{Elastic Net} en función de los errores de clasificación por validación cruzada de un modelo logístico construido con los coefficientes de cada método de reducción}
\label{fig:lasso_comp}
\end{figure}


<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cahce=TRUE, fig.height=3>>=

## Variables elegidas por enet 05
enet.05<-glmnet(stdLogLevels, labels, lambda = enet.05.cv$lambda.min ,alpha=0.5)
enet.05$df
#prediccion <- predict(lasso.cv, x=stdLogLevels, y=labels, newx=stdLogLevels,
#                    s = "lambda.min",
#                    type = "class")
#predichos = rep("No", length(prediccion))
#predichos[prediccion > 0.74] <- "Sí"
#observados <- labelsNames
#confusion <- table(predichos, observados)
#confusion

#library(pROC)
#prediccion
#labels
#r<-roc(labels,as.numeric(prediccion))
#plot(x=1-r$specificities,y=r$sensitivities,
#	col="cadetblue3",
#	lwd=1,
#	type="l",
#	xlab="FPR",
#    ylab="TPR",
#    cex.lab=0.75,
#	cex.axis=0.75)
#points(1-r$specificities[81], r$sensitivities[81])
#text(1-r$specificities[81]-0.03, r$sensitivities[81]+0.05, labels=round(r$thresholds[81],2), cex= 0.7)

@

\subsection{Selección por \textit{Nearest Shrunken Centroids}}

\textit{Nearest Shrunken Centroids} es un método desarrollado en el contexto de identificar clases de genes y selección automática de genes que caracterizan la clase. Es una modificación del \texit{nearest centroids} aplicando una contracción de los centroides de cada clase hacía cero, eliminando los que alcanzan el cero y clasificando sobre la nueva disposición por proximidad a los centroides que sobreviven. En la Figura ?? se muestran los resultados del error de clasificación por validación cruzada según el parámetro de umbral $\Delta$. Como clasificador el método no es notablemente mejor que las regresiones logísticas con reducción por \textit{Lasso} y \texit{Elastic Net}, pero puede ser un buen selector de variables. 

\begin{figure}[b]
<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cache=TRUE, fig.height=3, results='hide'>>=

library(pamr)

set.seed(1)
## Pamr pide una matrix de la forma
pamr.data <- list(x=t(stdLogLevels),y=factor(labelsNames))

# Entrenar NSC
nsc<-pamr.train(pamr.data)

# CV para el factor de shrinkage
nsc.cv<-pamr.cv(nsc, pamr.data)
nsc.cv$size
## Modifico la funcion pamr.plotcv
body(pamr.plotcv)[[3]] <- substitute(par(mfrow=c(1,2), mar=c(4, 4, 3, 1)))
body(pamr.plotcv)[[21]] <- substitute(print("ok"))
body(pamr.plotcv)[[2]] <- substitute(print("ok"))
body(pamr.plotcv)[[29]] <- substitute(print("ok"))
body(pamr.plotcv)[[28]] <- substitute(legend(0, 0.9, dimnames(table(y))[[1]], col = (2:(nc + 1)), lty = 1, cex = 0.75))
body(pamr.plotcv)[[15]] <- substitute(axis(3, at = fit$threshold, labels = paste(fit$size), srt = 90, 
    adj = 0,  cex.axis=0.75))
body(pamr.plotcv)[[17]] <- substitute(axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8),  cex.axis=0.75))
body(pamr.plotcv)[[25]] <- substitute(axis(3, at = fit$threshold, labels = paste(fit$size), srt = 90, 
    adj = 0,cex.axis=0.75))
body(pamr.plotcv)[[26]] <- substitute(axis(2, at = c(0, 0.2, 0.4, 0.6, 0.8), cex.axis=0.75))
# Or: 
#error.bars(fit$threshold, fit$err - se, fit$err + se)
body(pamr.plotcv)[[16]] <- substitute(mtext("", 3, 4, cex = 0.75))
body(pamr.plotcv)[[14]] <- substitute(
plot(fit$threshold, fit$error, ylim = c(-0.1, 0.8), xlab = "Valor de Umbral", 
    ylab = "Error de clasificación", type = "n", yaxt = "n", cex = 0.75, cex.lab=0.75,cex.axis=0.75))
body(pamr.plotcv)[[24]] <- substitute(plot(fit$threshold, err2[1, ], ylim = c(-0.1, 1.1), xlab = "Valor de Umbral", cex = 0.75,cex.lab=0.75,
  cex.axis=0.75,
    ylab = "Error de Clasificación", type = "n", yaxt = "n"))
##
as.list(body(pamr.plotcv))
## Grafico
p<-pamr.plotcv(nsc.cv)
@
\caption{Comparación de \textit{Lasso}, \textit{Ridge} y \textit{Elastic Net} en función de los errores de clasificación por validación cruzada de un modelo logístico contruido con los coefficientes de cada método de reducción}
\label{fig:lasso_comp}
\end{figure}

<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cache=TRUE, fig.height=4, results='hide'>>=

## Listado de genes sobrevivientes a cierto umbral de minimo error
# Asi lo modifica la documentacion, para agregar un nombre al gen
pamr.data2 <- list(x=t(stdLogLevels),y=factor(labelsNames), geneid=as.character(1:nrow(t(stdLogLevels))),
               genenames=paste("g",as.character(1:nrow(t(stdLogLevels))),sep=""))

# Estos son los genes que elige NSC
survival.genes <- pamr.listgenes(nsc, pamr.data2, threshold=1.635)
# Grafico de centroides
#nsc.cv$threshold
#pamr.plotcen(nsc, pamr.data2,threshold=1.635)
#pamr.geneplot(nsc, pamr.data2, threshold=2.3)
 
## Voy a usar estos genes pero para clasificar por logistico
## Armo un dataframe con solo estos genes y sus etiquetas
survival.data <- stdLogLevels[,as.integer(survival.genes[,1])]
survival.data <- as.data.frame(survival.data)

survival.data$Labels <- labelsNames
survival.data$Labels <- factor(survival.data$Labels)

## Si uso todas las variables de nsc da muy mal, estan correlacionadas
model <- glm(Labels ~ ., data = survival.data, family = binomial)

## entonces le passo un elastic net 0.25
#set.seed(1)
#lasso.cv = cv.glmnet(stdLogLevels[,as.integer(survival.genes[,1])], labels, family = "binomial", type.measure = "class", nfolds=10, alpha=0.25)
#plot(lasso.cv$cvm, x=log(lasso.cv$lambda), type="l", col = "cadetblue3",
#xlab=expression("Log Lambda"),
#ylab=expression(""), main="Ridge")
#abline(lty="dashed", v=log(lasso.cv$lambda.min))

## Coeficientes de alpha 0.25 elastic net con lambdamin
#enet.025<-glmnet(stdLogLevels[,as.integer(survival.genes[,1])], labels, lambda = lasso.cv$lambda.min )

@

\clearpage
\subsection{Resultados}


\section{Clasificación}

Sobre el conjunto de variables seleccionadas evaluamos distintos métodos de clasificación. Para ello 
\subsection{SVM}

<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cache=TRUE, fig.height=4, results='hide'>>=

library(e1071)
library(pamr)
train.x <- t(pamr.data2$x) 
train.y <- pamr.data2$y 

svm.trained <- svm(train.x, train.y, kernel="linear")
summary(svm.trained)

svm.trained.cv <- svm(train.x, train.y, kernel="linear", cross=10)
summary(svm.trained.cv)
@


<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cache=TRUE, fig.height=3>>=
CrossVal <- function(data, labels, k=10)
{
CV.error.single <- c()
heaps <- balanced.folds(labels,k)
# CV is just a FOR-loop:
for (i in 1:k)
{
# 1. split data for training and testing
data.test <- data[, heaps[[i]]]
data.train <- data[,-heaps[[i]]]
labels.test <- labels[ heaps[[i]]]
labels.train <- labels[-heaps[[i]]]

# 2. select genes only on training data
tscores <- mt.teststat(data.train, labels.train, test="t")
selection <- order(tscores, decreasing=TRUE)[1:100]

# Seleccion por enet 0.5 para cada test
# set.seed(1) #Para que sea reproducible
# Elige el lambda por cv
enet.05.cv = cv.glmnet(t(data.train), as.integer(labels.train), family = "binomial", type.measure = "class", nfolds=10, alpha=0.5)
enet.05<-glmnet(t(data.train), as.integer(labels.train), lambda = enet.05.cv$lambda.min ,alpha=0.5)
selection <- (enet.05$beta)@i
train.sel <- data.train[selection, ]

# 3. train the model
svm.model <- svm(t(train.sel), as.factor(labels.train), kernel="linear")

# 4. test the model
test.sel <- t(data.test[selection, ])
predicted <- predict(svm.model, test.sel)
svm.error <- sum(predicted != labels.test)*100/nrow(test.sel)
CV.error.single <- c(CV.error.single, svm.error)
}
CV.error.total <- mean(CV.error.single)
accuracy <- 100-CV.error.total
cat(k,"-fold cross validation achieves ",accuracy,"% accuracy\n",sep="")
}
library(quantileDA)
library(multtest)
train.x
CrossVal(t(train.x), train.y)

@


\end{figure}


<<echo=FALSE, fig.align='center', message=FALSE, warning=FALSE, cache=TRUE, fig.height=3>>=
CrossVal <- function(data, labels, k=10)
{
CV.error.single <- c()
heaps <- balanced.folds(labels,k)
# CV is just a FOR-loop:
for (i in 1:k)
{
# 1. split data for training and testing
data.test <- data[, heaps[[i]]]
data.train <- data[,-heaps[[i]]]
labels.test <- labels[ heaps[[i]]]
labels.train <- labels[-heaps[[i]]]

# 2. select genes only on training data
tscores <- mt.teststat(data.train, labels.train, test="t")
selection <- order(tscores, decreasing=TRUE)[1:100]

# Seleccion por enet 0.5 para cada test
# set.seed(1) #Para que sea reproducible
# Elige el lambda por cv
enet.05.cv = cv.glmnet(t(data.train), as.integer(labels.train), family = "binomial", type.measure = "class", nfolds=10, alpha=0.5)
enet.05<-glmnet(t(data.train), as.integer(labels.train), lambda = enet.05.cv$lambda.min ,alpha=0.5)
selection <- (enet.05$beta)@i
train.sel <- data.train[selection, ]

# 3. train the model
svm.model <- svm(t(train.sel), as.factor(labels.train), kernel="linear")

# 4. test the model
test.sel <- t(data.test[selection, ])
predicted <- predict(svm.model, test.sel)
svm.error <- sum(predicted != labels.test)*100/nrow(test.sel)
CV.error.single <- c(CV.error.single, svm.error)
}
CV.error.total <- mean(CV.error.single)
accuracy <- 100-CV.error.total
cat(k,"-fold cross validation achieves ",accuracy,"% accuracy\n",sep="")
}
library(quantileDA)
library(multtest)
train.x
CrossVal(t(train.x), train.y)

@

\clearpage
\section{Bibliografía}
\begin{enumerate}
\bibitem{paper} Alon U, Barkai N, Notterman DA, Gish K, Ybarra S, Mack D, et al. Broad patterns of gene expression revealed by clustering analysis of tumor and normal colon tissues probed by oligonucleotide arrays. Proceedings of the National Academy of Sciences of the United States of America. 1999.
\bibitem{link1}https://stats.stackexchange.com/questions/184029/what-is-elastic-net-regularization-and-how-does-it-solve-the-drawbacks-of-ridge
\end{enumerate}
%https://stats.stackexchange.com/questions/74013/when-would-i-choose-lasso-over-elastic-net
%https://stats.stackexchange.com/questions/264016/why-lasso-or-elasticnet-perform-better-than-ridge-when-the-features-are-correlat
%https://stats.stackexchange.com/questions/345343/any-disadvantages-of-elastic-net-over-lasso

%https://stats.stackexchange.com/questions/298/in-linear-regression-when-is-it-appropriate-to-use-the-log-of-an-independent-va
%https://stats.stackexchange.com/questions/18844/when-and-why-should-you-take-the-log-of-a-distribution-of-numbers
%file:///home/jose/Downloads/LogTransform.pdf
%https://www.researchgate.net/post/What_is_Log_transformation_and_why_do_we_do_it_in_gene_expression_analysis    !! muy buenoo
% https://stats.stackexchange.com/questions/407563/choosing-model-for-more-predictors-than-observations

%https://www.reddit.com/r/bioinformatics/comments/17my8g/plotting_microarray_data_in_r/
%http://bioconductor.org/packages/release/bioc/html/affy.html
%http://rstudio-pubs-static.s3.amazonaws.com/14438_729bb8430be242fba106c8ae3968458f.html
%http://manuals.bioinformatics.ucr.edu/home/R_BioCondManual#TOC-Bioconductor
\section{Anexo}

\torange{\textit{Gene Expression}} es un proceso en que la información contenida en algún gen de una célula se dispara para generar un producto, por ejemplo una proteina. \torange{\textit{Gene Expression Level}} es una medida de esta actividad de los genes computada a partir de la concentracion de \textit{mRNA - messenger Ribonucleioc acid}, que cuando es medida puede no representar la verdadera actividad del gen pues el \textit{mRNA} puede estar siendo regulado por otros procesos.

\end{document}