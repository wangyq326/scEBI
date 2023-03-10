#installed.packages('spatstat')
#install.packages('readxl')
library(spatstat)
library(readxl)


#set direction
setwd("")
######################################################################
#########G245
G245data<-read_xlsx('all-G245S.xlsx')

G245data$X=G245data$X*0.31
G245data$Y=G245data$Y*0.31


G245ppp<-ppp(G245data$X, G245data$Y, xrange=range(G245data$X), 
             yrange=range(G245data$Y), marks = G245data$markers)

plot(G245ppp, pch=c(20,20,20), cols=c("blue", "red"), cex=0.5, legend=F)
legend("bottom",  c("DC", "EPC"), text.font=10,
       pch=c(20,20,20), col=c("blue", "red"), cex=1.0, box.col = "white")
dev.off()
####
G245data1<-G245data
G245data1$markers[G245data1$markers=="CD11c+TIM3+"]<-2
G245data1$markers[G245data1$markers=="GYPA+GAL9+"]<-1
G245data1$markers=factor(G245data1$markers, levels=c(1,2), labels = c("EPC", "DC") )
G245ppp1<-ppp(G245data1$X, G245data1$Y, xrange=range(G245data1$X), 
              yrange=range(G245data1$Y), marks = G245data1$markers)
set.seed(1234)
plot(envelope(G245ppp1, fun="Gcross",correction="KM", r=seq(0,100,20)), main="G245S_EPC_DC")
dev.off()