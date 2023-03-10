#COX
library("survival")
library("survminer")

cox=read.csv("2016 JOH COX.csv")

cox$OUTCOME2[cox$OUTCOME2 == "RFS"] =1
cox$OUTCOME2[cox$OUTCOME2 == "Dead"] = 2
cox$OUTCOME2=as.numeric(cox$OUTCOME2)

cox$Stage[cox$Stage == "1"] ="pretext1-2"
cox$Stage[cox$Stage == "2"] ="pretext1-2"
cox$Stage[cox$Stage == "3"] ="pretext3"

cox$cluster[cox$cluster=="1"]="HB1-3"
cox$cluster[cox$cluster=="2"]="HB2"
cox$cluster[cox$cluster=="3"]="HB1-3"


cox$risk[cox$risk=="intermidiate"]="0risk 2-3"
cox$risk[cox$risk=="low"]="risk 1"
cox$risk[cox$risk==""]=NA
cox$risk[cox$risk=="high"]="0risk 2-3"

cox$group[cox$group=="low"]="epc low"

cox$Histology[cox$Histology%in%c("fetal and embryonal","fetal")]="epithelial"
cox$Histology[cox$Histology%in%c("epithelial with SC","mix with sc","atypical")]="others"

covariates <- c("Age", "Gender", "group", "Stage","Histology","cluster")
univ_formulas <- sapply(covariates,
                        function(x) as.formula(paste('Surv(OS, OUTCOME2)~', x)))

univ_models <- lapply( univ_formulas, function(x){coxph(x, data = cox)})

univ_results <- lapply(univ_models,
                       function(x){ 
                         x <- summary(x)
                         p.value<-signif(x$wald["pvalue"], digits=2)
                         HR <-signif(x$coef[2], digits=2);
                         HR.confint.lower <- signif(x$conf.int[,"lower .95"], 2)
                         HR.confint.upper <- signif(x$conf.int[,"upper .95"],2)
                         HR <- paste0(HR, " (", 
                                      HR.confint.lower, "-", HR.confint.upper, ")")
                         res<-c(p.value,HR)
                         names(res)<-c("p.value","HR (95% CI for HR)")
                         return(res)
                       })
res<- as.data.frame(univ_results, check.names = FALSE)
res=t(res)
res=res[,1:2]

as.data.frame(res)

###forest plot
library(tableone)  
library(survival)
##install.packages("forestplot")
library(forestplot)
library(stringr)

result<- read.csv('os.csv',header=T)
head(result)

forestplot(result[,c(1,5,6)],             
           mean=result[,2],   
           lower=result[,3],           
           upper=result[,4],      
           zero=1,            
           boxsize=0.3,      
           graph.pos=2) 

result1=rbind(c("Characteristics", NA, NA, NA, "HR(95%CI)","p"),result)

forestplot(result1[,c(1,5,6)], 
           mean=result1[,2],   
           lower=result1[,3],  
           upper=result1[,4], 
           zero=1,            
           boxsize=0.3,      
           clip=c(1,10),graph.pos= 2, 
           line.margin=unit(8, 'mm'),
           colgap=unit(6, 'mm'),
           fn.ci_norm="fpDrawDiamondCI", 
           col=fpColors(box ='#021eaa', 
                        lines ='#021eaa'),
           txt_gp=fpTxtGp(
             label=gpar(cex=1),
             ticks=gpar(cex=1), 
             xlab=gpar(cex=1)),
           lwd.zero=1,
           lwd.ci=2,
           lwd.xaxis=1, 
           lty.ci=1,
           ci.vertices =F,
           ineheight=unit(5, 'mm'), 
           xlab="Hazzard Ratio",#graphwidth = unit(.25,"npc"),
           hrzl_lines=list("1" = gpar(lty=1,lwd=2),
                           "2" = gpar(lty=2),
                           "8"= gpar(lwd=2,lty=1,columns=c(1:4))),
           xticks=c(0,1,3,5,7,10))
