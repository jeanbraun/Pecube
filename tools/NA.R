folder=readline("Enter run name (RUN00)> ")
fnme=paste(folder,"/NA/NA_results.csv",sep="")
if (file.exists(fnme)) {cat(fnme," exists","\n")} else{fnme=readline("Can't find input file, please enter new NA file name > ")}
x=read.table(fnme,sep=',',header=TRUE)
np=length(x)-1
mmax=median(x[,1])
option=1
title=c("param1","param2","param3","param4","param5","param6","param7","param8","param9","param10")
cat("Number of parameters detected: ",np,"\n")
pdff=FALSE

while (option>=1) {

cat("------------------------------------------------------------","\n")
cat("Options:","\n")
cat("0) Exit","\n")
cat("1) Show misfit evolution","\n")
cat("2) Show parameters evolution","\n")
cat("3) Show misfit evolution for each parameter","\n")
cat("4) Show misfit in binary plots","\n")
cat("5) Rename parameters ","\n")
cat("6) Change input file name - Current value is ",fnme,"\n")
if (pdff) cat("7) Turn PDF writing off","\n") else cat("7) Turn PDF writing on","\n")
nr=length(x[,1])
cat('Current number of forward runs > ',nr)
option=readline("Enter option > ")

# first bit of code to see misfit evolution

if (option=="1"){
cv=rainbow(256,start=0,end=4/6)
x=read.table(fnme,sep=',',header=TRUE)
mmax=median(x[,1])
par(mfrow=c(1,1),mar=c(5,5,5,5))
cl=pmin(1,(x[,1]-min(x[,1]))/(mmax-min(x[,1])))
plot(x[,1],pch=16,col=cv[1+255*cl],ylim=c(0,mmax),xlab="iteration",ylab="misfit")
if (pdff) {
pdf(paste(folder,"/NA/Misfit_evolution.pdf",sep=""))
par(mfrow=c(1,1),mar=c(5,5,5,5))
plot(x[,1],pch=16,col=cv[1+255*cl],ylim=c(0,mmax),xlab="iteration",ylab="misfit")
dev.off()
}
}

# second bit of code to see parameter evolution

if (option=="2") {
cv=rainbow(256,start=0,end=4/6)
x=read.table(fnme,sep=',',header=TRUE)
mmax=median(x[,1])
par(mfrow=c((np+1)/2,2),mar=c(2,2,2,2))
cl=pmin(1,(x[,1]-min(x[,1]))/(mmax-min(x[,1])))
for (i in 1:np){plot(x[,1+i],pch=16,col=cv[1+255*cl],xlab="iteration",ylab=title[i])}
if (pdff) {
pdf(paste(folder,"/NA/Parameters_evolution.pdf",sep=""),onefile=TRUE)
par(mfrow=c(np/2,2),mar=c(2,2,2,2))
for (i in 1:np){plot(x[,1+i],pch=16,col=cv[1+255*cl],xlab="iteration",ylab=title[i])}
dev.off()
}
}

# third bit of code to see how misfit varies with each parameter

if (option=="3") {
cv=rainbow(256,start=0,end=4/6)
x=read.table(fnme,sep=',',header=TRUE)
mmax=median(x[,1])
par(mfrow=c((np+1)/2,2),mar=c(2,2,2,2))
cl=pmin(1,(x[,1]-min(x[,1]))/(mmax-min(x[,1])))
for (i in 1:np){plot(x[,i+1],x[,1],ylim=c(0,mmax),ylab="misfit",xlab=title[i],pch=16,col=cv[1+255*cl])}
if (pdff) {
pdf(paste(folder,"/NA/Misfit_vs_Parameters.pdf",sep=""),onefile=TRUE)
par(mfrow=c(np/2,2),mar=c(2,2,2,2))
for (i in 1:np){plot(x[,i+1],x[,1],ylim=c(0,mmax),ylab="misfit",xlab=title[i],pch=16,col=cv[1+255*cl])}
dev.off()
}
}

# fourth bit of code to see how misfit varies per couple of parameters

if (option=="4") {
cv=rainbow(256,start=0,end=4/6)
x=read.table(fnme,sep=',',header=TRUE)
mmax=median(x[,1])
par(mfrow=c(np,np),mar=c(1,1,1,1))
cl=pmin(1,(x[,1]-min(x[,1]))/(mmax-min(x[,1])))
or=order(x[,1],decreasing=TRUE)
for (j in 1:np){for (i in 1:np) {plot(x[or,i+1],x[or,j+1],ylab=title[j],xlab=title[i],pch=16,col=cv[1+255*cl[or]])}}
if (pdff) {
pdf(paste(folder,"/NA/Parameters_2x2.pdf",sep=""),onefile=TRUE)
par(mfrow=c(np,np),mar=c(1,1,1,1))
for (j in 1:np){for (i in 1:np) {plot(x[or,i+1],x[or,j+1],ylab=title[j],xlab=title[i],pch=16,col=cv[1+255*cl[or]])}}
dev.off()
}
nmin=or[length(x[,1])]
cat("index of lowest misfit: ",nmin,"\n")
cat("lowest misfit: ",x[nmin,1],"\n")
for (i in 1:np){cat(title[i],":",x[nmin,i+1],"\n")}
}

if (option=="5") {
for(i in 1:np){title[i]=readline(paste("Title",i,": "))}
}

if (option=="6") {
fnme=readline("Enter input file name > ")
x=read.table(fnme,sep=',',header=TRUE)
np=length(x)-1
mmax=median(x[,1])
title=c("param1","param2","param3","param4","param5","param6","param7","param8","param9","param10")
cat("Number of parameters detected: ",np,"\n")
}

if (option=="7") {
pdff = ! pdff
}

}

