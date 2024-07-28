load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")
plot(System$Series$Time,System$Series$Termination)

load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_temp.RData")
lines(System$Series$Time,System$Series$Termination)
abline(v = Time_Critical, col = "red", lty = 2)

load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_temp3.RData")
lines(System$Series$Time,System$Series$Termination)







load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")
System$Series$PDI  = System$Series$X_w /System$Series$X_n
plot(System$Series$P[1:length(System$Series$PDI)],System$Series$PDI,type="l",xlab="Conv. (%)",ylab="PDI",main="Monomer conversion v.s. PDI",ylim=c(1,2),lwd=2)

load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_temp3.RData")
System$Series$PDI  = System$Series$X_w /System$Series$X_n
lines(System$Series$P[1:length(System$Series$PDI)],System$Series$PDI,lwd=2,col="red")

load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical2.RData")
System$Series$PDI  = System$Series$X_w /System$Series$X_n
lines(System$Series$P[1:length(System$Series$PDI)],System$Series$PDI,lwd=2, col="blue")

legend("bottomright",c("Relaxation(SSA)","SSA","Quasi-leaping"),lty=c(1,1,1),col=c("black","blue","red"),lwd=c(2,2,2))

##############


layout(1)
load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Leading_Leaping.RData")
plot(System$Series$Time,System$Series$P*100,type="l",lwd=2,col="black",xlab="Time",ylab="Conv. (%)",main="Time v.s. Monomer conversion",ylim=c(0,100))
logistic_model <- nls(P*100 ~ SSlogis(Time, Asym, xmid, scal), data = System$Series)
predicted_time <- seq(0, 21000, length.out = 200)
predicted_conversion <- predict(logistic_model, newdata = data.frame(Time = predicted_time))
lines(predicted_time, predicted_conversion, col = "blue", lwd = 2)

abline(h = 80, col = "blue", lty = 2)
load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")
sample_indices_FRP <- 0:80*100+1
sample_indices_RAFT_No <- 0:80*100+1
### Time v.s. Conv. ###
lines(0:10000,(1-exp(sqrt(4*k_p^2/(k_I*k_T)*I)*(exp(-k_I*0:10000/2)-1)))*100,
      lwd=2,type="l",col="orange")
points(System.FRP$Series$Time[sample_indices_FRP],System.FRP$Series$P[sample_indices_FRP]*100,col="blue",pch=3)
points(System.RAFT.No$Series$Time[sample_indices_RAFT_No],System.RAFT.No$Series$P[sample_indices_RAFT_No]*100,col="red",pch=4)
legend("bottomright",c("Theoretical","FRP","RAFT.No","RAFT.IS"),col=c("orange","blue","red","black")
       ,lty=c(1,0,0,1),lwd=c(1,1,1,2),pch=c(NA,3,4),cex=1.3)


load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")
System$Series$PDI  = System$Series$X_w /System$Series$X_n
plot(System$Series$P[1:1000]*100,c(System$Series$PDI,rep(NA,1000-length(System$Series$PDI))),type="l",lwd=2,col="blue",xlab="Conv. (%)",ylab="PDI",main="Time v.s. Monomer conversion", ylim=c(1,2))

load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_temp3.RData")
System$Series$PDI  = System$Series$X_w /System$Series$X_n
lines(System$Series$P[1:64]*100,System$Series$PDI[1:64],lwd=2)
lines(System$Series$P[64:1000]*100,System$Series$PDI[64:1000],lwd=2,col="red",lty=2)
abline(v = P_Relaxation*100, col = "red", lty = 2)

load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_temp.RData")
System$Series$PDI  = System$Series$X_w /System$Series$X_n
lines(System$Series$P[64:length(System$Series$PDI)]*100,System$Series$PDI[64:length(System$Series$PDI)],lwd=2, col="red",lty=1)

mtext(paste(round(P_Relaxation, 5)*100,"%  →"), side=3, adj=0, cex=0.8, col="red")
legend("bottomright",c("Relaxation(SSA)","SSA","Leaping_attempt_1","Leaping_attempt_2"),lty=c(1,1,2,1),col=c("black","blue","red","red"),lwd=c(2,2,2,2),cex=1.8)

load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Leading_Leaping.RData")
plot(System$Series$Time,System$Series$Radical,
     col=gray(0.6),xlab="Time",ylab="[R·]",
     main="Time v.s. [R·]")
abline(v = Time_Critical, col = "red", lty = 2)
# Smooth.line = filter(System$Series$Radical, sides=2, rep(1,100)/100)
loess_model = loess(Radical ~ Time, data = System$Series, span = 0.1)
Smooth.line = predict(loess_model, System$Series$Time)
lines(System$Series$Time,
      Smooth.line,col="blue")
lines(System$Series$Time,
      sqrt(k_I/k_T*I*exp(-k_I*System$Series$Time)),
      lwd=2,col="red",lty=2)

loess_model <- loess(Radical ~ Time, data = System$Series, span = 0.1)
smoothed_values <- predict(loess_model, System$Series$Time)
lines(System$Series$Time,smoothed_values)


plot(Leading_System$Series$Time,log(1/(1-Leading_System$Series$P)),type="l",lwd=2,col="black",xlab="Time",ylab="Conv.",main="Time v.s. Conv.",ylim=c(0,1))
points(System.FRP$Series$Time[sample_indices_FRP],log(1/(1-System.FRP$Series$P[sample_indices_FRP])),col="blue",pch=3)
