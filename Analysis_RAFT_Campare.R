# load & set ###################################################################
rm(list = ls())
options(max.print = 20)
load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical_7.05.RData")
# System.RAFT = System
# layout(matrix(c(1,1,2,2,7,7,1,1,2,2,8,8,3,4,5,6,9,9),6,3))
layout(matrix(c(1,2,5,1,2,6,3,4,7),3,3))
par(mar=c(4,4,2,2))

###
# test = c(Radical[Radical!=0],CTAgent[CTAgent!=0],Termination[Termination!=0])
# (sd(test)^2+mean(test)^2)/mean(test)
# (sd(test)^2)/mean(test)+mean(test)
# mean(test^2)/mean(test)
# System.RAFT.No$Series$X_w[8001]
###
# sqrt(k_I/k_T*I+k_CTAgent^2/(16*k_T^2)*CTA^2)-k_CTAgent/(4*k_T)*CTA
###

# GPC Function #################################################################
transparent_colors = sapply(1:8, function(i) {
  rgb_val = col2rgb(i)
  rgb(rgb_val[1], rgb_val[2], rgb_val[3], maxColorValue = 255, alpha = 255*0.5)  # 设置 alpha = 255*0.7 (70% 透明度)
})

k1 = 1
k2 = 1
Gap = 2
Critical = c("10 %","20 %","30 %","40 %","50 %","60 %","70 %","80 %")


GPC = function(Sys){
  i = Critical[8]
  Weight = Sys[[i]]$Data
  
  Max_MW = max(which(Weight != 0))
  Breaks = c(seq(0, Max_MW, by=Gap))
  Intervals = cut(1:Max_MW, breaks=Breaks, include.lowest = TRUE, right = TRUE)
  Temp = Weight[1:Max_MW]*1:Max_MW
  Weight_intervals = as.numeric(tapply(Temp, Intervals, sum))
  Numbers = as.numeric(tapply(Weight[1:Max_MW], Intervals, sum))
  GPC_Result = matrix(0, ncol = 2, nrow = length(Numbers))
  non_zero_indices <- which(Weight_intervals != 0)
  GPC_Result[non_zero_indices, 1] <- log10(Weight_intervals[non_zero_indices] / Numbers[non_zero_indices]) * k1 + k2
  GPC_Result[non_zero_indices, 2] <- (-0.4228 * GPC_Result[non_zero_indices, 1] + 10.38) * Weight_intervals[non_zero_indices]
  GPC_Result = GPC_Result[apply(GPC_Result != 0, 1, any), , drop = FALSE]
  # GPC_Result[,2] = GPC_Result[,2]/max(GPC_Result[,2])
  plot(GPC_Result[,1],GPC_Result[,2],type="l",col=8,main="GPC",lwd=2,xlim = rev(range(GPC_Result[,1])))
  
  for (i in Critical[-8]){
    Weight = Sys[[i]]$Data
    Max_MW = max(which(Weight != 0))
    Breaks = c(seq(0, Max_MW, by=Gap))
    Intervals = cut(1:Max_MW, breaks=Breaks, include.lowest = TRUE, right = TRUE)
    Temp = Weight[1:Max_MW]*1:Max_MW
    Weight_intervals = as.numeric(tapply(Temp, Intervals, sum))
    Numbers = as.numeric(tapply(Weight[1:Max_MW], Intervals, sum))
    GPC_Result = matrix(0, ncol = 2, nrow = length(Numbers))
    non_zero_indices <- which(Weight_intervals != 0)
    GPC_Result[non_zero_indices, 1] <- log10(Weight_intervals[non_zero_indices] / Numbers[non_zero_indices]) * k1 + k2
    GPC_Result[non_zero_indices, 2] <- (-0.4228 * GPC_Result[non_zero_indices, 1] + 10.38) * Weight_intervals[non_zero_indices]
    GPC_Result = GPC_Result[apply(GPC_Result != 0, 1, any), , drop = FALSE]
    # GPC_Result[,2] = GPC_Result[,2]/max(GPC_Result[,2])
    lines(GPC_Result[,1],GPC_Result[,2],col=which(i == Critical),lwd=2)
  }
  legend("topright",Critical,col=1:8,border = NA,lwd=2,ncol = 2, cex=1.3)
}

#################################################################################################################################
### FRP #########################################################################################################################
#################################################################################################################################
### One ########################################################################
### Distribution ###
nonzero_indices = which(System.FRP$Critical$`80 %`$Data != 0)
plot(1:5000,System.FRP$Critical$`80 %`$Data,#filter(System.FRP$Critical$`80 %`$Data, sides=2, rep(1,10)/10),
     col=1,lwd=2,type="l",xlab="Chain Length",ylab="Amount",main="Chain Length Distribution--FRP",
     xlim = rev(range(1:nonzero_indices[length(nonzero_indices)])),
     ylim = range(System.FRP$Critical$`80 %`$Data)
)
polygon(c(1:5000, rev(1:5000)), 
        c(System.FRP$Critical$`80 %`$Termination[-1], 
          rep(0, length(System.FRP$Critical$`80 %`$Termination[-1]))), 
        col = transparent_colors[8],border = NA)
for (i in rev(Critical[1:7])){
  lines(1:5000,System.FRP$Critical[[i]]$Data,
        # filter(System.FRP$Critical[[i]]$Data, sides=2, rep(1,10)/10),
        col=which(i == Critical),lwd=2)
  polygon(c(1:5000, rev(1:5000)), 
          c(System.FRP$Critical[[i]]$Termination[-1], 
            rep(0, length(System.FRP$Critical[[i]]$Termination[-1]))), 
          col = transparent_colors[which(i == Critical)],border = NA)
}
legend("topleft",Critical,col=1:8,fill=transparent_colors[1:8],lwd=2,ncol = 2,cex=1.3)

### GPC ###
GPC(System.FRP$Critical)

nonzero_indices = which(System.FRP$Critical$`80 %`$Data != 0)
plot(1:5000,System.FRP$Critical$`80 %`$Data/sum(System.FRP$Critical$`80 %`$Data),#filter(System.FRP$Critical$`80 %`$Data, sides=2, rep(1,10)/10),
     col=1,lwd=2,type="l",xlab="Chain Length",ylab="Rate",main="Relative Frequency Distribution",
     xlim = rev(range(1:1300)),
     ylim = range(System.FRP$Critical$`80 %`$Data/sum(System.FRP$Critical$`80 %`$Data))
)
# polygon(c(1:5000, rev(1:5000)), 
#         c(System.FRP$Critical$`80 %`$Termination[-1], 
#           rep(0, length(System.FRP$Critical$`80 %`$Termination[-1]))), 
#         col = transparent_colors[8],border = NA)
for (i in rev(Critical[1:7])){
  lines(1:5000,System.FRP$Critical[[i]]$Data/sum(System.FRP$Critical[[i]]$Data),
        # filter(System.FRP$Critical[[i]]$Data, sides=2, rep(1,10)/10),
        col=which(i == Critical),lwd=2)
  # polygon(c(1:5000, rev(1:5000)), 
  #         c(System.FRP$Critical[[i]]$Termination[-1], 
  #           rep(0, length(System.FRP$Critical[[i]]$Termination[-1]))), 
  #         col = transparent_colors[which(i == Critical)],border = NA)
}
legend("topleft",Critical,col=1:8,lwd=2,ncol = 2,cex=1.3)

nonzero_indices = which(System.FRP$Critical$`80 %`$Data != 0)
plot(1:5000,System.FRP$Critical$`80 %`$Data,#filter(System.FRP$Critical$`80 %`$Data, sides=2, rep(1,10)/10),
     col=1,lwd=2,type="n",xlab="Chain Length",ylab="Amount",main="Termination",
     xlim = rev(range(1:1500)),
     ylim = range(System.FRP$Critical$`80 %`$Data)
)
polygon(c(1:5000, rev(1:5000)), 
        c(System.FRP$Critical$`80 %`$Termination[-1], 
          rep(0, length(System.FRP$Critical$`80 %`$Termination[-1]))), 
        col = transparent_colors[8],border = NA)
for (i in rev(Critical[1:7])){
  # lines(1:5000,System.FRP$Critical[[i]]$Data,
  #       # filter(System.FRP$Critical[[i]]$Data, sides=2, rep(1,10)/10),
  #       col=which(i == Critical),lwd=2)
  polygon(c(1:5000, rev(1:5000)), 
          c(System.FRP$Critical[[i]]$Termination[-1], 
            rep(0, length(System.FRP$Critical[[i]]$Termination[-1]))), 
          col = transparent_colors[which(i == Critical)],border = NA)
}
legend("topleft",Critical,col=1:8,fill=transparent_colors[1:8],border = NA,ncol = 2,cex=1.3)
# ### Two ########################################################################
# ### Time v.s. [R·] ###
# plot(System.FRP$Series$Time,System.FRP$Series$Radical,
#      col=gray(0.6),xlab="Time",ylab="[R·]",
#      main="Time v.s. [R·]")
# # Smooth.line = filter(System.FRP$Series$Radical, sides=2, rep(1,100)/100)\
# loess_model = loess(Radical ~ Time, data = System.FRP$Series, span = 0.1)
# Smooth.line = predict(loess_model, System.FRP$Series$Time)
# lines(System.FRP$Series$Time,
#       Smooth.line,col="blue")
# lines(System.FRP$Series$Time,
#       sqrt(k_I/k_T*I*exp(-k_I*System.FRP$Series$Time)),
#       lwd=2,col="red",lty=2)
# legend("bottomright",c("Simulation","Smoothed","Theorical"),
#        col=c(gray(0.6),"blue","red"),lty=c(0,1,2),lwd=c(0,1,2),pch=c(1,NA,NA),cex=0.5)
# 
# ### Time v.s. [T] ###
# plot(1:(max(System.FRP$Series$Time)/250)*250,I*(1-exp(-k_I*1:(max(System.FRP$Series$Time)/250)*250)),col=gray(0.8),
#      lwd=1,pch=3,
#      xlab="Time",ylab="[T]",
#      main="Time v.s. [T]")
# lines(System.FRP$Series$Time,System.FRP$Series$Termination,col="blue")
# points(1:(max(System.FRP$Series$Time)/250)*250,M*(1-exp(-k_I*1:(max(System.FRP$Series$Time)/250)*250)),col=gray(0.8),
#        lwd=1,pch=3)
# legend("bottomright",c("Theorical","FRP"),col=c(col=gray(0.8),"blue")
#        ,lty=c(0,1,1),lwd=c(1,1,1),pch=c(3,NA,NA),cex=0.5)
# 
# ### Conv. v.s. [R·] ###
# plot(System.FRP$Series$P,System.FRP$Series$Radical,
#      col=gray(0.6),xlab="Conv.",ylab="[R·]",
#      main="Conv. v.s. [R·]")
# # Smooth.line = filter(System.FRP$Series$Radical, sides=2, rep(1,100)/100)
# loess_model = loess(Radical ~ P, data = System.FRP$Series, span = 0.1)
# Smooth.line = predict(loess_model, System.FRP$Series$P)
# lines(System.FRP$Series$P,Smooth.line,col="blue")
# 
# ### Conv. v.s. [T] ###
# plot(System.FRP$Series$P,System.FRP$Series$Termination,
#      col="blue",type="l",
#      xlab="Conv.",ylab="[R·]",main="Conv. v.s. [T]")

### Three ######################################################################
### X_n ###
plot(System.FRP$Series$P*100,System.FRP$Series$X_n,
     col="blue",type="l",
     xlab="Conv. (%)",ylab="X_n",main="Monomer conversion v.s. Number average DP")

### X_w ###
plot(System.FRP$Series$P*100,System.FRP$Series$X_w,
     col="blue",type="l",
     xlab="Conv. (%)",ylab="X_w",main="Monomer conversion v.s. Weight average DP")

### PDI ###
plot(System.FRP$Series$P*100,System.FRP$Series$PDI,
     col="blue",type="l",ylim=c(1,2),
     xlab="Conv. (%)",ylab="PDI",main="Monomer conversion v.s. PDI")

#################################################################################################################################
### RAFT without Intermediate ###################################################################################################
#################################################################################################################################
### One ########################################################################
### Distribution ###
nonzero_indices = which(System.RAFT.No$Critical$`80 %`$Data != 0)
plot(1:500,System.RAFT.No$Critical$`80 %`$Data,
     col=1,lwd=2,type="l",xlab="Chain Length",ylab="Amount",main="Chain Length Distribution--RAFT NO Intermediate",
     xlim = rev(range(1:nonzero_indices[length(nonzero_indices)])),
     ylim = range(System.RAFT.No$Critical$`10 %`$Data)
)
polygon(c(1:500, rev(1:500)), 
        c(System.RAFT.No$Critical$`80 %`$Termination[-1], 
          rep(0, length(System.RAFT.No$Critical$`80 %`$Termination[-1]))), 
        col = transparent_colors[8])
for (i in rev(Critical[1:7])){
  lines(1:500,
        System.RAFT.No$Critical[[i]]$Data,
        col=which(i == Critical),lwd=2)
  polygon(c(1:500, rev(1:500)), 
          c(System.RAFT.No$Critical[[i]]$Termination[-1], 
            rep(0, length(System.RAFT.No$Critical[[i]]$Termination[-1]))), 
          col = transparent_colors[which(i == Critical)])
}
legend("topleft",Critical,col=1:8,fill=transparent_colors[1:8],lwd=2,ncol = 2,cex=1.3)

### GPC ###
GPC(System.RAFT.No$Critical)

nonzero_indices = which(System.RAFT.No$Critical$`80 %`$Data != 0)
plot(1:500,System.RAFT.No$Critical$`80 %`$Data/sum(System.RAFT.No$Critical$`80 %`$Data),
     col=1,lwd=2,type="l",xlab="Chain Length",ylab="Rate",main="Relative Frequency Distribution",
     xlim = rev(range(1:nonzero_indices[length(nonzero_indices)])),
     ylim = range(System.RAFT.No$Critical$`10 %`$Data/sum(System.RAFT.No$Critical$`10 %`$Data))
)
# polygon(c(1:500, rev(1:500)), 
#         c(System.RAFT.No$Critical$`80 %`$Termination[-1], 
#           rep(0, length(System.RAFT.No$Critical$`80 %`$Termination[-1]))), 
#         col = transparent_colors[8])
for (i in rev(Critical[1:7])){
  lines(1:500,
        System.RAFT.No$Critical[[i]]$Data/sum(System.RAFT.No$Critical[[i]]$Data),
        col=which(i == Critical),lwd=2)
  # polygon(c(1:500, rev(1:500)), 
  #         c(System.RAFT.No$Critical[[i]]$Termination[-1], 
  #           rep(0, length(System.RAFT.No$Critical[[i]]$Termination[-1]))), 
  #         col = transparent_colors[which(i == Critical)])
}
legend("topleft",Critical,col=1:8,lwd=2,ncol = 2,cex=1.3)

nonzero_indices = which(System.RAFT.No$Critical$`80 %`$Data != 0)
plot(1:500,System.RAFT.No$Critical$`80 %`$Data,
     col=1,lwd=2,type="n",xlab="Chain Length",ylab="Amount",main="Termination",
     xlim = rev(range(1:nonzero_indices[length(nonzero_indices)])),
     ylim = range(System.RAFT.No$Critical$`80 %`$Termination)
)
polygon(c(1:500, rev(1:500)), 
        c(System.RAFT.No$Critical$`80 %`$Termination[-1], 
          rep(0, length(System.RAFT.No$Critical$`80 %`$Termination[-1]))), 
        col = transparent_colors[8],border = NA)
for (i in rev(Critical[1:7])){
  # lines(1:500,
  #       System.RAFT.No$Critical[[i]]$Data,
  #       col=which(i == Critical),lwd=2)
  polygon(c(1:500, rev(1:500)), 
          c(System.RAFT.No$Critical[[i]]$Termination[-1], 
            rep(0, length(System.RAFT.No$Critical[[i]]$Termination[-1]))), 
          col = transparent_colors[which(i == Critical)],border = NA)
}
legend("topleft",Critical,col=1:8,fill=transparent_colors[1:8],border = NA,ncol = 2,cex=1.1)
# ### Two ########################################################################
# ### Time v.s. [R·] ###
# plot(System.RAFT.No$Series$Time,System.RAFT.No$Series$Radical,
#      col=gray(0.6),xlab="Time",ylab="[R·]",
#      main="Time v.s. [R·]")
# # Smooth.line = filter(System.RAFT.No$Series$Radical, sides=2, rep(1,100)/100)
# loess_model = loess(Radical ~ Time, data = System.RAFT.No$Series, span = 0.1)
# Smooth.line = predict(loess_model, System.RAFT.No$Series$Time)
# lines(System.RAFT.No$Series$Time,
#       Smooth.line,col="blue")
# lines(System.RAFT.No$Series$Time,
#       sqrt(k_I/k_T*I*exp(-k_I*System.RAFT.No$Series$Time)),
#       lwd=2,col="red",lty=2)
# legend("bottomright",c("Simulation","Smoothed","Theorical"),
#        col=c(gray(0.6),"blue","red"),lty=c(0,1,2),lwd=c(0,1,2),pch=c(1,NA,NA),cex=0.5)
# 
# ### Time v.s. [T] ###
# plot(1:(max(System.RAFT.No$Series$Time)/250)*250,I*(1-exp(-k_I*1:(max(System.RAFT.No$Series$Time)/250)*250)),col=gray(0.8),
#      lwd=1,pch=3,
#      xlab="Time",ylab="[T]",
#      main="Time v.s. [T]")
# lines(System.RAFT.No$Series$Time,System.RAFT.No$Series$Termination,col="blue")
# points(1:(max(System.RAFT.No$Series$Time)/250)*250,M*(1-exp(-k_I*1:(max(System.RAFT.No$Series$Time)/250)*250)),col=gray(0.8),
#        lwd=1,pch=3)
# legend("bottomright",c("Theorical","FRP"),col=c(col=gray(0.8),"blue")
#        ,lty=c(0,1,1),lwd=c(1,1,1),pch=c(3,NA,NA),cex=0.5)
# 
# ### Conv. v.s. [R·] ###
# plot(System.RAFT.No$Series$P,System.RAFT.No$Series$Radical,
#      col=gray(0.6),xlab="Conv.",ylab="[R·]",
#      main="Conv. v.s. [R·]")
# # Smooth.line = filter(System.RAFT.No$Series$Radical, sides=2, rep(1,100)/100)
# loess_model = loess(Radical ~ P, data = System.RAFT.No$Series, span = 0.1)
# Smooth.line = predict(loess_model, System.RAFT.No$Series$P)
# lines(System.RAFT.No$Series$P,Smooth.line,col="blue")
# 
# ### Conv. v.s. [T] ###
# plot(System.RAFT.No$Series$P,System.RAFT.No$Series$Termination,
#      col="blue",type="l",
#      xlab="Conv.",ylab="[R·]",main="Conv. v.s. [T]")

### Three ######################################################################
### X_n ###
plot(System.RAFT.No$Series$P*100,System.RAFT.No$Series$X_n,
     col="blue",type="l",
     xlab="Conv. (%)",ylab="X_n",main="Monomer conversion v.s. Number average DP")

### X_w ###
plot(System.RAFT.No$Series$P*100,System.RAFT.No$Series$X_w,
     col="blue",type="l",
     xlab="Conv. (%)",ylab="X_w",main="Monomer conversion v.s. Weight average DP")

### PDI ###
plot(System.RAFT.No$Series$P*100,System.RAFT.No$Series$PDI,
     col="blue",type="l",ylim=c(1,2),
     xlab="Conv. (%)",ylab="PDI",main="Monomer conversion v.s. PDI")


### Four #######################################################################
### Distribution ###
nonzero_indices = which(System.RAFT.No$Critical$`80 %`$Data != 0)
plot(1:500,System.RAFT.No$Critical$`80 %`$Data,
     col=1,lwd=2,type="l",xlab="Chain Length",ylab="amount",main="Chain Length Distribution--RAFT NO Intermediate",
     xlim = rev(range(1:nonzero_indices[length(nonzero_indices)])),
     ylim = range(System.RAFT.No$Critical$`10 %`$Data)
)
polygon(c(1:500, rev(1:500)), 
        c(System.RAFT.No$Critical$`80 %`$Termination[-1], 
          rep(0, length(System.RAFT.No$Critical$`80 %`$Termination[-1]))), 
        col = transparent_colors[8])
for (i in rev(Critical[1:7])){
  lines(1:500,
        System.RAFT.No$Critical[[i]]$Data,
        col=which(i == Critical),lwd=2)
  polygon(c(1:500, rev(1:500)), 
          c(System.RAFT.No$Critical[[i]]$Termination[-1], 
            rep(0, length(System.RAFT.No$Critical[[i]]$Termination[-1]))), 
          col = transparent_colors[which(i == Critical)])
}
legend("topleft",Critical,col=1:8,fill=transparent_colors[1:8],lwd=2,ncol = 2)

### GPC ###
GPC(System.RAFT.No$Critical)

### Time v.s. [R·] ###
plot(System.RAFT.No$Series$Time,System.RAFT.No$Series$Radical,
     col=gray(0.6),xlab="Time",ylab="[R·]",
     main="Time v.s. [R·]")
Smooth.line = filter(System.RAFT.No$Series$Radical, sides=2, rep(1,100)/100)
lines(System.RAFT.No$Series$Time,
      Smooth.line,col="blue")
lines(System.RAFT.No$Series$Time,
      sqrt(k_I/k_T*I*exp(-k_I*System.RAFT.No$Series$Time)),
      lwd=2,col="red",lty=2)
legend("topright",c("Simulation","Smoothed","Theorical"),
       col=c(gray(0.6),"blue","red"),lty=c(0,1,2),lwd=c(0,1,2),pch=c(1,NA,NA),cex=0.5)

### Time v.s. [T] ###
plot(1:(max(System.RAFT.No$Series$Time)/250)*250,I*(1-exp(-k_I*1:(max(System.RAFT.No$Series$Time)/250)*250)),col=gray(0.8),
     lwd=1,pch=3,
     xlab="Time",ylab="[T]",
     main="Time v.s. [T]")
lines(System.RAFT.No$Series$Time,System.RAFT.No$Series$Termination,col="blue")
points(1:(max(System.RAFT.No$Series$Time)/250)*250,M*(1-exp(-k_I*1:(max(System.RAFT.No$Series$Time)/250)*250)),col=gray(0.8),
       lwd=1,pch=3)
legend("bottomright",c("Theorical","FR"),col=c(col=gray(0.8),"blue")
       ,lty=c(0,1,1),lwd=c(1,1,1),pch=c(3,NA,NA),cex=0.5)

### Conv. v.s. [R·] ###
plot(System.RAFT.No$Series$P,System.RAFT.No$Series$Radical,
     col=gray(0.6),xlab="P",ylab="[R·]",
     main="Conv. v.s. [R·]")
Smooth.line = filter(System.RAFT.No$Series$Radical, sides=2, rep(1,100)/100)
lines(System.RAFT.No$Series$P ,Smooth.line,col="blue")

### Conv. v.s. [T] ###
plot(System.RAFT.No$Series$P,System.RAFT.No$Series$Termination,
     col="blue",type="l",
     xlab="Conv.",ylab="[R·]",main="Conv. v.s. [T]")

### X_n ###
plot(System.RAFT.No$Series$P,System.RAFT.No$Series$X_n2,
     col="blue",type="l",
     xlab="Conv.",ylab="X_n",main="Conv. v.s. X_n")

### X_w ###
plot(System.RAFT.No$Series$P,System.RAFT.No$Series$X_w2,
     col="blue",type="l",
     xlab="Conv.",ylab="X_w",main="Conv. v.s. X_w")

### PDI ###
plot(System.RAFT.No$Series$P,System.RAFT.No$Series$PDI2,
     col="blue",type="l",
     xlab="Conv.",ylab="PDI",main="Conv. v.s. PDI")

#################################################################################################################################
### RAFT with Intermediate ######################################################################################################
#################################################################################################################################
# Distribution #################################################################
System$Series$PDI  = System$Series$X_w /System$Series$X_n
System.RAFT = System
### One ########################################################################
### Distribution ###
nonzero_indices = which(System.RAFT$Critical$`80 %`$Data != 0)
plot(1:1000,System.RAFT$Critical$`80 %`$Data,
     col=1,lwd=2,type="l",xlab="Chain Length",ylab="amount",main="Chain Length Distribution--RAFT NO Intermediate",
     xlim = rev(range(1:200)),
     ylim = range(System.RAFT$Critical$`10 %`$Data)
)
# polygon(c(1:1000, rev(1:1000)), 
#         c(System.RAFT$Critical$`80 %`$Termination[-1], 
#           rep(0, length(System.RAFT$Critical$`80 %`$Termination[-1]))), 
#         col = transparent_colors[8])
for (i in rev(Critical[1:7])){
  lines(1:1000,
        System.RAFT$Critical[[i]]$Data,
        col=which(i == Critical),lwd=2)
  # polygon(c(1:1000, rev(1:1000)), 
  #         c(System.RAFT$Critical[[i]]$Termination[-1], 
  #           rep(0, length(System.RAFT$Critical[[i]]$Termination[-1]))), 
  #         col = transparent_colors[which(i == Critical)])
}
legend("topleft",Critical,col=1:8,fill=transparent_colors[1:8],lwd=2,ncol = 2)

### GPC ###
GPC(System.RAFT$Critical)

### Two ########################################################################
### Time v.s. [R·] ###
plot(System.RAFT$Series$Time,System.RAFT$Series$Radical,
     col=gray(0.6),xlab="Time",ylab="[R·]",
     main="Time v.s. [R·]")
abline(v = Time_Critical, col = "red", lty = 2)
Smooth.line = filter(System.RAFT$Series$Radical, sides=2, rep(1,100)/100)
lines(System.RAFT$Series$Time,
      Smooth.line,col="blue")
lines(System.RAFT$Series$Time,
      sqrt(k_I/k_T*I*exp(-k_I*System.RAFT$Series$Time)),
      lwd=2,col="red",lty=2)
legend("topright",c("Simulation","Smoothed","Theorical"),
       col=c(gray(0.6),"blue","red"),lty=c(0,1,2),lwd=c(0,1,2),pch=c(1,NA,NA),cex=0.5)
mtext(round(Time_Critical, 5), side=3, adj=0, cex=0.8, col="red")
### Time v.s. [T] ###
# plot(1:(max(System.RAFT$Series$Time)/250)*250,I*(1-exp(-k_I*1:(max(System.RAFT$Series$Time)/250)*250)),col=gray(0.8),
#      lwd=1,pch=3,
#      xlab="Time",ylab="[T]",
#      main="Time v.s. [T]")
plot(System.RAFT$Series$Time,System.RAFT$Series$Termination,col="blue",xlab="Time",ylab="[T]",main="Time v.s. [T]")
abline(v = Time_Critical, col = "red", lty = 2)
points(1:(max(System.RAFT$Series$Time)/250)*250,M*(1-exp(-k_I*1:(max(System.RAFT$Series$Time)/250)*250)),col=gray(0.8),
       lwd=1,pch=3)
legend("bottomright",c("Theorical","FR"),col=c(col=gray(0.8),"blue")
       ,lty=c(0,1,1),lwd=c(1,1,1),pch=c(3,NA,NA),cex=0.5)
mtext(round(Time_Critical, 5), side=3, adj=0, cex=0.8, col="red")
### Conv. v.s. [R·] ###
plot(System.RAFT$Series$P[1:length(System.RAFT$Series$Radical)],System.RAFT$Series$Radical,
     col=gray(0.6),xlab="Conv.",ylab="[R·]",
     main="Conv. v.s. [R·]")
abline(v = P_Relaxation, col = "red", lty = 2)
Smooth.line = filter(System.RAFT$Series$Radical, sides=2, rep(1,100)/100)
lines(System.RAFT$Series$P[1:length(System.RAFT$Series$Radical)],Smooth.line,col="blue")
mtext(round(P_Relaxation, 5), side=3, adj=0, cex=0.8, col="red")

### Conv. v.s. [T] ###
plot(System.RAFT$Series$P[1:length(System.RAFT$Series$Termination)],System.RAFT$Series$Termination,
     col="blue",type="l",
     xlab="Conv.",ylab="[R·]",main="Conv. v.s. [T]")
abline(v = P_Relaxation, col = "red", lty = 2)
mtext(round(P_Relaxation, 5), side=3, adj=0, cex=0.8, col="red")
### Three ######################################################################
### X_n ###
plot(System.RAFT$Series$P[1:length(System.RAFT$Series$X_n)],System.RAFT$Series$X_n,
     col="blue",type="l",
     xlab="Conv.",ylab="X_n",main="Conv. v.s. X_n")
abline(v = P_Relaxation, col = "red", lty = 2)
### X_w ###
plot(System.RAFT$Series$P[1:length(System.RAFT$Series$X_w)],System.RAFT$Series$X_w,
     col="blue",type="l",
     xlab="Conv.",ylab="X_w",main="Conv. v.s. X_w")
abline(v = P_Relaxation, col = "red", lty = 2)
### PDI ###
plot(System.RAFT$Series$P[1:length(System.RAFT$Series$PDI)],System.RAFT$Series$PDI,
     col="blue",type="l",
     xlab="Conv.",ylab="PDI",main="Conv. v.s. PDI")
abline(v = P_Relaxation, col = "red", lty = 2)

### Four #######################################################################
### Distribution ###
plot(1:500,System.RAFT$Critical$`80 %`$Data,
     col=1,lwd=2,type="l",xlab="Chain Length",ylab="amount",main="Chain Length Distribution--RAFT NO Intermediate",
     xlim = rev(range(1:200)),
     ylim = range(System.RAFT$Critical$`10 %`$Data)
)
polygon(c(1:500, rev(1:500)), 
        c(System.RAFT$Critical$`80 %`$Termination[-1], 
          rep(0, length(System.RAFT$Critical$`80 %`$Termination[-1]))), 
        col = transparent_colors[8])
for (i in rev(Critical[1:7])){
  lines(1:500,
        System.RAFT$Critical[[i]]$Data,
        col=which(i == Critical),lwd=2)
  polygon(c(1:500, rev(1:500)), 
          c(System.RAFT$Critical[[i]]$Termination[-1], 
            rep(0, length(System.RAFT$Critical[[i]]$Termination[-1]))), 
          col = transparent_colors[which(i == Critical)])
}
legend("topleft",Critical,col=1:8,fill=transparent_colors[1:8],lwd=2,ncol = 2)

### GPC ###
GPC(System.RAFT$Critical)

### Time v.s. [R·] ###
plot(System.RAFT$Series$Time,System.RAFT$Series$Radical,
     col=gray(0.6),xlab="Time",ylab="[R·]",
     main="Time v.s. [R·]")
Smooth.line = filter(System.RAFT$Series$Radical, sides=2, rep(1,100)/100)
lines(System.RAFT$Series$Time,
      Smooth.line,col="blue")
lines(System.RAFT$Series$Time,
      sqrt(k_I/k_T*I*exp(-k_I*System.RAFT$Series$Time)),
      lwd=2,col="red",lty=2)
legend("topright",c("Simulation","Smoothed","Theorical"),
       col=c(gray(0.6),"blue","red"),lty=c(0,1,2),lwd=c(0,1,2),pch=c(1,NA,NA),cex=0.5)

### Time v.s. [T] ###
plot(1:(max(System.RAFT$Series$Time)/250)*250,I*(1-exp(-k_I*1:(max(System.RAFT$Series$Time)/250)*250)),col=gray(0.8),
     lwd=1,pch=3,
     xlab="Time",ylab="[T]",
     main="Time v.s. [T]")
lines(System.RAFT$Series$Time,System.RAFT$Series$Termination,col="blue")
points(1:(max(System.RAFT$Series$Time)/250)*250,M*(1-exp(-k_I*1:(max(System.RAFT$Series$Time)/250)*250)),col=gray(0.8),
       lwd=1,pch=3)
legend("bottomright",c("Theorical","FR"),col=c(col=gray(0.8),"blue")
       ,lty=c(0,1,1),lwd=c(1,1,1),pch=c(3,NA,NA),cex=0.5)

### Conv. v.s. [R·] ###
plot(System.RAFT$Series$P,System.RAFT$Series$Radical,
     col=gray(0.6),xlab="Conv.",ylab="[R·]",
     main="Conv. v.s. [R·]")
Smooth.line = filter(System.RAFT$Series$Radical, sides=2, rep(1,100)/100)
lines(System.RAFT$Series$P,Smooth.line,col="blue")

### Conv. v.s. [T] ###
plot(System.RAFT$Series$P,System.RAFT$Series$Termination,
     col="blue",type="l",
     xlab="Conv.",ylab="[R·]",main="Conv. v.s. [T]")

### X_n ###
plot(System.RAFT$Series$P,System.RAFT$Series$X_n2,
     col="blue",type="l",
     xlab="Conv.",ylab="X_n",main="Conv. v.s. X_n")

### X_w ###
plot(System.RAFT$Series$P,System.RAFT$Series$X_w2,
     col="blue",type="l",
     xlab="Conv.",ylab="X_w",main="Conv. v.s. X_w")

### PDI ###
plot(System.RAFT$Series$P,System.RAFT$Series$PDI2,
     col="blue",type="l",
     xlab="Conv.",ylab="PDI",main="Conv. v.s. PDI")

############


sample_indices_FRP <- 0:80*100+1
sample_indices_RAFT_No <- 0:80*100+1
### Time v.s. [T] ###
plot(1:(max(System.FRP$Series$Time)/250)*250,I*(1-exp(-k_I*1:(max(System.FRP$Series$Time)/250)*250)),
     lwd=1,type="l",
     xlab="Time",ylab="[T]",
     main="Time v.s. [T]")
points(System.FRP$Series$Time[sample_indices_FRP],System.FRP$Series$Termination[sample_indices_FRP],col="blue",pch=3)
points(System.RAFT.No$Series$Time[sample_indices_RAFT_No],System.RAFT.No$Series$Termination[sample_indices_RAFT_No],col="red",pch=4)
legend("bottomright",c("Theorical","FRP","without Intermediate"),col=c("black","blue","red")
       ,lty=c(1,0,0),lwd=c(1,1,1),pch=c(NA,3,4))


layout(1)
sample_indices_FRP <- 0:80*100+1
sample_indices_RAFT_No <- 0:80*100+1
### Time v.s. Conv. ###
plot(0:10000,(1-exp(sqrt(4*k_p^2/(k_I*k_T)*I)*(exp(-k_I*0:10000/2)-1)))*100,
     lwd=1,type="l",
     xlab="Time",ylab="Conv. (%)",ylim=c(0,100),
     main="Time v.s. Monomer conversion")
abline(h=80,lty=2)
points(System.FRP$Series$Time[sample_indices_FRP],System.FRP$Series$P[sample_indices_FRP]*100,col="blue",pch=3)
points(System.RAFT.No$Series$Time[sample_indices_RAFT_No],System.RAFT.No$Series$P[sample_indices_RAFT_No]*100,col="red",pch=4)
legend("bottomright",c("Theoretical","FRP","RAFT.No"),col=c("black","blue","red")
       ,lty=c(1,0,0),lwd=c(1,1,1),pch=c(NA,3,4),cex=1.5)


# Compare ######################################################################
# Table.FR = table(System.RAFT$Critical$`80 %`$Termination)
# 
# plot(as.numeric(names(Table.FR)),Table.FR,type="l",col="blue")
# 
# plot(density(System.RAFT$Critical$`80 %`$Termination))
# barplot(Table.FR[1:1250], col = gray(0.8), border = NA, 
#         xlab="Chain Length [1:1250]", 
#         main="Chain Length Distribution for Free Radical Process")
# Smooth.line = filter(Table.FR, sides=2, rep(1,50)/30)
# abline(v = which.max(Smooth.line), col = "red", lty = 2,lwd=3)
# text(which.max(Smooth.line)+30, 
#      Smooth.line[which.max(Smooth.line)]/2, 
#      which.max(Smooth.line), pos = 3, col = "red")
# 
# Table.RAFT = table(c(System.RAFT$Termination,
#                      rowSums(System.RAFT$Intermediate),
#                      System.RAFT$CTAgent[System.RAFT$CTAgent!=0]))
# barplot(Table.RAFT[1:1250], col = gray(0.8), border = NA, 
#         xlab="Chain Length [1:1250]", 
#         main="Chain Length Distribution for RAFT Process (Only From 1:1250)")
# Smooth.line = filter(Table.RAFT, sides=2, rep(1,50)/30)
# abline(v = which.max(Smooth.line), col = "red", lty = 2,lwd=3)
# text(which.max(Smooth.line)+30, 
#      Smooth.line[which.max(Smooth.line)]/2, 
#      which.max(Smooth.line), pos = 3, col = "red")
library(ggplot2)
library(reshape2)

X_df <- melt(Intermediate[1:50,1:50])

# 绘制热力图
ggplot(X_df, aes(Var1, Var2, fill = value)) +
  geom_tile() +
  scale_fill_gradient(low = "white", high = "red") +
  labs(x = "Column", y = "Row", fill = "Value") +
  theme_minimal() +
  theme(axis.text.x = element_blank(), axis.text.y = element_blank())


library(readxl)
RAFT_Real_Data2 <- read_excel("D:/Lloyd Yin/Desktop/Dissertation/RAFT_Real_Data2.xlsx")
RAFT_Real_Data2=RAFT_Real_Data2[-4,]
Real_Conv = as.numeric(RAFT_Real_Data2$Conv)
Real_Time = as.numeric(RAFT_Real_Data2$Time) 
legend_text <- c(paste0("Conv.(%)", "   " ,"Time(mins)"),
                 paste0(round(Real_Conv[1], 2), "           ", round(Real_Time[1], 2)),
                 paste0(round(Real_Conv[2], 2), "           ", round(Real_Time[2], 2)),
                 paste0(round(Real_Conv[3], 2), "           ", round(Real_Time[3], 2)),
                 paste0(round(Real_Conv[4], 2), "           ", round(Real_Time[4], 2)))
RAFT_Real_Data <- read_excel("D:/Lloyd Yin/Desktop/Dissertation/RAFT_Real_Data.xlsx")
# View(RAFT_Real_Data)
layout(matrix(c(1,2,1,3),2,2))
Real_Time = names(RAFT_Real_Data)
plot(RAFT_Real_Data$`RT (mins)`,RAFT_Real_Data$`5h 38min`,type="l",col=5,lwd=2,
     ylim = c(min(RAFT_Real_Data[,-1]),max(RAFT_Real_Data[,-1])),
     main = "Real GPC Result")
for (i in Real_Time[c(-1,-5)]){
  lines(RAFT_Real_Data$`RT (mins)`,RAFT_Real_Data[[i]],col=which(i == Real_Time),lwd=2)
}
legend("topleft", legend = legend_text, col = c(NA, 2:5), lwd = c(NA, rep(2, 4)), text.font = c(2, rep(1, 4)), cex = 1.4)

plot(RAFT_Real_Data$`RT (mins)`,RAFT_Real_Data$`5h 38min`,type="l",col=5,lwd=2,
     xlim = c(14,20), ylim = c(0,290),main = "Zoom")
for (i in Real_Time[c(-1,-5)]){
  lines(RAFT_Real_Data$`RT (mins)`,RAFT_Real_Data[[i]],col=which(i == Real_Time),lwd=2)
}
legend("topleft", legend = legend_text, col = c(NA, 2:5), lwd = c(NA, rep(2, 4)), text.font = c(2, rep(1, 4)), cex = 0.8)

temp = which(RAFT_Real_Data$`RT (mins)`<=19 & RAFT_Real_Data$`RT (mins)`>14)
standareize = sum(RAFT_Real_Data$`5h 38min`[temp])
plot(RAFT_Real_Data$`RT (mins)`,RAFT_Real_Data$`5h 38min`/standareize,type="l",col=5,lwd=2,
     xlim = c(14,19), ylim = c(0,0.025),main = "Zoom (Standardize)")
for (i in Real_Time[c(-1,-5)]){
  standareize = sum(RAFT_Real_Data[[i]][temp])
  lines(RAFT_Real_Data$`RT (mins)`,RAFT_Real_Data[[i]]/standareize,col=which(i == Real_Time),lwd=2)
}
legend("topleft", legend = legend_text, col = c(NA, 2:5), lwd = c(NA, rep(2, 4)), text.font = c(2, rep(1, 4)), cex = 0.8)




library(readxl)
RAFT_Real_Data2 <- read_excel("D:/Lloyd Yin/Desktop/Dissertation/RAFT_Real_Data2.xlsx")
zero_row <- data.frame(Time = 0, Mn = 0, Mw = 0, PDI = NA, Conv = 0)
RAFT_Real_Data2=RAFT_Real_Data2[-4,]
RAFT_Real_Data2 = rbind(new_row, RAFT_Real_Data2)
# View(RAFT_Real_Data2)
# layout(1)
layout(matrix(1:4,2,2,byrow = TRUE))

plot(RAFT_Real_Data2$Time, RAFT_Real_Data2$Conv, xlab="Time (mins)", ylab = "Conv. (%)", main="Time v.s. Monomer conversion")
# logistic_model <- nls(Conv ~ SSlogis(Time, Asym, xmid, scal), data = RAFT_Real_Data2)
# predicted_time <- seq(0, 350, length.out = 200)
# predicted_conversion <- predict(logistic_model, newdata = data.frame(Time = predicted_time))
# lines(predicted_time, predicted_conversion, col = "blue", lwd = 1)
# lines(0:350,(1-exp(100*(exp(-10*10^(-5)*0:350)-1)))*100)
spline_model <- smooth.spline(RAFT_Real_Data2$Time,RAFT_Real_Data2$Conv)
predicted_time <- seq(0, 350, length.out = 200)
predicted_conversion <- predict(spline_model, predicted_time)$y
lines(predicted_time, predicted_conversion, col = "blue", lwd = 1)

plot(RAFT_Real_Data2$Conv, RAFT_Real_Data2$Mn, xlab="Conv. (%)", ylab = "Mn (Da)", main="Time v.s. Number-average molecular weight")
regression <- lm(Mn ~ 0 + Conv, data = RAFT_Real_Data2)
abline(regression, lwd = 1.5, col = "blue") 

plot(RAFT_Real_Data2$Conv, RAFT_Real_Data2$Mw, xlab="Conv. (%)", ylab = "Mw (Da)", main="Monomer conversion v.s. Weight-average molecular weight")
regression <- lm(Mw ~ 0 + Conv, data = RAFT_Real_Data2)
abline(regression, lwd = 1.5, col = "blue") 

plot(RAFT_Real_Data2$Conv, RAFT_Real_Data2$PDI, xlab="Conv. (%)", ylab = "PDI", main="Monomer conversion v.s. PDI",ylim=c(1,2))
neg_log_model <- nls(PDI ~ a - b * log(Conv + 1), data = RAFT_Real_Data2, start = list(a = 2, b = 0.5))
predicted_conv <- seq(0, 100, length.out = 200)
predicted_pdi <- predict(neg_log_model, newdata = data.frame(Conv = predicted_conv))
lines(predicted_conv, predicted_pdi, col = "blue", lwd = 1.5)