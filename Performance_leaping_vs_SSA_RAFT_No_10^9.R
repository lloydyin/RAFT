# I       -> 2R
# R + M   -> R
# R + CTA -> CTA + R
# R + R   -> T
Start_time = Sys.time()

k_I       = 3.78*10^(-5)
k_p       = 5.20*10^(-6)
k_T       = 3.16*10^(-1)
k_CTAgent = 5.00*10^(-4)
k_CTADE   = 0

# Current Stage
Monomer      = M
Initiator    = I-1
Radical      = c(0,0)
CTAgent      = c(CTA,rep(0,500))
Intermediate = NULL
Termination  = rep(0,500+1) # Chain Length x: 0:5000
P_Current    = 0
Time_Current = rexp(1,k_I*Initiator)

System = list(
  k = list(
    I           = k_I,
    p           = k_p,
    Termination = k_T,
    CTAgent     = k_CTAgent,
    CTADE       = k_CTADE
  ),
  Initial  = list(
    Monomer   = M,
    Initiator = I,
    CTAgent   = NULL
  ),
  Critical = vector("list", 8),
  Series   = list(
    P           = 0:8000/10000,
    Time        = 0,
    Radical     = 0,
    Termination = 0,
    X_n         = 1,
    X_w         = 1,
    PDI         = NA,
    X_n2        = 1,
    X_w2        = 1,
    PDI2        = NA
  )
)
names(System$Critical) <- paste(c(1:8*10),"%")

for(i in 1:(length(P_critical)-1)){
  cat("\n",P_critical[i]*100,"% ~",P_critical[i+1]*100,"%\n")   # "\014",
  pb=txtProgressBar(min=1,max=length(P_lists),style=3)
  for (j in 1:length(P_lists)){
    P_threshold = P_critical[i]+P_lists[j]
    setTxtProgressBar(pb, j)
    while (P_Current < P_threshold){
      Length_Radical = length(Radical)
      Temp.react = c(k_I*Initiator,
                     k_T*(Length_Radical)^2)
      Temp.time = rexp(1,sum(Temp.react))
      Time_Current = Time_Current+Temp.time
      
      Propagation  = rpois(1, k_p*Length_Radical*Monomer*Temp.time)
      CTAChange    = rpois(1,k_CTAgent*Length_Radical*CTA*Temp.time)
      
      if (Propagation!=0 & CTAChange!=0){
        Monomer = Monomer - Propagation
        P_Current = 1-Monomer/M
        # 选择并删除参与交换的自由基
        Selected_CTA = sample(0:500,CTAChange,prob=CTAgent,replace=TRUE)
        All_elements = c(Radical,as.numeric(Selected_CTA))
        Selected_CTA = table(Selected_CTA)
        Index = as.numeric(names(Selected_CTA))+1
        CTAgent[Index] = CTAgent[Index] - as.numeric(Selected_CTA)
        
        # 增长
        Index = sample(length(All_elements),Propagation,replace = TRUE)
        Index = factor(Index, levels = 1:length(All_elements))
        Counts = table(Index)
        Index = 1:(Length_Radical+CTAChange)
        All_elements[Index] = All_elements[Index] + as.numeric(Counts)
        
        # 返回自由基
        Radical = All_elements[(CTAChange+1):(length(All_elements))]
        # 返回CTAngent
        CTAgent_Back = table(All_elements[1:CTAChange])
        Index = as.numeric(names(CTAgent_Back))+1
        CTAgent[Index] = CTAgent[Index] + as.numeric(CTAgent_Back)
      } else if (Propagation!=0){
        Monomer = Monomer - Propagation
        P_Current = 1-Monomer/M
        Index = sample(Length_Radical, Propagation, replace = TRUE)
        Index = factor(Index, levels = 1:Length_Radical)
        Counts = table(Index)
        Index = 1:Length_Radical
        Radical[Index] = Radical[Index] + as.numeric(Counts)
      } else if (CTAChange!=0){
        Selected_CTA = sample(0:500,CTAChange,prob=CTAgent,replace=TRUE)
        All_elements = c(Radical,Selected_CTA)
        Selected_CTA = table(Selected_CTA)
        Index = as.numeric(names(Selected_CTA))+1
        CTAgent[Index] = CTAgent[Index] - as.numeric(Selected_CTA)
        
        Radical = All_elements[(CTAChange+1):(length(All_elements))]
        CTAgent_Back = table(All_elements[1:CTAChange])
        Index = as.numeric(names(CTAgent_Back))+1
        CTAgent[Index] = CTAgent[Index] + as.numeric(CTAgent_Back)
        
      }
      
      Temp.react = Temp.react+min(Temp.react)*0.00001
      Temp.react = cumsum(Temp.react)/sum(Temp.react)
      Index.01 = runif(1)
      if (Index.01 < Temp.react[1]){
        Initiator = Initiator-1
        Radical = c(Radical,0,0)
      } else if (Index.01 < Temp.react[2] & Length_Radical>=2){
        Index = cpp_sample_two(Length_Radical) # sample(Length_Radical,2,replace=FALSE)
        Length = sum(Radical[Index])+1
        Termination[Length] = Termination[Length] + 1
        Radical = Radical[-Index]
      }
      # print(Length_Radical)
    }
    Index.State = 1000*(i-1)+j+1
    System$Series$Time[Index.State] = Time_Current
    System$Series$Radical[Index.State] = length(Radical)
    System$Series$Termination[Index.State] = sum(Termination)
    Data = Termination[-1]
    Counts = table(Radical)[-1]
    Index.Data = as.numeric(names(Counts))
    Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
    Counts = CTAgent[-1]
    Index.Data = which(Counts!=0)
    Counts = Counts[Counts != 0]
    Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
    System$Series$X_n [Index.State] = sum(Data*1:500)/sum(Data)
    System$Series$X_w [Index.State] = sum(Data*(1:500)^2)/sum(Data*1:500)
    System$Series$X_n2[Index.State] = sum(Termination[-1]*1:500)/sum(Termination[-1])
    System$Series$X_w2[Index.State] = sum(Termination[-1]*(1:500)^2)/sum(Termination[-1]*1:500)
  }
  System$Critical[[i]] = list(
    P=P_Current,
    Monomer=Monomer,
    Initiator=Initiator,
    Radical=Radical,
    CTAgent=CTAgent,
    Intermediate=Intermediate,
    Termination=Termination,
    Data=Data)
  End_time = Sys.time()
  Time$RAFT.No.Leaping[i] = difftime(End_time, Start_time, units = "secs")
}
System$Series$PDI  = System$Series$X_w /System$Series$X_n
System$Series$PDI2 = System$Series$X_w2/System$Series$X_n2

System.RAFT.No = System
End_time = Sys.time()
Time$RAFT.No.Leaping.Total= difftime(End_time, Start_time, units = "mins")



Start_time = Sys.time()


# Current Stage
Monomer      = M
Initiator    = I-1
Radical      = c(0,0)
CTAgent      = c(CTA,rep(0,500))
Intermediate = NULL
Termination  = rep(0,500+1) # Chain Length x: 0:5000
P_Current    = 0
Time_Current = rexp(1,k_I*Initiator)

System = list(
  k = list(
    I           = k_I,
    p           = k_p,
    Termination = k_T,
    CTAgent     = k_CTAgent,
    CTADE       = k_CTADE
  ),
  Initial  = list(
    Monomer   = M,
    Initiator = I,
    CTAgent   = NULL
  ),
  Critical = vector("list", 8),
  Series   = list(
    P           = 0:8000/10000,
    Time        = 0,
    Radical     = 0,
    Termination = 0,
    X_n         = 1,
    X_w         = 1,
    PDI         = NA,
    X_n2        = 1,
    X_w2        = 1,
    PDI2        = NA
  )
)
names(System$Critical) <- paste(c(1:8*10),"%")

for(i in 1:(length(P_critical)-1)){
  cat("\n",P_critical[i]*100,"% ~",P_critical[i+1]*100,"%\n")   # "\014",
  pb=txtProgressBar(min=1,max=length(P_lists),style=3)
  for (j in 1:length(P_lists)){
    P_threshold = P_critical[i]+P_lists[j]
    setTxtProgressBar(pb, j)
    while (P_Current < P_threshold){
      Length_Radical = length(Radical)
      Temp.react = c(k_I*Initiator,
                     k_p*Length_Radical*Monomer,
                     k_CTAgent*Length_Radical*CTA,
                     k_T*(Length_Radical)^2)
      Temp.time = rexp(1,sum(Temp.react))
      Time_Current = Time_Current+Temp.time
      
      Temp.react = Temp.react+min(Temp.react)*0.00001
      Temp.react = cumsum(Temp.react)/sum(Temp.react)
      Index.01 = runif(1)
      if (Index.01 < Temp.react[1]){
        Initiator = Initiator-1
        Radical = c(Radical,0,0)
      } else if (Index.01 < Temp.react[2]){
        Monomer = Monomer-1
        Index = cpp_sample_one(Length_Radical)
        Radical[Index] = Radical[Index]+1
        P_Current = 1-Monomer/M
      } else if (Index.01 < Temp.react[3]){
        Index.1 = cpp_sample_one(Length_Radical)
        Index.2 = cpp_weighted_sample_one(CTAgent)
        CTAgent[Index.2] = CTAgent[Index.2] - 1
        CTAgent[Radical[Index.1]+1] = CTAgent[Radical[Index.1]+1] + 1
        Radical[Index.1] = Index.2 - 1
      } else if (Index.01 < Temp.react[4] & Length_Radical>=2){
        Index = cpp_sample_two(Length_Radical) # sample(Length_Radical,2,replace=FALSE)
        Length = sum(Radical[Index])+1
        Termination[Length] = Termination[Length] + 1
        Radical = Radical[-Index]
      }
      # print(Length_Radical)
    }
    Index.State = 1000*(i-1)+j+1
    System$Series$Time[Index.State] = Time_Current
    System$Series$Radical[Index.State] = length(Radical)
    System$Series$Termination[Index.State] = sum(Termination)
    Data = Termination[-1]
    Counts = table(Radical)[-1]
    Index.Data = as.numeric(names(Counts))
    Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
    Counts = CTAgent[-1]
    Index.Data = which(Counts!=0)
    Counts = Counts[Counts != 0]
    Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
    System$Series$X_n [Index.State] = sum(Data*1:500)/sum(Data)
    System$Series$X_w [Index.State] = sum(Data*(1:500)^2)/sum(Data*1:500)
    System$Series$X_n2[Index.State] = sum(Termination[-1]*1:500)/sum(Termination[-1])
    System$Series$X_w2[Index.State] = sum(Termination[-1]*(1:500)^2)/sum(Termination[-1]*1:500)
  }
  System$Critical[[i]] = list(
    P=P_Current,
    Monomer=Monomer,
    Initiator=Initiator,
    Radical=Radical,
    CTAgent=CTAgent,
    Intermediate=Intermediate,
    Termination=Termination,
    Data=Data)
  End_time = Sys.time()
  Time$RAFT.No.SSA[i] = difftime(End_time, Start_time, units = "secs")
}
System$Series$PDI  = System$Series$X_w /System$Series$X_n
System$Series$PDI2 = System$Series$X_w2/System$Series$X_n2

System.RAFT.No = System
End_time = Sys.time()
Time$RAFT.No.SSA.Total = difftime(End_time, Start_time, units = "mins")

sum(Radical)+Monomer+sum(Termination*0:500)+sum(CTAgent*0:500)-10^9
sum(CTAgent)
Initiator+length(Radical)/2+sum(Termination)
