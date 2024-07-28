# RAFT Process ###################################################################
# I       -> 2R
# R + M   -> R
# R + CTA -> CTA-R
# CTA-R   -> R + CTA
# R + R   -> T
rm(list = ls())
options(max.print = 10)

# Base ###########################################################################
# M 1mol/L
#k_M = 3.78e-05
# M = 1039
# I = 3.16e-07
# CTA = 10^7

# 100:1:1
M = 10^9
I = 10^7
CTA = 10^7

P_lists = 1:1000/10000
P_critical = 0:8/10

Time = list()
# Na = 6.022*10^23

# Function #####################################################################
## Functions for C++ ############################################################
library(Rcpp)

cppFunction('
int cpp_sample_one(int n) {
  return floor(R::runif(0, 1) * n) + 1;
}')

cppFunction('
IntegerVector cpp_sample_two(int n) {
  if (n < 2) stop("n must be at least 2");

  IntegerVector result(2);
  result[0] = floor(R::runif(0, 1) * n) + 1;
  
  int second;
  do {
    second = floor(R::runif(0, 1) * n) + 1;
  } while (second == result[0]);

  result[1] = second;

  return result;
}')

cppFunction('
int cpp_weighted_sample_one(NumericVector prob) {
  int n = prob.size();
  double sum_prob = sum(prob);
  
  NumericVector normalized_prob = prob / sum_prob;
  double u = R::runif(0, 1);
  double cumulative_prob = 0.0;
  for (int i = 0; i < n; ++i) {
    cumulative_prob += normalized_prob[i];
    if (u < cumulative_prob) {
      return i + 1;
    }
  }
  return n; // Return the last element if not returned in the loop (this handles rounding errors)
}')

cppFunction('
NumericVector weighted_random_sample(NumericMatrix mat) {
  int nrow = mat.nrow();
  int ncol = mat.ncol();
  std::vector<double> weights;
  std::vector<int> row_indices;
  std::vector<int> col_indices;
  
  // Gather weights and their corresponding indices
  for (int i = 0; i < nrow; ++i) {
    for (int j = 0; j < ncol; ++j) {
      if (mat(i, j) > 0) {
        weights.push_back(mat(i, j));
        row_indices.push_back(i + 1); // R indices are 1-based
        col_indices.push_back(j + 1);
      }
    }
  }
  
  // Normalize weights
  double sum_weights = std::accumulate(weights.begin(), weights.end(), 0.0);
  std::vector<double> normalized_weights(weights.size());
  std::transform(weights.begin(), weights.end(), normalized_weights.begin(), 
                 [sum_weights](double w) { return w / sum_weights; });
  
  // Generate a random number
  double u = R::runif(0, 1);
  double cumulative_prob = 0.0;
  
  // Find the index based on the random number
  for (size_t k = 0; k < normalized_weights.size(); ++k) {
    cumulative_prob += normalized_weights[k];
    if (u < cumulative_prob) {
      return NumericVector::create(row_indices[k], col_indices[k]);
    }
  }
  
  return NumericVector::create(row_indices.back(), col_indices.back());
}')

## Functions for R #############################################################

# Free Radical System ##########################################################
# I       -> 2R
# R + M   -> R
# R + R   -> T
Start_time = Sys.time()

k_I       = 3.78*10^(-5)
k_p       = 5.20*10^(-6)
k_T       = 3.16*10^(-1)
k_CTAgent = 0
k_CTADE   = 0

# Current Stage
Monomer      = M
Initiator    = I
Radical      = NULL
CTAgent      = NULL
Intermediate = NULL
Termination  = rep(0,5000+1)  # Chain Length x: 0:5000
P_Current    = 0
Time_Current = 0

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
    Time        = rep(0,8001),
    Radical     = rep(0,8001),
    Termination = rep(0,8001),
    X_n         = rep(1,8001),
    X_w         = rep(1,8001),
    PDI         = NA,
    X_n2        = rep(1,8001),
    X_w2        = rep(1,8001),
    PDI2        = NA
  )
)
names(System$Critical) <- paste(c(1:8*10),"%")

# cat("\014")
for(i in 1:(length(P_critical)-1)){
  cat("\n",P_critical[i]*100,"% ~",P_critical[i+1]*100,"%\n")   # "\014",
  pb=txtProgressBar(min=1,max=length(P_lists),style=3)
  for (j in 1:length(P_lists)){
    P_threshold = P_critical[i]+P_lists[j]
    setTxtProgressBar(pb, j)
    while (P_Current < P_threshold){
      length_Radical = length(Radical)
      Temp.react = c(k_I*Initiator,
                     k_T*(length_Radical)^2)
      Temp.time = rexp(1,sum(Temp.react))
      Time_Current = Time_Current+Temp.time
      
      Propagation  = rpois(1, k_p*length_Radical*Monomer*Temp.time)
      Monomer = Monomer - Propagation
      Index = sample(length_Radical, Propagation, replace = TRUE)
      Index = factor(Index, levels = 1:length_Radical)
      Counts = table(Index)
      Index = 1:length_Radical
      Radical[Index] = Radical[Index] + as.numeric(Counts)
      P_Current = 1-Monomer/M
      
      Temp.react = Temp.react+min(Temp.react)*0.00001
      Temp.react = cumsum(Temp.react)/sum(Temp.react)
      Index.01 = runif(1)
      if (Index.01 < Temp.react[1]){
        Initiator = Initiator-1
        Radical = c(Radical,0,0)
      } else if (Index.01 < Temp.react[2] & length_Radical>=2){
        Index = cpp_sample_two(length_Radical) # sample(length_Radical,2,replace=FALSE)
        Length = sum(Radical[Index])
        Termination[Length+1] = Termination[Length+1] + 1
        Radical = Radical[-Index]
      }
      # print(length_Radical)
    }
    Index.State = 1000*(i-1)+j+1
    System$Series$Time       [Index.State] = Time_Current
    System$Series$Radical    [Index.State] = length(Radical)
    System$Series$Termination[Index.State] = sum(Termination)
    Data = Termination[-1]
    Counts = table(Radical)[-1]
    Index.Data = as.numeric(names(Counts))
    Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
    System$Series$X_n        [Index.State] = sum(Data*1:5000)/sum(Data)
    System$Series$X_w        [Index.State] = sum(Data*(1:5000)^2)/sum(Data*1:5000)
    System$Series$X_n2       [Index.State] = sum(Termination[-1]*1:5000)/sum(Termination[-1])
    System$Series$X_w2       [Index.State] = sum(Termination[-1]*(1:5000)^2)/sum(Termination[-1]*1:5000)
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
  Time$FRP[i] = difftime(End_time, Start_time, units = "secs")
}
System$Series$PDI  = System$Series$X_w /System$Series$X_n
System$Series$PDI2 = System$Series$X_w2/System$Series$X_n2

System.FRP = System
End_time = Sys.time()
Time$FRP.Total = difftime(End_time, Start_time, units = "mins")
save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_10^9_7.1.RData")

# RAFT.No System ###############################################################
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

# cat("\014")
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
        # # 选择并删除参与交换的自由基
        # Selected_CTA = sample(0:5000,CTAChange,prob=CTAgent,replace=TRUE)
        # All_elements = c(Radical,as.numeric(Selected_CTA))
        # Selected_CTA = table(Selected_CTA)
        # Index = as.numeric(names(Selected_CTA))
        # CTAgent[Index+1] = CTAgent[Index+1] - as.numeric(Selected_CTA)
        # 
        # # 四边形法则1, 构建反应顺序
        # Index.Propagation = sort(sample(1:(Propagation+CTAChange), CTAChange, replace=F))      # 转移次数位置   行
        # Index.CTAChange   =      sample(1:Length_Radical        , CTAChange, replace=T)       # 替换自由基坐标 列
        # Temp.Propagation  = c(rep(0,Length_Radical), Index.Propagation , rep(Propagation,Length_Radical))
        # Temp.CTAChange    = c(1:Length_Radical     , Index.CTAChange   , 1:Length_Radical)
        # 
        # # 处理顺序
        # Unique_elements <- unique(Temp.CTAChange)
        # S_mat     = sapply(Unique_elements, function(xxx) Temp.Propagation[which( Temp.CTAChange == xxx)[1:max(table( Temp.CTAChange))]]) # Survival_matrix
        # Index.mat = sapply(Unique_elements, function(xxx)                  which(Index.CTAChange == xxx)[1:max(table(Index.CTAChange))])
        # Index.mat[is.na(Index.mat)] = 0
        # S_mat[c(-1,-nrow(S_mat)), ] = S_mat[c(-1,-nrow(S_mat)), ] - Index.mat
        # S_mat                       = S_mat[-1,]-S_mat[-nrow(S_mat),]
        # S_Series = as.vector(na.omit(as.vector(t(S_mat)))) # 前Length_Radical个为原始自由基, 后面CTAChange个为Agent, 大小为存在时间(概率) #  Survival Series
        # 
        # # 增长
        # Index = sample(length(S_Series),Propagation,prob=S_Series,replace = TRUE)
        # Index = factor(Index, levels = 1:(Length_Radical+CTAChange))
        # Counts = table(Index)
        # Index = 1:(Length_Radical+CTAChange)
        # All_elements[Index] = All_elements[Index] + as.numeric(Counts)
        # 
        # # 返回自由基
        # Index = rev(c(1:Length_Radical,Index.CTAChange))
        # Index.Remain = Length_Radical+CTAChange+1-sapply(unique(Index), function(x) which(Index == x)[1])
        # Radical = All_elements[Index.Remain]
        # 
        # # 返回CTAngent
        # CTAgent_Back = table(All_elements[-Index.Remain])
        # Index = as.numeric(names(CTAgent_Back))
        # CTAgent[Index+1] = CTAgent[Index+1] + as.numeric(CTAgent_Back)
        
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
        CTAgent_Back = table(All_elements[1:(CTAChange)])
        Index = as.numeric(names(CTAgent_Back))+1
        CTAgent[Index] = CTAgent[Index] + as.numeric(CTAgent_Back)
      } else if (Propagation==0){
          for (k in 1:CTAChange){
            Index.1 = cpp_sample_one(Length_Radical)
            Index.2 = cpp_weighted_sample_one(CTAgent)
            CTAgent[Radical[Index.1]+1] = CTAgent[Radical[Index.1]+1] + 1
            Radical[Index.1] = Index.2 - 1
            CTAgent[Index.2] = CTAgent[Index.2] - 1
          }
      } else { # CTAChange==0
        Monomer = Monomer - Propagation
        P_Current = 1-Monomer/M
          for (k in 1:Propagation){
            Index.1 = cpp_sample_one(Length_Radical)
            Radical[Index.1] = Radical[Index.1] + 1
          }
      }
      
      Temp.react = Temp.react+min(Temp.react)*0.00001
      Temp.react = cumsum(Temp.react)/sum(Temp.react)
      Index.01 = runif(1)
      if (Index.01 < Temp.react[1]){
        Initiator = Initiator-1
        Radical = c(Radical,0,0)
      } else if (Index.01 < Temp.react[2] & Length_Radical>=2){
        Index = cpp_sample_two(Length_Radical) # sample(Length_Radical,2,replace=FALSE)
        Length = sum(Radical[Index])
        Termination[Length+1] = Termination[Length+1] + 1
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
  Time$RAFT.No[i] = difftime(End_time, Start_time, units = "secs")
}
System$Series$PDI  = System$Series$X_w /System$Series$X_n
System$Series$PDI2 = System$Series$X_w2/System$Series$X_n2

System.RAFT.No = System
End_time = Sys.time()
Time$RAFT.No.Total = difftime(End_time, Start_time, units = "mins")

save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_10^9_7.1.RData")


# sum(Data*1:500)+Monomer-10^9
# sum(CTAgent)-10^7
# Initiator+sum(Termination)+length(Radical)/2-10^7
# # RAFT System ################################################################
# I       -> 2R
# R + M   -> R
# R + CTA -> CTA-R
# CTA-R   -> R + CTA
# R + R   -> T
# Start_time = Sys.time()

k_I       = 3.78*10^(-5)
k_p       = 5.20*10^(-6)
k_T       = 3.16*10^(-1)
k_CTAgent = 5.00*10^(-4)
k_CTADE   = 5.00*10^(-4)

# Current Stage
Monomer      = M
Initiator    = I-1
Radical      = c(0,0)
CTAgent      = c(CTA,rep(0,5000))
Intermediate = NULL
Termination  = rep(0,5000+1) # Chain Length x: 0:5000
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
    X_n         = 0,
    X_w         = 0,
    PDI         = NA
  )
)
names(System$Critical) <- paste(c(1:8*10),"%")

# cat("\014")
for(i in 1:(length(P_critical)-1)){
  cat("\n",P_critical[i]*100,"% ~",P_critical[i+1]*100,"%\n")   # "\014",
  pb=txtProgressBar(min=1,max=length(P_lists),style=3)
  for (j in 1:length(P_lists)){
    P_threshold = P_critical[i]+P_lists[j]
    setTxtProgressBar(pb, j)
    while (P_Current < P_threshold){
      Temp.react = c(k_I*Initiator,
                     k_CTAgent*length(Radical)*sum(CTAgent),
                     k_CTADE*nrow(Intermediate),
                     k_T*(length(Radical))^2)
      Temp.time = rexp(1,sum(Temp.react))
      Time_Current = Time_Current+Temp.time

      Propagation  = rpois(1, k_p*length(Radical)*Monomer*Temp.time)
      Monomer = Monomer - Propagation
      for (k in 1:Propagation){
        Index.1 = cpp_sample_one(length(Radical))
        Radical[Index.1] = Radical[Index.1] + 1
      }
      P_Current = 1-Monomer/M
      print(length(Radical))
      Temp.react = Temp.react+min(Temp.react)*0.00001
      Temp.react = cumsum(Temp.react)/sum(Temp.react)
      Index.01 = runif(1)
      if (Index.01 < Temp.react[1]){
        Initiator = Initiator-1
        Radical = c(Radical,0,0)
      } else if (index < temp[2]){
        Index.1 = cpp_sample_one(length(Radical))
        Index.2 = cpp_sample_one(length(CTAgent))
        Intermediate = rbind(Intermediate,matrix(c(Radical[Index.1],CTAgent[Index.2]),1,2))
        Radical = Radical[-Index.1]
        CTAgent = CTAgent[-Index.2]
      } else if (index < temp[3]){
        Index.1 = cpp_sample_one(nrow(Intermediate))
        Index.2 = rbinom(1,1,0.5)
        Radical = c(Radical,Intermediate[Index.1,Index.2])
        CTAgent = c(CTAgent,Intermediate[Index.1,-Index.2])
        Intermediate = Intermediate[-Index.1, ,drop=FALSE]
      } else if (index < temp[4] & length(Radical)>=2){
        Index = cpp_sample_two(length(Radical)) # sample(length(Radical),2,replace=FALSE)
        Termination = c(Termination,Radical[Index[1]]+Radical[Index[2]])
        Radical = Radical[-Index]
      }
      # print(length(Radical))
    }
    System$Series$Time = c(System$Series$Time,Time_Current)
    System$Series$Radical = c(System$Series$Radical,length(Radical))
    System$Series$Termination = c(System$Series$Termination,sum(Termination))

    Data = Termination[-1]
    Counts = table(Radical)[-1]
    Index = as.numeric(names(Counts))
    Data[Index] = Data[Index] + as.numeric(Counts)
    Counts = CTAgent[-1]
    Index = which(Counts!=0)
    Counts = Counts[Counts != 0]
    Data[Index] = Data[Index] + as.numeric(Counts)
    System$Series$X_n = c(System$Series$X_n,sum(Data*1:5000)/sum(Data))
    System$Series$X_w = c(System$Series$X_w,sum(Data*(1:5000)^2)/sum(Data*1:5000))
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
}
System$Series$PDI = System$Series$X_w/System$Series$X_n

System.RAFT.No = System
End_time = Sys.time()
Time$RAFT.No = End_time - Start_time

x.1 = 0
Start_time = Sys.time()
for (i in 1:100000){
  x.1 = x.1 + sample(mat,1,prob=mat)
}
End_time = Sys.time()
End_time - Start_time


Start_time = Sys.time()
for (i in 1:100000){
  x.1 = x.1 + weighted_random_sample(mat)
}
End_time = Sys.time()
End_time - Start_time