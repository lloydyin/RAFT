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
M = 10^10
I = 10^8
CTA = 10^8

P_lists = 1:1000/10000
P_critical = 0:8/10

Time = list()
# Na = 6.022*10^23

# Function #####################################################################
# 均匀选一个         cpp_sample_one(v)
# 均匀选两个         cpp_sample_two(v)
# 按照向量权重选一个 cpp_weighted_sample_one(v)
# 按矩阵权重选一个   samplePoint(m)              # 1*2  matrix
# 按照矩阵权重选多个 sampleMultiplePoints(m,n)   #
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
  NumericMatrix samplePoint(NumericMatrix mat) {
    int rows = mat.nrow();
    int cols = mat.ncol();
    int total_size = rows * cols;

    // Flatten the matrix and calculate cumulative weights
    NumericVector mat_vector(total_size);
    NumericVector cumulative_weights(total_size);
    double sum_weights = 0;

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        int index = i + j * rows;
        mat_vector[index] = mat(i, j);
        sum_weights += mat(i, j);
        cumulative_weights[index] = sum_weights;
      }
    }

    // Generate a random number in the range of cumulative weights
    double rand_value = R::runif(0, sum_weights);

    // Perform a binary search to find the corresponding index
    int left = 0, right = total_size - 1, mid;
    while (left < right) {
      mid = (left + right) / 2;
      if (cumulative_weights[mid] < rand_value)
        left = mid + 1;
      else
        right = mid;
    }

    // Convert flat index back to matrix coordinates
    int selected_row = left % rows + 1;
    int selected_col = left / rows + 1;

    NumericMatrix result(1, 2);
    result(0, 0) = selected_row;
    result(0, 1) = selected_col;

    return result;
  }
')

cppFunction('
  NumericMatrix sampleMultiplePoints(NumericMatrix mat, int n_samples) {
    int rows = mat.nrow();
    int cols = mat.ncol();
    int total_size = rows * cols;
    
    // Flatten the matrix and calculate cumulative weights
    NumericVector mat_vector(total_size);
    NumericVector cumulative_weights(total_size);
    double sum_weights = 0;

    for (int i = 0; i < rows; ++i) {
      for (int j = 0; j < cols; ++j) {
        int index = i + j * rows;
        mat_vector[index] = mat(i, j);
        sum_weights += mat(i, j);
        cumulative_weights[index] = sum_weights;
      }
    }

    NumericMatrix selected_points(n_samples, 2);

    for (int sample = 0; sample < n_samples; ++sample) {
      // If sum_weights is zero, stop sampling
      if (sum_weights <= 0) {
        NumericMatrix result(sample, 2);
        for (int i = 0; i < sample; ++i) {
          result(i, 0) = selected_points(i, 0);
          result(i, 1) = selected_points(i, 1);
        }
        return result;
      }

      // Generate a random number in the range of cumulative weights
      double rand_value = R::runif(0, sum_weights);

      // Perform a binary search to find the corresponding index
      int left = 0, right = total_size - 1, mid;
      while (left < right) {
        mid = (left + right) / 2;
        if (cumulative_weights[mid] < rand_value)
          left = mid + 1;
        else
          right = mid;
      }

      // Convert flat index back to matrix coordinates
      int selected_row = left % rows;
      int selected_col = left / rows;

      selected_points(sample, 0) = selected_row + 1; // Convert to 1-based index
      selected_points(sample, 1) = selected_col + 1; // Convert to 1-based index

      // Update weights and cumulative weights
      mat_vector[left] -= 1.0; 
            if (mat_vector[left] < 0) {
                mat_vector[left] = 0; // Ensure weight does not go below 0
            }
            
            sum_weights = 0;
            for (int i = 0; i < total_size; ++i) {
                sum_weights += mat_vector[i];
                cumulative_weights[i] = sum_weights;
            }
            }

return selected_points;
}')


## Functions for R #############################################################
Leaping_Probablity = function(x, y) {
  total_length = x + y
  
  if (x < y) {
    # Ascend from 1 to x, then stay at x for (y - x) times, then descend back to 1
    sequence <- c(1:x, rep(x, y - x), x:1)
  } else {
    # Calculate the ascend and descend lengths
    ascend_length = (total_length + 1) %/% 2
    descend_length = total_length %/% 2
    
    # Generate the sequence
    sequence = c(1:ascend_length, descend_length :1)
  }
  
  return(sequence)
}

# Free Radical System ##########################################################
# I       -> 2R
# R + M   -> R
# R + R   -> T
Start_time = Sys.time()

k_I       = 3.78*10^(-5)
k_p       = 5.20*10^(-7)
k_T       = 3.16*10^(-2)
k_CTAgent = 0
k_CTADE   = 0

# Current Stage
Monomer      = M
Initiator    = I
Radical      = NULL
CTAgent      = NULL
Intermediate = NULL
Termination  = rep(0,5000+1) # Chain Length x: 0:5000
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
      Temp.react = c(k_I*Initiator,
                     k_T*(length(Radical))^2)
      Temp.time = rexp(1,sum(Temp.react))
      Time_Current = Time_Current+Temp.time
      
      Propagation  = rpois(1, k_p*length(Radical)*Monomer*Temp.time)
      Monomer = Monomer - Propagation
      Index = sample(length(Radical), Propagation, replace = TRUE)
      Index = factor(Index, levels = 1:length(Radical))
      Counts = table(Index)
      Index = 1:length(Radical)
      Radical[Index] = Radical[Index] + as.numeric(Counts)
      P_Current = 1-Monomer/M
      
      Temp.react = Temp.react+min(Temp.react)*0.00001
      Temp.react = cumsum(Temp.react)/sum(Temp.react)
      Index.01 = runif(1)
      if (Index.01 < Temp.react[1]){
        Initiator = Initiator-1
        Radical = c(Radical,0,0)
      } else if (Index.01 < Temp.react[2] & length(Radical)>=2){
        Index = cpp_sample_two(length(Radical)) # sample(length(Radical),2,replace=FALSE)
        Length = sum(Radical[Index])+1
        Termination[Length] = Termination[Length] + 1
        Radical = Radical[-Index]
      }
      # print(length(Radical))
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
save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_10^10_7.05_2.RData")

# RAFT.No System ###############################################################
# I       -> 2R
# R + M   -> R
# R + CTA -> CTA + R
# R + R   -> T
Start_time = Sys.time()

k_I       = 3.78*10^(-5)
k_p       = 5.20*10^(-7)
k_T       = 3.16*10^(-2)
k_CTAgent = 2.00*10^(-5)
k_CTADE   = 0

# Current Stage
Monomer      = M
Initiator    = I-1
Radical      = c(0,0)
CTAgent      = c(CTA,rep(0,500))
Intermediate = NULL
Termination  = rep(0,500+1) # Chain Length x: 0:500
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
        Monomer = Monomer - Propagation
        P_Current = 1-Monomer/M
        # Select and delete exchange Radical
        Selected_CTA = sample(0:500,CTAChange,prob=CTAgent,replace=TRUE)
        All_elements = c(Radical,as.numeric(Selected_CTA))
        Selected_CTA = table(Selected_CTA)
        Index = as.numeric(names(Selected_CTA))+1
        CTAgent[Index] = CTAgent[Index] - as.numeric(Selected_CTA)
        
        # Propagation
        Index = sample(length(All_elements),Propagation,replace = TRUE)
        Index = factor(Index, levels = 1:length(All_elements))
        Counts = table(Index)
        Index = 1:(Length_Radical+CTAChange)
        All_elements[Index] = All_elements[Index] + as.numeric(Counts)
        
        # Return Free Radical
        Radical = All_elements[(CTAChange+1):(length(All_elements))]
        # Return CTAngent
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

save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_10^10_7.05_2.RData")

# sum(Data*1:500)+Monomer-10^10
# sum(CTAgent)-10^8
# Initiator+sum(Termination)+length(Radical)/2-10^8
# RAFT System ##################################################################
# I       -> 2R
# R + M   -> R
# R + CTA -> CTA + R
# R + R   -> T
Start_time = Sys.time()

k_I       = 3.78*10^(-5)
k_p       = 5.20*10^(-7)
k_T       = 3.16*10^(-2)
k_CTAgent = 2.00*10^(-5)
k_CTADE   = 2.00*10^(-3)

# Current Stage
Monomer      = M
Initiator    = I-1
Radical      = c(0,0)
CTAgent      = c(CTA,rep(0,500))             # 0:500
Intermediate = matrix(rep(0,501^2),501,501)  # 0:500,0:500
Termination  = rep(0,1000+1)                  # Chain Length x: 0:500
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

# 3592.528 77497.493 27027.021 19819.858     7.110
  
Leaping_condition = FALSE
Exit_outer_loop = FALSE
for(i in 1:(length(P_critical)-1)){
  cat("\n",P_critical[i]*100,"% ~",P_critical[i+1]*100,"%\n")   # "\014",
  pb=txtProgressBar(min=1,max=length(P_lists),style=3)
  if (Leaping_condition) {
    for (j in 1:length(P_lists)){
      P_threshold = P_critical[i]+P_lists[j]
      setTxtProgressBar(pb, j)
      while (P_Current < P_threshold){
        # Leaping
        Length_Radical = length(Radical)
        Temp.react = c(k_I*Initiator,
                       k_T*(Length_Radical)^2)
        Temp.time = rexp(1,sum(Temp.react))
        Time_Current = Time_Current+Temp.time
        
        Propagation  = rpois(1, k_p       *Length_Radical *Monomer     *Temp.time)
        CTAChange    = rpois(1, k_CTAgent *Length_Radical *sum(CTAgent)*Temp.time)
        # print(c(Length_Radical,Propagation,CTAChange))
        cat("\r",P_Current)
        # cat("\r",Length_Radical)
        if (Propagation!=0 & CTAChange!=0){
          Monomer = Monomer - Propagation
          P_Current = 1-Monomer/M
          # Select and delete exchange Radical
          Selected_Intermediate = sampleMultiplePoints(Intermediate,CTAChange)      # a n*2 matrix, every colmun is an index 
          # Intermediate = Selected_Intermediate$Updated_Matrix
          Selected_Intermediate = Selected_Intermediate$Selected_Index
          Index = rbinom(CTAChange,1,0.5) + 1                                       # 1 is return to Radical, 2 is return to CTAgent
          CTA_Return = table(Selected_Intermediate[matrix(c(1:CTAChange,Index),,2)])
          # return to CTAgent
          
          # 
          All_elements = c(Radical,Selected_Intermediate[matrix(c(1:CTAChange,3-Index),,2)]-1)
          Index = sample(length(All_elements),Propagation,
                         prob = Leaping_Probablity(Length_Radical,CTAChange),
                         replace = TRUE)
          Index = factor(Index, levels = 1:length(All_elements))
          Counts = table(Index)
          Index = as.numeric(names(Counts)) 
          All_elements[Index] = All_elements[Index] + as.numeric(Counts)
          
          # Return Radical
          Radical = All_elements[(CTAChange+1):(length(All_elements))]
          
          # Deal CTAgent
          Selected_CTA = sample(1:501,CTAChange,prob=CTAgent,replace=TRUE)
          # Delete CTA
          Selected_CTA_table = table(Selected_CTA)
          Index = as.numeric(names(Selected_CTA_table))
          CTAgent[Index] = CTAgent[Index] - as.numeric(Selected_CTA_table)
          # Return CTA
          Index = as.numeric(names(CTA_Return)) # (1:501)
          CTAgent[Index] = CTAgent[Index] + as.numeric(CTA_Return)
          
          # Return Intermediate
          Intermediate_return = matrix(c(All_elements[1:(CTAChange)]+1,Selected_CTA),,2)
          incrementMatrix(Intermediate, Intermediate_return)
        } else if (CTAChange == 0){
          Monomer = Monomer - Propagation
          P_Current = 1-Monomer/M
          Index = sample(length(Radical), Propagation, replace = TRUE)
          Index = factor(Index, levels = 1:length(Radical))
          Counts = table(Index)
          Index = 1:length(Radical)
          Radical[Index] = Radical[Index] + as.numeric(Counts)
        } else if (Propagation == 0){
          Selected_Intermediate = sampleMultiplePoints(Intermediate,CTAChange)      # a n*2 matrix, every colmun is an index 
          Selected_Intermediate = Selected_Intermediate$Selected_Index
          Index = rbinom(CTAChange,1,0.5) + 1                                       # 1 is return to Radical, 2 is return to CTAgent
          CTA_Return = table(Selected_Intermediate[matrix(c(1:CTAChange,Index),,2)])
          All_elements = c(Radical,Selected_Intermediate[matrix(c(1:CTAChange,3-Index),,2)]-1)
          Radical = All_elements[(CTAChange+1):(length(All_elements))]
          Selected_CTA = sample(1:501,CTAChange,prob=CTAgent,replace=TRUE)
          Selected_CTA_table = table(Selected_CTA)
          Index = as.numeric(names(Selected_CTA_table))
          CTAgent[Index] = CTAgent[Index] - as.numeric(Selected_CTA_table)
          Index = as.numeric(names(CTA_Return)) # (1:501)
          CTAgent[Index] = CTAgent[Index] + as.numeric(CTA_Return)
          Intermediate_return = matrix(c(All_elements[1:(CTAChange)]+1,Selected_CTA),,2)
          incrementMatrix(Intermediate, Intermediate_return)
        }
        Temp.react = Temp.react+min(Temp.react)*0.00001
        Temp.react = cumsum(Temp.react)/sum(Temp.react)
        index = runif(1)
        
        if (index < Temp.react[1]){
          Index.1 = cpp_sample_one(Length_Radical)
          Index.2 = sample(0:500,1,prob=CTAgent) + 1 # Chain Length 
          Choosen_Radical =  Radical[Index.1] + 1
          Intermediate[Choosen_Radical,Index.2] = Intermediate[Choosen_Radical,Index.2] + 1
          Radical = Radical[-Index.1]
          CTAgent[Index.2] = CTAgent[Index.2] - 1
          
          Index.1 = cpp_sample_one(Length_Radical)
          Index.2 = sample(0:500,1,prob=CTAgent) + 1 # Chain Length 
          Choosen_Radical =  Radical[Index.1] + 1
          Intermediate[Choosen_Radical,Index.2] = Intermediate[Choosen_Radical,Index.2] + 1
          Radical = Radical[-Index.1]
          CTAgent[Index.2] = CTAgent[Index.2] - 1
          
          Initiator = Initiator-1
          Radical = c(Radical,0,0)
        } else if (index < Temp.react[2] & Length_Radical>=2){
          Index.1 = samplePoint(Intermediate)
          Intermediate[Index.1] = Intermediate[Index.1] - 1
          Index.2 = rbinom(1,1,0.5)+1
          Radical = c(Radical,Index.1[Index.2] - 1)
          CTAgent[Index.1[3-Index.2]] = CTAgent[Index.1[3-Index.2]] + 1
          
          Index.1 = samplePoint(Intermediate)
          Intermediate[Index.1] = Intermediate[Index.1] - 1
          Index.2 = rbinom(1,1,0.5)+1
          Radical = c(Radical,Index.1[Index.2] - 1)
          CTAgent[Index.1[3-Index.2]] = CTAgent[Index.1[3-Index.2]] + 1
          
          Index = cpp_sample_two(Length_Radical) # sample(Length_Radical,2,replace=FALSE)
          Length = sum(Radical[Index]) + 1
          Termination[Length+1] = Termination[Length+1] + 1
          Radical = Radical[-Index]
        }
        
        for (k in 1:10){
          Length_Radical = length(Radical)
          Temp.react = c(k_I*Initiator,
                         k_CTAgent*Length_Radical*sum(CTAgent),
                         k_CTADE  *sum(Intermediate),
                         k_T*(Length_Radical)^2)
          Temp.time = rexp(1,sum(Temp.react))
          Time_Current = Time_Current+Temp.time
          
          Propagation  = rpois(1, k_p*Length_Radical*Monomer*Temp.time)
          Monomer = Monomer - Propagation
          Index = sample(Length_Radical, Propagation, replace = TRUE)
          Index = factor(Index, levels = 1:Length_Radical)
          Counts = table(Index)
          Index = 1:Length_Radical
          Radical[Index] = Radical[Index] + as.numeric(Counts)
          P_Current = 1-Monomer/M
          
          Temp.react = Temp.react+min(Temp.react)*0.00001
          Temp.react = cumsum(Temp.react)/sum(Temp.react)
          index = runif(1)
          
          if (index < Temp.react[1]){
            Initiator = Initiator-1
            Radical = c(Radical,0,0)
          } else if (index < Temp.react[2]){
            Index.1 = cpp_sample_one(Length_Radical)
            Index.2 = sample(0:500,1,prob=CTAgent) + 1 # Chain Length 
            Choosen_Radical =  Radical[Index.1] + 1
            Intermediate[Choosen_Radical,Index.2] = Intermediate[Choosen_Radical,Index.2] + 1
            Radical = Radical[-Index.1]
            CTAgent[Index.2] = CTAgent[Index.2] - 1
          } else if (index < Temp.react[3]){
            Index.1 = samplePoint(Intermediate)
            Intermediate[Index.1] = Intermediate[Index.1] - 1
            Index.2 = rbinom(1,1,0.5)+1
            Radical = c(Radical,Index.1[Index.2] - 1)
            CTAgent[Index.1[3-Index.2]] = CTAgent[Index.1[3-Index.2]] + 1
          } else if (index < Temp.react[4] & Length_Radical>=2){
            Index = cpp_sample_two(Length_Radical) # sample(Length_Radical,2,replace=FALSE)
            Length = sum(Radical[Index]) + 1
            Termination[Length+1] = Termination[Length+1] + 1
            Radical = Radical[-Index]
          }
        }
      }
      Index.State = 1000*(i-1)+j+1
      System$Series$Time[Index.State] = Time_Current
      System$Series$Radical[Index.State] = length(Radical)
      System$Series$Termination[Index.State] = sum(Termination)
      Data = c(Termination[-c(1,2)],0)
      Counts = table(Radical)[-1]
      Index.Data = as.numeric(names(Counts))
      Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
      Counts = CTAgent[-1]
      Index.Data = which(Counts!=0)
      Counts = Counts[Counts != 0]
      Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
      Intermediate.table = sapply(2:1001, function(d) {sum(Intermediate[row(Intermediate) + col(Intermediate) == d + 1])}) # ingore the first term length=(0,0)
      Data = Data + Intermediate.table
      System$Series$X_n [Index.State] = sum(Data*1:1000)/sum(Data)
      System$Series$X_w [Index.State] = sum(Data*(1:1000)^2)/sum(Data*1:1000)
      System$Series$X_n2[Index.State] = sum(Termination[-1]*1:1000)/sum(Termination[-1])
      System$Series$X_w2[Index.State] = sum(Termination[-1]*(1:1000)^2)/sum(Termination[-1]*1:1000)
    }
  } else {
    for (j in 1:length(P_lists)){
      P_threshold = P_critical[i]+P_lists[j]
      setTxtProgressBar(pb, j)
      
      while (P_Current < P_threshold){
        Length_Radical = length(Radical)
        # print(Length_Radical)
        A.CTAgent = k_CTAgent*Length_Radical*sum(CTAgent)
        A.CTADE   = k_CTADE  *sum(Intermediate)
        if (A.CTAgent >= A.CTADE || Length_Radical <= 10) {
          # normal SSA
          Temp.react = c(k_I*Initiator,
                         k_p*Length_Radical*Monomer,
                         A.CTAgent,
                         A.CTADE,
                         k_T*(Length_Radical)^2)
          Temp.time = rexp(1,sum(Temp.react))
          Time_Current = Time_Current+Temp.time

          Temp.react = Temp.react+min(Temp.react)*0.00001
          Temp.react = cumsum(Temp.react)/sum(Temp.react)
          index = runif(1)

          if (index < Temp.react[1]){
            Initiator = Initiator-1
            Radical = c(Radical,0,0)
          } else if (index < Temp.react[2]){
            Monomer = Monomer-1
            Index = cpp_sample_one(length(Radical))
            Radical[Index] = Radical[Index] + 1
            P_Current = 1-Monomer/M
          } else if (index < Temp.react[3]){
            Index.1 = cpp_sample_one(length(Radical))
            Index.2 = sample(0:500,1,prob=CTAgent) + 1 # Chain Length 
            Choosen_Radical =  Radical[Index.1] + 1
            Intermediate[Choosen_Radical,Index.2] = Intermediate[Choosen_Radical,Index.2] + 1
            Radical = Radical[-Index.1]
            CTAgent[Index.2] = CTAgent[Index.2] - 1
          } else if (index < Temp.react[4]){
            Index.1 = samplePoint(Intermediate)
            Intermediate[Index.1] = Intermediate[Index.1] - 1
            Index.2 = rbinom(1,1,0.5)+1
            Radical = c(Radical,Index.1[Index.2] - 1)
            CTAgent[Index.1[3-Index.2]] = CTAgent[Index.1[3-Index.2]] + 1
          } else if (index < Temp.react[5] & length(Radical)>=2){
            Index = cpp_sample_two(length(Radical)) # sample(length(Radical),2,replace=FALSE)
            Length = sum(Radical[Index]) + 1
            Termination[Length+1] = Termination[Length+1] + 1
            Radical = Radical[-Index]
          }
        } else {
          # get leaping
          Leaping_condition = TRUE
          Exit_outer_loop = TRUE
          P_Relaxation = P_Current
          Time_Relaxation = Time_Current
          save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical_Point.RData")
          break
        }
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
      Intermediate.table = sapply(2:1001, function(d) {sum(Intermediate[row(Intermediate) + col(Intermediate) == d + 1])}) # ingore the first term length=(0,0)
      Data = Data + Intermediate.table
      System$Series$X_n [Index.State] = sum(Data*1:1000)/sum(Data)
      System$Series$X_w [Index.State] = sum(Data*(1:1000)^2)/sum(Data*1:1000)
      System$Series$X_n2[Index.State] = sum(Termination[-1]*1:1000)/sum(Termination[-1])
      System$Series$X_w2[Index.State] = sum(Termination[-1]*(1:1000)^2)/sum(Termination[-1]*1:1000)
      if (Exit_outer_loop) {
        break
      }
    }
    if (Leaping_condition) {
      for (j in j:length(P_lists)) {
        P_threshold = P_critical[i]+P_lists[j]
        setTxtProgressBar(pb, j)
        while (P_Current < P_threshold){
          # Leaping
          Length_Radical = length(Radical)
          Temp.react = c(k_I*Initiator,
                         k_T*(Length_Radical)^2)
          Temp.time = rexp(1,sum(Temp.react))
          Time_Current = Time_Current+Temp.time
          
          Propagation  = rpois(1, k_p       *Length_Radical *Monomer     *Temp.time)
          CTAChange    = rpois(1, k_CTAgent *Length_Radical *sum(CTAgent)*Temp.time)
          # print(c(Length_Radical,Propagation,CTAChange))
          cat("\r",P_Current)
          # cat("\r",Length_Radical)
          if (Propagation!=0 & CTAChange!=0){
            Monomer = Monomer - Propagation
            P_Current = 1-Monomer/M
            # 选择并删除参与交换的自由基
            Selected_Intermediate = sampleMultiplePoints(Intermediate,CTAChange)      # a n*2 matrix, every colmun is an index 
            # Intermediate = Selected_Intermediate$Updated_Matrix
            Selected_Intermediate = Selected_Intermediate$Selected_Index
            Index = rbinom(CTAChange,1,0.5) + 1                                       # 1 is return to Radical, 2 is return to CTAgent
            CTA_Return = table(Selected_Intermediate[matrix(c(1:CTAChange,Index),,2)])
            # return to CTAgent
            
            # 
            All_elements = c(Radical,Selected_Intermediate[matrix(c(1:CTAChange,3-Index),,2)]-1)
            Index = sample(length(All_elements),Propagation,
                           prob = Leaping_Probablity(Length_Radical,CTAChange),
                           replace = TRUE)
            Index = factor(Index, levels = 1:length(All_elements))
            Counts = table(Index)
            Index = as.numeric(names(Counts)) 
            All_elements[Index] = All_elements[Index] + as.numeric(Counts)
            
            # Return Radical
            Radical = All_elements[(CTAChange+1):(length(All_elements))]
            
            # Deal CTAgent
            Selected_CTA = sample(1:501,CTAChange,prob=CTAgent,replace=TRUE)
            # Delete CTA
            Selected_CTA_table = table(Selected_CTA)
            Index = as.numeric(names(Selected_CTA_table))
            CTAgent[Index] = CTAgent[Index] - as.numeric(Selected_CTA_table)
            # Return CTA
            Index = as.numeric(names(CTA_Return)) # (1:501)
            CTAgent[Index] = CTAgent[Index] + as.numeric(CTA_Return)
            
            # Return Intermediate
            Intermediate_return = matrix(c(All_elements[1:(CTAChange)]+1,Selected_CTA),,2)
            incrementMatrix(Intermediate, Intermediate_return)
          } else if (CTAChange == 0){
            Monomer = Monomer - Propagation
            P_Current = 1-Monomer/M
            Index = sample(length(Radical), Propagation, replace = TRUE)
            Index = factor(Index, levels = 1:length(Radical))
            Counts = table(Index)
            Index = 1:length(Radical)
            Radical[Index] = Radical[Index] + as.numeric(Counts)
          } else if (Propagation == 0){
            Selected_Intermediate = sampleMultiplePoints(Intermediate,CTAChange)      # a n*2 matrix, every colmun is an index 
            Selected_Intermediate = Selected_Intermediate$Selected_Index
            Index = rbinom(CTAChange,1,0.5) + 1                                       # 1 is return to Radical, 2 is return to CTAgent
            CTA_Return = table(Selected_Intermediate[matrix(c(1:CTAChange,Index),,2)])
            All_elements = c(Radical,Selected_Intermediate[matrix(c(1:CTAChange,3-Index),,2)]-1)
            Radical = All_elements[(CTAChange+1):(length(All_elements))]
            Selected_CTA = sample(1:501,CTAChange,prob=CTAgent,replace=TRUE)
            Selected_CTA_table = table(Selected_CTA)
            Index = as.numeric(names(Selected_CTA_table))
            CTAgent[Index] = CTAgent[Index] - as.numeric(Selected_CTA_table)
            Index = as.numeric(names(CTA_Return)) # (1:501)
            CTAgent[Index] = CTAgent[Index] + as.numeric(CTA_Return)
            Intermediate_return = matrix(c(All_elements[1:(CTAChange)]+1,Selected_CTA),,2)
            incrementMatrix(Intermediate, Intermediate_return)
          }
          Temp.react = Temp.react+min(Temp.react)*0.00001
          Temp.react = cumsum(Temp.react)/sum(Temp.react)
          index = runif(1)
          
          if (index < Temp.react[1]){
            Index.1 = cpp_sample_one(Length_Radical)
            Index.2 = sample(0:500,1,prob=CTAgent) + 1 # Chain Length 
            Choosen_Radical =  Radical[Index.1] + 1
            Intermediate[Choosen_Radical,Index.2] = Intermediate[Choosen_Radical,Index.2] + 1
            Radical = Radical[-Index.1]
            CTAgent[Index.2] = CTAgent[Index.2] - 1
            
            Index.1 = cpp_sample_one(Length_Radical)
            Index.2 = sample(0:500,1,prob=CTAgent) + 1 # Chain Length 
            Choosen_Radical =  Radical[Index.1] + 1
            Intermediate[Choosen_Radical,Index.2] = Intermediate[Choosen_Radical,Index.2] + 1
            Radical = Radical[-Index.1]
            CTAgent[Index.2] = CTAgent[Index.2] - 1
            
            Initiator = Initiator-1
            Radical = c(Radical,0,0)
          } else if (index < Temp.react[2] & Length_Radical>=2){
            Index.1 = samplePoint(Intermediate)
            Intermediate[Index.1] = Intermediate[Index.1] - 1
            Index.2 = rbinom(1,1,0.5)+1
            Radical = c(Radical,Index.1[Index.2] - 1)
            CTAgent[Index.1[3-Index.2]] = CTAgent[Index.1[3-Index.2]] + 1
            
            Index.1 = samplePoint(Intermediate)
            Intermediate[Index.1] = Intermediate[Index.1] - 1
            Index.2 = rbinom(1,1,0.5)+1
            Radical = c(Radical,Index.1[Index.2] - 1)
            CTAgent[Index.1[3-Index.2]] = CTAgent[Index.1[3-Index.2]] + 1
            
            Index = cpp_sample_two(Length_Radical) # sample(Length_Radical,2,replace=FALSE)
            Length = sum(Radical[Index]) + 1
            Termination[Length+1] = Termination[Length+1] + 1
            Radical = Radical[-Index]
          }
          
          for (k in 1:10){
            Length_Radical = length(Radical)
            Temp.react = c(k_I*Initiator,
                           k_CTAgent*Length_Radical*sum(CTAgent),
                           k_CTADE  *sum(Intermediate),
                           k_T*(Length_Radical)^2)
            Temp.time = rexp(1,sum(Temp.react))
            Time_Current = Time_Current+Temp.time
            
            Propagation  = rpois(1, k_p*Length_Radical*Monomer*Temp.time)
            Monomer = Monomer - Propagation
            Index = sample(Length_Radical, Propagation, replace = TRUE)
            Index = factor(Index, levels = 1:Length_Radical)
            Counts = table(Index)
            Index = 1:Length_Radical
            Radical[Index] = Radical[Index] + as.numeric(Counts)
            P_Current = 1-Monomer/M
            
            Temp.react = Temp.react+min(Temp.react)*0.00001
            Temp.react = cumsum(Temp.react)/sum(Temp.react)
            index = runif(1)
            
            if (index < Temp.react[1]){
              Initiator = Initiator-1
              Radical = c(Radical,0,0)
            } else if (index < Temp.react[2]){
              Index.1 = cpp_sample_one(Length_Radical)
              Index.2 = sample(0:500,1,prob=CTAgent) + 1 # Chain Length 
              Choosen_Radical =  Radical[Index.1] + 1
              Intermediate[Choosen_Radical,Index.2] = Intermediate[Choosen_Radical,Index.2] + 1
              Radical = Radical[-Index.1]
              CTAgent[Index.2] = CTAgent[Index.2] - 1
            } else if (index < Temp.react[3]){
              Index.1 = samplePoint(Intermediate)
              Intermediate[Index.1] = Intermediate[Index.1] - 1
              Index.2 = rbinom(1,1,0.5)+1
              Radical = c(Radical,Index.1[Index.2] - 1)
              CTAgent[Index.1[3-Index.2]] = CTAgent[Index.1[3-Index.2]] + 1
            } else if (index < Temp.react[4] & Length_Radical>=2){
              Index = cpp_sample_two(Length_Radical) # sample(Length_Radical,2,replace=FALSE)
              Length = sum(Radical[Index]) + 1
              Termination[Length+1] = Termination[Length+1] + 1
              Radical = Radical[-Index]
            }
          }
        }
        Index.State = 1000*(i-1)+j+1
        System$Series$Time[Index.State] = Time_Current
        System$Series$Radical[Index.State] = length(Radical)
        System$Series$Termination[Index.State] = sum(Termination)
        Data = c(Termination[-c(1,2)],0)
        Counts = table(Radical)[-1]
        Index.Data = as.numeric(names(Counts))
        Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
        Counts = CTAgent[-1]
        Index.Data = which(Counts!=0)
        Counts = Counts[Counts != 0]
        Data[Index.Data] = Data[Index.Data] + as.numeric(Counts)
        Intermediate.table = sapply(2:1001, function(d) {sum(Intermediate[row(Intermediate) + col(Intermediate) == d + 1])}) # ingore the first term length=(0,0)
        Data = Data + Intermediate.table
        System$Series$X_n [Index.State] = sum(Data*1:1000)/sum(Data)
        System$Series$X_w [Index.State] = sum(Data*(1:1000)^2)/sum(Data*1:1000)
        System$Series$X_n2[Index.State] = sum(Termination[-1]*1:1000)/sum(Termination[-1])
        System$Series$X_w2[Index.State] = sum(Termination[-1]*(1:1000)^2)/sum(Termination[-1]*1:1000)
      }
    }
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
  Time$RAFT[i] = difftime(End_time, Start_time, units = "secs")
  save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_10^10.RData")
}
System$Series$PDI  = System$Series$X_w /System$Series$X_n
System$Series$PDI2 = System$Series$X_w2/System$Series$X_n2

System.RAFT = System
End_time = Sys.time()
Time$RAFT.Total = difftime(End_time, Start_time, units = "mins")
save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_10^10.RData")

##########
k_p*Length_Radical*Monomer,
A.CTAgent,
A.CTADE,

for(i in 1:(length(P_critical)-1)){
  cat("\n",P_critical[i]*100,"% ~",P_critical[i+1]*100,"%\n")   # "\014",
  pb=txtProgressBar(min=1,max=length(P_lists),style=3)
  for (j in 1:length(P_lists)){
    P_threshold = P_critical[i]+P_lists[j]
    setTxtProgressBar(pb, j)
    while (P_Current < P_threshold){
      print(Length_Radical)
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
  Time$RAFT[i] = difftime(End_time, Start_time, units = "secs")
}
System$Series$PDI  = System$Series$X_w /System$Series$X_n
System$Series$PDI2 = System$Series$X_w2/System$Series$X_n2

System.RAFT = System
End_time = Sys.time()
Time$RAFT.Total = difftime(End_time, Start_time, units = "mins")
save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_10^10_6.29.RData")

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
      
      Propagation  = rpois(1,k_p*Length_Radical*Monomer*Temp.time)
      CTAChange    = rpois(1,k_CTAgent*Length_Radical*CTA*Temp.time)
      Monomer = Monomer - Propagation
      
      if (Propagation!=0 & CTAChange!=0){
        # # Select and delete exchange Radical
        # Selected_CTA = sample(0:5000,CTAChange,prob=CTAgent,replace=TRUE)
        # All_elements = c(Radical,as.numeric(Selected_CTA))
        # Selected_CTA = table(Selected_CTA)
        # Index = as.numeric(names(Selected_CTA))
        # CTAgent[Index+1] = CTAgent[Index+1] - as.numeric(Selected_CTA)
        # 
        # 
        # Index.Propagation = sort(sample(1:(Propagation+CTAChange), CTAChange, replace=F))     
        # Index.CTAChange   =      sample(1:Length_Radical         , CTAChange, replace=T)      
        # Temp.Propagation  = c(rep(0,Length_Radical), Index.Propagation , rep(Propagation,Length_Radical))
        # Temp.CTAChange    = c(1:Length_Radical     , Index.CTAChange   , 1:Length_Radical)
        # 
        # 
        # Unique_elements <- unique(Temp.CTAChange)
        # S_mat     = sapply(Unique_elements, function(xxx) Temp.Propagation[which( Temp.CTAChange == xxx)[1:max(table( Temp.CTAChange))]]) # Survival_matrix
        # Index.mat = sapply(Unique_elements, function(xxx)                  which(Index.CTAChange == xxx)[1:max(table(Index.CTAChange))])
        # Index.mat[is.na(Index.mat)] = 0
        # S_mat[c(-1,-nrow(S_mat)), ] = S_mat[c(-1,-nrow(S_mat)), ] - Index.mat
        # S_mat                       = S_mat[-1,]-S_mat[-nrow(S_mat),]
        # S_Series = as.vector(na.omit(as.vector(t(S_mat)))) # 前Length_Radical个为原始自由基, 后面CTAChange个为Agent, 大小为存在时间(概率) #  Survival Series
        # 
        # # Propagation
        # Index = sample(length(S_Series),Propagation,prob=S_Series,replace = TRUE)
        # Index = factor(Index, levels = 1:(Length_Radical+CTAChange))
        # Counts = table(Index)
        # Index = 1:(Length_Radical+CTAChange)
        # All_elements[Index] = All_elements[Index] + as.numeric(Counts)
        # 
        # # Return Radical
        # Index = rev(c(1:Length_Radical,Index.CTAChange))
        # Index.Remain = Length_Radical+CTAChange+1-sapply(unique(Index), function(x) which(Index == x)[1])
        # Radical = All_elements[Index.Remain]
        # 
        # # Return CTAagent
        # CTAgent_Back = table(All_elements[-Index.Remain])
        # Index = as.numeric(names(CTAgent_Back))
        # CTAgent[Index+1] = CTAgent[Index+1] + as.numeric(CTAgent_Back)
        
        # Select and delete exchange Radical
        Selected_CTA = sample(0:5000,CTAChange,prob=CTAgent,replace=TRUE)
        All_elements = c(Radical,as.numeric(Selected_CTA))
        Selected_CTA = table(Selected_CTA)
        Index = as.numeric(names(Selected_CTA))+1
        CTAgent[Index] = CTAgent[Index] - as.numeric(Selected_CTA)
        
        # Propagation
        Index = sample(length(All_elements),Propagation,replace = TRUE)
        Index = factor(Index, levels = 1:length(All_elements))
        Counts = table(Index)
        Index = 1:(Length_Radical+CTAChange)
        All_elements[Index] = All_elements[Index] + as.numeric(Counts)
        
        # Return Radicals
        Radical = All_elements[(CTAChange+1):(length(All_elements))]
        # Return CTAngent
        CTAgent_Back = table(All_elements[1:(CTAChange+1)])
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
        for (k in 1:Propagation){
          Index.1 = cpp_sample_one(Length_Radical)
          Radical[Index.1] = Radical[Index.1] + 1
        }
      }
      
      P_Current = 1-Monomer/M
      
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
      print(Length_Radical)
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
    System$Series$X_n [Index.State] = sum(Data*1:5000)/sum(Data)
    System$Series$X_w [Index.State] = sum(Data*(1:5000)^2)/sum(Data*1:5000)
    System$Series$X_n2[Index.State] = sum(Termination[-1]*1:5000)/sum(Termination[-1])
    System$Series$X_w2[Index.State] = sum(Termination[-1]*(1:5000)^2)/sum(Termination[-1]*1:5000)
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
  Time$RAFT[i] = difftime(End_time, Start_time, units = "secs")
}
System$Series$PDI  = System$Series$X_w /System$Series$X_n
System$Series$PDI2 = System$Series$X_w2/System$Series$X_n2

System.RAFT = System
End_time = Sys.time()
Time$RAFT.Total = difftime(End_time, Start_time, units = "mins")

save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Quasi_Leaping_10^10_6.29.RData")



Leaping_condition = FALSE
for(i in 1:(length(P_critical)-1)){
  cat("\n",P_critical[i]*100,"% ~",P_critical[i+1]*100,"%\n")   # "\014",
  pb=txtProgressBar(min=1,max=length(P_lists),style=3)
  if (Leaping_condition) {
    for (j in 1:length(P_lists)){
      P_threshold = P_critical[i]+P_lists[j]
      setTxtProgressBar(pb, j)
      while (P_Current < P_threshold){
        Length_Radical = length(Radical)
        print(Length_Radical)
      }
    }
  } else {
    for (j in 1:length(P_lists)){
      P_threshold = P_critical[i]+P_lists[j]
      setTxtProgressBar(pb, j)
      Length_Radical = length(Radical)
      if (data[i] > 50) {
        program_A(data[i])
      } else {
        Leaping_condition = TRUE
        program_B(data[i])
        break
      }
    }
    if (Leaping_condition) {
      for (k in (j+1):length(data)) {
        program_B(data[j])
      }
    }
  }
}



# generate a 500*500 large matrix
set.seed(123)
matrix_size <- 500
mat <- matrix(runif(matrix_size^2, min = 0.1, max = 10), nrow = matrix_size)

Start_time = Sys.time()
mat_vector <- as.vector(mat)
x = replicate(1000,weighted_random_sample(mat))
End_time = Sys.time()
End_time-Start_time

Start_time = Sys.time()
y = replicate(1000,samplePoint(mat))
End_time = Sys.time()
End_time-Start_time
