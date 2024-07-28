# rm(list = ls())
load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical_7.09.RData")
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
# RAFT System ##################################################################
# I       -> 2R
# R + M   -> R
# R + CTA -> CTA + R
# R + R   -> T
options(digits = 3)

# 3592.528 77497.493 27027.021 19819.858     7.110
Start_time = Sys.time()
# RAFT System ##################################################################
load("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")
for (j in j:length(P_lists)) {
  P_threshold = P_critical[i]+P_lists[j]
  while (P_Current < P_threshold){
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
    cat("\r",P_Current)
    Temp.react = Temp.react+min(Temp.react)*0.00001
    Temp.react = cumsum(Temp.react)/sum(Temp.react)
    index = runif(1)
    
    if (index < Temp.react[1]){
      Initiator = Initiator-1
      Radical = c(Radical,0,0)
      # cat("\nI")
    } else if (index < Temp.react[2]){
      Index.1 = cpp_sample_one(Length_Radical)
      Index.2 = sample(0:500,1,prob=CTAgent) + 1 # Chain Length 
      Choosen_Radical =  Radical[Index.1] + 1
      Intermediate[Choosen_Radical,Index.2] = Intermediate[Choosen_Radical,Index.2] + 1
      Radical = Radical[-Index.1]
      CTAgent[Index.2] = CTAgent[Index.2] - 1
      # cat("\nC")
    } else if (index < Temp.react[3]){
      Index.1 = samplePoint(Intermediate)
      Intermediate[Index.1] = Intermediate[Index.1] - 1
      Index.2 = rbinom(1,1,0.5)+1
      Radical = c(Radical,Index.1[Index.2] - 1)
      CTAgent[Index.1[3-Index.2]] = CTAgent[Index.1[3-Index.2]] + 1
      # cat("\nD")
    } else if (index < Temp.react[4] & Length_Radical>=2){
      Index = cpp_sample_two(Length_Radical) # sample(Length_Radical,2,replace=FALSE)
      Length = sum(Radical[Index]) + 1
      Termination[Length+1] = Termination[Length+1] + 1
      Radical = Radical[-Index]
      # cat("\nT")
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
  save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")
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
save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")

for(i in 2:(length(P_critical)-1)){
  for (j in 1:length(P_lists)){
    P_threshold = P_critical[i]+P_lists[j]
    
    while (P_Current < P_threshold){
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
      cat("\r",P_Current)
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
    save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")
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
  save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")
}
System$Series$PDI  = System$Series$X_w /System$Series$X_n

System.RAFT = System
End_time = Sys.time()
Time$RAFT.Total = difftime(End_time, Start_time, units = "mins")
save.image("D:/Study/Code-R/Monte Carlo Simulation/Free Radical/RAFT/Critical.RData")