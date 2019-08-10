#This code has been developed by Arkasama Bandyopadhyay at the University of Texas at Austin in 2019
#This code represents a convex optimization model which minimizes costs and discomfort/inconvenience for a residential 
#customer by reducing and/or shifting loads from four end-use appliances- HVAC systems, electric water heaters, 
#electric vehicles, and pool pumps on the summer peak day under the effect of four different time-varying pricing structures
#This code accompanies the paper: "As one falls, another rises? Residential peak load reduction through
#electricity rate structures" by Arkasama Bandyopadhyay, Benjamin D. Leibowicz, Emily A. Beagle, and Michael E. Webber

#---------Input parameters----------
#Analysis is conducted at time steps of 15 minutes
delta_t <- (15*60) # seconds 

#--------Additional costs (costs other than the energy charge) ----------
customer_charge<- 10 #$/month
#Amortization is performed since analysis is for 1 day only
i <- 7/(365*100) #amortization rate
customer_charge_day <- 10*i*(1+i)^30/((1+i)^30-1) # in cents per day
PSA<-2.895# cents/kWh # power supply adjustment
CBC <- 0.154+0.124+0.335 #cents/kWh #community benefit charges
RC<-1.342 #cents/kWh #regulatory charges

#Overall cost for the case with real-time prices
#Note: Replace "RTP_summer_15_min" with constant rates/time-of-use rates/critical peak pricing as necessary
cost=RTP_summer_15_min*2.68+PSA+CBC+RC


#Rated power of each appliance-------
P_HVAC <- 3.5 #kW
P_EWH <- 4.5 #kW 
P_EV <- 6.6 #kW
P_PP <- 1.1 #kW 
#----------------

#-----Distribution line capacity bound----------
P_g_max <- 20 #kW #maximum capacity of the wires
#------------------------

#--------HVAC thermal properties---------
T_r_min <- (21+273) #K #customer-specified minimum room temperature
T_r_max <- (24+273) #K #customer-specified maximum room temperature
MaCa <- 30*10^6 #J/K #heat capacity of air inside the room
R_wall <- 0.005 #K/W #thermal resistance of the building envelop
T_r_in <- (22.5+273) #K #initial room temperature
COP <- 2.5 #coefficient of performance of the HVAC
#----------------------------------

#---------Discomfort parameter for HVAC-----------------
#We start analysis one time step before and continue analysis till one time step after what is needed
#to eliminate boundary effects (thus there are 98 intervals instead of 96)
no_of_dec_var <- 98 
alpha_HVAC <- matrix(0, nrow=no_of_dec_var, ncol=1)
for (i in 1:37){
  alpha_HVAC[i,1] <-100*291*4*0.25 *(9/5)^2/10^6 #$/degree C^2 or $/K^2  
}
for (i in 38:65){ #from 9 am to 4 pm
  alpha_HVAC[i,1] <-100*155*4*0.25*(9/5)^2/10^6 #$/degree C^2 or $/K^2 
}
for (i in 66:98){
  alpha_HVAC[i,1] <-100*291*4*0.25*(9/5)^2/10^6 #$/degree C^2 or $/K^2 
}
#-----------------------------------

#------------EWH thermal properties------------------
T_min_EWH <- (40+273) #K #customer-specified minimum water temperature
T_max_EWH <- (45+273) #K #customer-specified maximum water temperature
Cp <- 4186 #J/kgK # specific heat capacity of water 
rho <- 1000 #kg/m^3 #density of water
V <- 0.19 #m^3 #volume of EWH tank
const_F <- 2.6*10^-5 #m^3/s #flow rate of hot water
SA <- 0.18 #m^2 #surface area of tank
R <- 2.64 #m^2 K/W #thermal resistance of tank
T_water_in <- (32+273) #K #temperature of cold water coming in to the tank
T_EWH_in <-42.5+273 #initial temperature of water in EWH tank
alpha_EWH <- alpha_HVAC #same discomfort parameter as HVAC
#--------------------------------

#---------Iconvenience parameter for EV----------------------
alpha_EV <- matrix(0, nrow=no_of_dec_var, ncol=1)
for (i in 30:73){ #from 7 am to 6 pm
  alpha_EV[i,1] <-0.55 #$ 
}
for (i in 74:85){ #from 6 pm to 9 pm
  alpha_EV[i,1] <-0.55/3 #$ 
}
#-----------Efficiency and energy consumed for EV-----------------
eff_EV <-0.864 #roundtrip charging and discharging efficiency of EV 
E_EV_avg <- 11.79 #KWh #this value if from empirical data

#-----Inconvenience parameter for PP---------------------
alpha_PP<- matrix(0, nrow=no_of_dec_var, ncol=1)
for (i in 1:33){ #from 12 am to 8 am
  alpha_PP[i,1] <-0.55/3 #$ 
}
for (i in 90:98){ #from 10 pm to midnight
  alpha_PP[i,1] <-0.55/3 #$ 
}
#-----------Efficiency and energy consumed for PP-----------------
eff_PP <- 0.7 #efficiency of PP operation
E_consumed_avg <- 13.41 #kWh #House 5357
#-----------------------------

#customer-specified set point temperatures for HVAC and EWH
T_sp <- matrix((22.2+273), nrow=no_of_dec_var, ncol=1) # HVAC set point temperature
T_sp_EWH <- matrix((42+273), nrow=no_of_dec_var, ncol=1) #EWH set point temperature


#--------------------------------------------
#Objective function
#--------------------------------------------
#install.packages("CVXR")
library(CVXR)
#no of decision variables
no_of_dec_var <-(98) #no of 15 min prices for each day...98 instead of 96 because we start and end with an extra
#time step to eliminate boundary effects
# P_bought 98 times, [(T_room-T_sp), (T_water-T_sp_EWH), S_EV, S_PP each 98 times]
#objective function
mu<- matrix(0, nrow=5*no_of_dec_var, ncol=1)
#Linear portion of the objective function
#Use the following expression for mu if cost is varying over time
mu[,1] <- c(cost*delta_t/(100*3600), rep(0, (no_of_dec_var)), rep(0, (no_of_dec_var)), alpha_EV, alpha_PP)
#Use the following expression for mu if cost is constant over time
#mu[,1] <- c(rep((cost*delta_t/(100*3600)), no_of_dec_var), rep(0, (no_of_dec_var)), rep(0, (no_of_dec_var)), alpha_EV, alpha_PP)
#Quadratic portion of the objective function
sigma_interm<- matrix(0, nrow=5*no_of_dec_var, ncol=1)
sigma_interm <- c(rep(0, no_of_dec_var), alpha_HVAC, alpha_EWH, rep(0, (no_of_dec_var)), rep(0, (no_of_dec_var)))
sigma <- diag(sigma_interm)
w=Variable(5*no_of_dec_var)
part1<- t(mu) %*% w
part2 <- quad_form(w, sigma)
obj <- part1+part2

#------Upper and lower bounds for the decision variables-------------
upper_bound <- matrix(0, nrow=5*no_of_dec_var, ncol=1)
upper_bound[,1] <- c(rep(P_g_max, no_of_dec_var), (rep(T_r_max, (no_of_dec_var))-T_sp), (rep(T_max_EWH, (no_of_dec_var))-T_sp_EWH), rep(1, (no_of_dec_var)), rep(1, (no_of_dec_var)))
lower_bound <- matrix(0, nrow=5*no_of_dec_var, ncol=1)
lower_bound[,1] <- c(rep(0, no_of_dec_var), (rep(T_r_min, (no_of_dec_var))-T_sp), (rep(T_min_EWH, (no_of_dec_var))-T_sp_EWH), rep(0, (no_of_dec_var)), rep(0, (no_of_dec_var)))

#---------loading required packages--------------------
#install.packages("readxl")
library(readxl)

#Uncontrollable power profile is obtained from  empirical 1-minute interval data 
P_uncontrol_15_min <- matrix(0, nrow=no_of_dec_var, ncol=1)
k=1
#Conversion to 15-minute interval data
for(i in 1:no_of_dec_var){
  P_uncontrol_15_min[i,1] <-round(sum(P_uncontrol[k:(15*i)]/15),2) #in KW
  k=k+15
}

#Solar generation profile is obtained from empirical 1-minute interval data 
P_solar_gen_15_min <- matrix(0, nrow=no_of_dec_var, ncol=1)
k=1
#Conversion to 15-minute interval data
for(i in 1:no_of_dec_var){
  P_solar_gen_15_min[i,1] <-round(sum(P_solar_gen[k:(15*i)]/15),2) #in KW
  k=k+15
}
# 

#--------------------------------------------
#Constraints
#--------------------------------------------


#--------------------------------------------
#Constraint 1 - load constraint (sum of power flowing into the home from the grid and solar generation
#must be greater than or equal to power usage in the home at each time step)
#--------------------------------------------
#1 matrix for each kind of decision variable
mat_cost <- diag (-1, nrow= no_of_dec_var, ncol = no_of_dec_var)

#---------Switching from S_HVAC to (T_r - T_sp) as the decision variable----------
first_term = (1-(delta_t/(MaCa*R_wall))) #A
second_term = delta_t/(MaCa*R_wall) #B
third_term <- COP*P_HVAC*1000*delta_t/MaCa #C

mat_HVAC <- diag ((-P_HVAC/third_term), nrow= no_of_dec_var, ncol = no_of_dec_var)
for (i in 2:no_of_dec_var){
  mat_HVAC[i, (i-1)] <- first_term*P_HVAC/third_term
}

#Note: Obtain T_amb (ambient temperature) from empirical data
#matrix for min and max
mat_amb <- matrix(0, nrow <- no_of_dec_var, ncol <- 1)
mat_amb[1,1] <- (-P_HVAC*first_term*T_r_in/third_term) -P_HVAC*second_term *T_amb[1]/third_term
for (i in 2:no_of_dec_var){
  mat_amb[i,1] <- -P_HVAC*second_term *T_amb[i]/third_term
}
mat_amb_2 <- mat_amb -mat_HVAC%*%T_sp


#------Switching from S_EWH to (T_EWH-T_sp_EWH) as the decision variable-----------
#-----one-element EWH model--------
C=Cp*rho*V
B= rho*const_F*Cp
G= SA/R
R_prime <- 1/(G+B)

first_term_EWH <- G*R_prime #B
second_term_EWH <- B*R_prime*T_water_in #C
third_term_EWH <- 0.9*P_EWH *1000* R_prime#D
fourth_term_EWH <- 1-exp(-delta_t/(R_prime*C)) #E
fifth_term_EWH <- exp(-delta_t/(R_prime*C)) #A

mat_EWH <- diag(P_EWH/(third_term_EWH*fourth_term_EWH), nrow <- no_of_dec_var, ncol <- no_of_dec_var)

for (i in 2:no_of_dec_var){
  mat_EWH[i, (i-1)] <- -P_EWH*fifth_term_EWH/(third_term_EWH*fourth_term_EWH)
}
mat_amb_EWH <- matrix(0, nrow <- no_of_dec_var, ncol <- 1)
mat_amb_EWH[1,1] <-(P_EWH*T_EWH_in* fifth_term_EWH/(third_term_EWH*fourth_term_EWH))+(P_EWH*first_term_EWH*T_amb[1]/third_term_EWH)+(P_EWH*second_term_EWH/third_term_EWH)
for (i in 2:no_of_dec_var){ 
  mat_amb_EWH[i,1] <- (P_EWH*first_term_EWH*T_amb[i]/third_term_EWH)+(P_EWH*second_term_EWH/third_term_EWH)
}

mat_amb_EWH_2 <- mat_amb_EWH -mat_EWH%*%T_sp_EWH
#--------------------------------------------------------------
mat_EV <- diag (P_EV, nrow= no_of_dec_var, ncol = no_of_dec_var)
mat_PP <- diag (P_PP, nrow= no_of_dec_var, ncol = no_of_dec_var)

#LHS of Constraint 1
mat_constraint_1 <- cbind(mat_cost, mat_HVAC,mat_EWH, mat_EV,mat_PP)

#RHS of Constraint 1
mat_1_RHS <- -P_uncontrol_15_min + mat_amb_2+mat_amb_EWH_2+P_solar_gen_15_min
#Note: Neglect P_solar_gen_15_min if the sample home does not have solar panels

#--------------------------------------------
#Constraint 2 - one parameter thermal model HVAC
#--------------------------------------------

first_term = (1-(delta_t/(MaCa*R_wall))) 
second_term = delta_t/(MaCa*R_wall)
third_term <- COP*P_HVAC*1000*delta_t/MaCa 

mat_HVAC <- diag ((-1/third_term), nrow= no_of_dec_var, ncol = no_of_dec_var)
for (i in 2:no_of_dec_var){
  mat_HVAC[i, (i-1)] <- first_term/third_term
}

#matrix for min and max
mat_amb <- matrix(0, nrow <- no_of_dec_var, ncol <- 1)
mat_amb[1,1] <- (-first_term*T_r_in/third_term) -second_term *T_amb[1]/third_term
for (i in 2:no_of_dec_var){
  mat_amb[i,1] <- -second_term *T_amb[i]/third_term
}
mat_amb_2 <- mat_amb -mat_HVAC%*%T_sp

mat_0 =matrix(0, ncol=no_of_dec_var, nrow=no_of_dec_var)
mat_HVAC_interm <- cbind(mat_0, mat_HVAC, mat_0,mat_0,mat_0)
mat_HVAC<- cbind(mat_0, mat_HVAC, mat_0,mat_0,mat_0)

mat_LHS <- matrix(0, nrow <- no_of_dec_var, ncol <-1)+mat_amb_2

mat_RHS <- matrix(1, nrow <- no_of_dec_var, ncol <-1)+mat_amb_2

#--------------------------------------------
#Constraint 3 - one-element EWH model
#--------------------------------------------
#-----EWH properties--------
C=Cp*rho*V
B= rho*const_F*Cp
G= SA/R
R_prime <- 1/(G+B)

first_term_EWH <- G*R_prime #B
second_term_EWH <- B*R_prime*T_water_in #C
third_term_EWH <- 0.9*P_EWH*1000* R_prime#D
fourth_term_EWH <- 1-exp(-delta_t/(R_prime*C)) #E
fifth_term_EWH <- exp(-delta_t/(R_prime*C)) #A

mat_S_EWH <- diag((1/(third_term_EWH*fourth_term_EWH)), nrow <- no_of_dec_var, ncol <- no_of_dec_var)

for (i in 2:no_of_dec_var){
  mat_S_EWH[i, (i-1)] <- -fifth_term_EWH/(third_term_EWH*fourth_term_EWH)
}

mat_EWH_interm <- cbind(mat_0, mat_0, mat_S_EWH,mat_0,mat_0)

mat_EWH <-cbind(mat_0, mat_0, mat_S_EWH,mat_0,mat_0)

mat_amb_EWH <- matrix(0, nrow <- no_of_dec_var, ncol <- 1)
mat_amb_EWH[1,1] <-(T_EWH_in* fifth_term_EWH/(third_term_EWH*fourth_term_EWH))+(first_term_EWH*T_amb[1]/third_term_EWH)+(second_term_EWH/third_term_EWH)
for (i in 2:no_of_dec_var){ 
  mat_amb_EWH[i,1] <- (first_term_EWH*T_amb[i]/third_term_EWH)+(second_term_EWH/third_term_EWH)
}

mat_amb_EWH_2 <- mat_amb_EWH-mat_S_EWH%*%T_sp_EWH

mat_T_w_min <- matrix(0, nrow<- no_of_dec_var, ncol <-1)

mat_LHS_EWH <- mat_T_w_min+mat_amb_EWH_2

mat_RHS_EWH <- matrix(1, nrow <- no_of_dec_var, ncol <-1)+mat_amb_EWH_2

#--------------------------------------------
#Constraint 4: EV charging model
#--------------------------------------------

mat_S_EV <- matrix(1, nrow = 1, ncol=no_of_dec_var)

mat_0 = matrix(0, nrow=1, ncol=no_of_dec_var)
mat_constraint_4_interm <-cbind(mat_0,mat_0, mat_0, mat_S_EV, mat_0)

mat_RHS_EV <- (E_EV_avg+0.1*E_EV_avg)/(P_EV*eff_EV*delta_t/3600)#10% tolerance interval
mat_LHS_EV <- (E_EV_avg-0.1*E_EV_avg)/(P_EV*eff_EV*delta_t/3600)#10% tolerance interval


#--------------------------------------------
#Constraint 5: PP operation model
#--------------------------------------------

mat_S_PP <- matrix(1, nrow = 1, ncol=no_of_dec_var)

mat_0 = matrix(0, nrow=1, ncol=no_of_dec_var)
mat_constraint_5_interm <-cbind(mat_0,mat_0, mat_0, mat_0, mat_S_PP)

mat_PP_RHS <- (E_consumed_avg+0.1*E_consumed_avg)/(P_PP*eff_PP*delta_t/3600)
mat_PP_LHS <- (E_consumed_avg-0.1*E_consumed_avg)/(P_PP*eff_PP*delta_t/3600)

#list of all constraints
constr <- list(w >= lower_bound, w<=upper_bound, mat_constraint_1%*%w<=mat_1_RHS, mat_HVAC%*%w<=mat_RHS, mat_HVAC%*%w>=mat_LHS,mat_EWH%*%w<=mat_RHS_EWH, mat_EWH%*%w>=mat_LHS_EWH,mat_constraint_4_interm %*%w>=mat_LHS_EV, mat_constraint_4_interm %*%w<=mat_RHS_EV, mat_constraint_5_interm%*%w>=mat_PP_LHS, mat_constraint_5_interm%*%w<=mat_PP_RHS)
#Minimize the objective function
prob <- Problem(Minimize(obj), constr)
result <- solve(prob)

#-----Analyzing the solutions-----------
solution=result$getValue(w)

#sol_P_bought <- solution[1:98]
sol_P_bought <- solution[2:97]

sol_T_r_original <- solution[99:196]+ T_sp
sol_T_r <- solution[100:195]+T_sp[2:97]

sol_T_EWH_original<- solution[197:294]+T_sp_EWH
sol_T_EWH <- solution[198:293]+T_sp_EWH[2:97]

#sol_S_EV_RTP <- solution[295:392]
sol_S_EV_RTP <- solution[296:391]

#sol_S_PP_RTP <- solution[393:490]
sol_S_PP_RTP <- solution[394:489]

#------Plotting the solutions-----
plot(sol_P_bought) 
plot(sol_T_r-273)
plot(sol_T_EWH-273)
plot(sol_S_EV_RTP)
plot(sol_S_PP_RTP)


#-----Operational level of HVAC---------
sol_S_HVAC_RTP= matrix(0, nrow=98, ncol=1)
sol_S_HVAC_RTP[1]= (first_term*T_r_in + T_amb[1]*second_term-sol_T_r_original[1])/third_term
for (i in 2:98){
  sol_S_HVAC_RTP[i]=(first_term*sol_T_r_original[i-1] + T_amb[i]*second_term-sol_T_r_original[i])/third_term
}


#------Operational level of EWH---------
sol_S_EWH_RTP= matrix(0, nrow=98, ncol=1)
sol_S_EWH_RTP[1]= (sol_T_EWH_original[1] - fifth_term_EWH*T_EWH_in - first_term_EWH*T_amb[1]*fourth_term_EWH-second_term_EWH*fourth_term_EWH)/(third_term_EWH*fourth_term_EWH)
for (i in 2:98){
  sol_S_EWH_RTP[i]= (sol_T_EWH_original[i] - fifth_term_EWH*sol_T_EWH_original[i-1] - first_term_EWH*T_amb[i]*fourth_term_EWH-second_term_EWH*fourth_term_EWH)/(third_term_EWH*fourth_term_EWH)
  
}


