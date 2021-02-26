#remove previous items
rm(list= ls())

#for nice plotting
#library(plotly)

### Parameters
# Elasticty parameter mmust be less than 1
a <- 0.5
#Number of cities
V <- 30
#Number of total steps
T_max <- 1000
#Number of final industries
N <- 100
#sensitivity to profit
m <- 0.1
#Geographical distance cost coefficient
g <- 1
# Circle length
L <- 2
#fixed cost
fixed <- 5
#technology elasticity
t_e <- 1
#propensity to R&D
prop_rd <- 0.05
#probability of discovering a new technology
p_disc <- 0.2
#Wages
w <- rep(1, V)


### Variables

#distance matrices (length of circle = L)
#distance between cities
geo_distance <- pmin(abs(matrix(rep(seq(0, L - L/V, length.out = V), V), nrow = V, ncol = V) - matrix(rep(seq(0, L - L/V, length.out = V), V), nrow = V, ncol = V, byrow = TRUE)),
                 L - abs(matrix(rep(seq(0, L - L/V, length.out = V), V), nrow = V, ncol = V) - matrix(rep(seq(0, L - L/V, length.out = V), V), nrow = V, ncol = V, byrow = TRUE)))


#distance between industries

ind_prox <- matrix(0, N, N)
'
N_mod <- 20
for (x in (0:(N_mod-1))) {
  
  ind_prox[ seq(x* N/N_mod - (2 * (x != 0)), (x+1)* N/N_mod, 1) ,
            seq(x* N/N_mod - (2 * (x != 0)) , (x+1)* N/N_mod, 1) ] <- 1
  
  
}
'
ind_prox <- mapply(function(x,y) 1*((y-x) <= 2),
                   expand.grid(x = 1:N, y= 1:N)$x,
                   expand.grid(x = 1:N, y= 1:N)$y,
                   SIMPLIFY = TRUE)
ind_prox <- pmin(matrix(ind_prox, nrow = N, ncol = N), 1)
ind_prox[lower.tri(ind_prox)] <- 0
ind_prox <- ind_prox + t(ind_prox)
diag(ind_prox) <- 1
sum(ind_prox - t(ind_prox))

#colnames(ind_prox) <- rownames(ind_prox) <- paste("Ind", 1:N)

#write.table(ind_prox,
#          "/Volumes/Donnees/switchdrive/Mehdi-Celine/Chapitre - Theories and models of urbanization/technology_space/adj_mat.csv",
#          row.names = TRUE, col.names = TRUE, sep = ";", quote = FALSE)

image(ind_prox)

#ind_prox <- matrix(runif(N^2), nrow = N, ncol = N)
#ind_prox[upper.tri(ind_prox)] <- ind_prox[lower.tri(ind_prox)]
#ind_prox <- ind_prox + diag(1, N) - diag(diag(ind_prox))

#Indicates whether the industry was discovered or not
#discovered <- c(rep(TRUE,V), rep(FALSE, N-V))
discovered <- c(TRUE, rep(FALSE, N-1))
##incomes
#Initial Incomes of each city
inc <- rep(500, V)
#Income share of each industry, changes with the entry of new industries
ind_inc_share <- rep(0, N) ##For the moment no industries

#indicates the presence of an industry in a city, to use to calculate probabilites of entries
#we approximate the max with an n norm, with n very big
ind_city <- matrix(0, nrow = N, ncol = V)
diag(ind_city) <- 1

#matrix indicating for each city and technology, the proximity of the closest tech
#we take the max with 1 to correct approximation errors
max_ind_proximity_city <- pmin((ind_prox^150 %*% ind_city)^(1/150) , 1)

##Creating vectors of firms (for each city there is initially a local industry)

#Initializing vector of firms activities 
active <- c(rep(TRUE, V), rep(FALSE, N*V + as.integer(N*V*100)))

#Firms industries
firm_ind <- c(rep(1,V),
              rep(0, N*V + as.integer(N*V*100))) ## as.vector collapses by columns by default

#Firms cities
firm_city <- c(1:V, rep(0, N*V + as.integer(N*V*100)))

#Initializing firms technologies
tech <- c(rep(1, V), rep(0, N*V + as.integer(N*V*100)))

##Matrix giving the entry time of each industry
entry_time <- c(rep(1, V), rep(Inf, N*V + as.integer(N*V*100)))


#Creating the distance cost between firms and cities matrix
firm_var_cost <- matrix(10^50, nrow = length(active), ncol = V)
firm_var_cost[active,] <-  matrix(rep(w[firm_city[active]], V), nrow= length(firm_city[active]), ncol = V) * 
                  ((1+geo_distance[firm_city[active], 1:V])^g - 1 + matrix(rep(1/tech[active]^t_e, V), nrow= length(firm_city[active]), ncol = V))

##Creating the matrix of prices (producers in lines, (industries, cities))
price <- matrix(10^50, nrow = length(active), ncol = V)
## Creating the matix of the quantities
quantity <- matrix(0, nrow = length(active), ncol = V)

##Initializing the variation index for firms
var_index <- matrix(0, nrow = length(active), ncol = V)

###saving previous states
#firm's variable cost
firm_var_cost_prec <- firm_var_cost
#profit
profit_prec <-  matrix(0, nrow = length(active), ncol = V)  ##Only for computation
#price
price_prec <- price

####Optimization and testing

##For testing only
j <- 1
t_stop <- 0
#test speed
ptm <- proc.time()

plot(1:V, inc, xlim = c(1,V), ylim = c(10, 10000), log= "xy",
     xaxt = "s", yaxt = "s", type = "l", col = rainbow(20)[15])

##########Starting the iterations
while (j < T_max & t_stop < 50) {
  
  if(all(discovered)) {t_stop <- t_stop + 1}
  if(sum(inc > 0.1) <= 1) {t_stop <- 150}
  
  ###Deciding on new entrants
  if (j > 1) {
  
  #Reinitialize entry matrix
  entry_prob <- matrix(1, nrow = N, ncol = V) * matrix(rep(discovered, V), nrow = N, ncol = V, byrow = FALSE)
    
  #Compute closest industry to each city
  max_ind_proximity_city <- pmin((ind_prox^150 %*% ind_city)^(1/150) , 1)
  
  #Compute entry probs
  entry_prob <- max_ind_proximity_city
  
  entry_prob[!discovered,] <- p_disc*max_ind_proximity_city[!discovered,]
  
  entring <- matrix(sapply(as.vector(entry_prob), rbinom, n=1, size=1), nrow = N, ncol = V, byrow = FALSE)
  
  if(sum(entring) > 0) {
  
  #update the firms vectors
  entry_time[entry_time == Inf][1:sum(entring)] <- j
  
  firm_ind[entry_time == j] <- unlist(mapply(rep, 1:N, apply(entring, 1, sum), USE.NAMES = FALSE))
  
  firm_city[entry_time == j] <- Filter(function(x) x>0,
                                       as.vector(t(entring*matrix(rep(1:V,N), nrow = N, ncol = V, byrow = TRUE))))
  
  
  ##compute technologies of entrants
  if (any(active)) {
  max_tech_city <- aggregate.data.frame(tech[active], 
                             by = list(city = firm_city[active], ind = firm_ind[active]),
                             FUN = max)
  max_tech_city <- mapply(function(ind, city) max_tech_city[max_tech_city$city == city & max_tech_city$ind == ind,3],
                          rep(1:N, V), unlist(lapply(1:V, rep, N)), SIMPLIFY = TRUE)
  max_tech_city[sapply(max_tech_city, length) == 0] <- 1 ####Here careful
  max_tech_city <- matrix(unlist(max_tech_city), nrow = N, ncol = V, byrow = FALSE)
  
  ##Determine the for each city and technology the closest industry in the city (returns index)
  closest_ind_city <- mapply(function(ind, city) which.max(ind_city[,city] * ind_prox[,ind]),
                             rep(1:N, V), unlist(lapply(1:V, rep, N)))
  closest_ind_city[sapply(closest_ind_city, length) == 0] <- 0
  closest_ind_city <- matrix(unlist(closest_ind_city), nrow = N, ncol = V, byrow = FALSE)
  
  #For each city, ind, gives the technology of the closest present industry
  max_tech_closest_city <-  mapply(function(ind, city) max_tech_city[ind, city],
                                   as.vector(closest_ind_city), unlist(lapply(1:V, rep, N)))
  max_tech_closest_city[sapply(max_tech_closest_city, length) == 0] <- 1
  max_tech_closest_city <- matrix(unlist(max_tech_closest_city), nrow = N, ncol = V, byrow = FALSE)
  
  #computing the mode of the tech
  entrants_tech <- 1 + mapply(rbinom, as.integer(max_tech_closest_city) - 1, as.vector(max_ind_proximity_city), n = 1)
  entrants_tech <- matrix(entrants_tech, nrow = N, ncol = V) * entring
  entrants_tech[!discovered,] <- entring[!discovered,]
        
  #computing the techs
  tech[entry_time == j] <- Filter(function(x) x > 0, as.vector(t(entrants_tech)))
                                 
  
  } else {
    tech[entry_time == j] <- 1
  }
  }
  discovered <- discovered | (rowSums(entring) > 0)
  }
  
  #Initialize the prices for entrants
  if (any(entry_time ==  j)) {
  
  #Change their state to active  
  active[entry_time ==  j] <- TRUE
  
  ##Change of the income share to account for the new entring indusries
  ind_inc_share[unique(firm_ind[active])] <- 1/length(unique(firm_ind[active]))
  
  firm_var_cost[entry_time ==  j,] <-  matrix(rep(w[firm_city[entry_time ==  j]], V), nrow= length(firm_city[entry_time ==  j]), ncol = V) * 
    ((1+geo_distance[firm_city[entry_time ==  j], 1:V])^g - 1 + matrix(rep(1/tech[entry_time ==  j]^t_e, V), nrow= length(firm_city[entry_time ==  j]), ncol = V))
  
  
  price[entry_time == j,] <-  firm_var_cost[entry_time == j,] * (1+m)
                                
  }
  
  #set quantities demanded
  #######Remark : rowsum sums the columns follwing the groups, so the vector of group has length nrow
  quantity[active,] <- price[active,]^(a-1) / 
                       rowsum(price[active,]^a, group = as.character(firm_ind[active]), reorder = FALSE)[as.character(firm_ind[active]),] *
                       (ind_inc_share[firm_ind[active]] %*% t(inc))
  
  #Calculate variation of price index
  var_index[active,] <- (quantity[active,]*(price[active,] - firm_var_cost_prec[active,]) - profit_prec[active,]) *
                        log(price[active,]/price_prec[active,])
  var_index[entry_time == j,] <- 100
  
  ##Storing old variables
  #Store old profit (it is not exactly the profit, its actually profit + fixed cost )
  profit_prec[active,] <-  quantity[active,]*(price[active,] - firm_var_cost[active,])
  #store old cost
  firm_var_cost_prec[active,] <- firm_var_cost[active,]
  
  #update incomes
  inc[firm_city[active][!duplicated(firm_city[active])]] <- rowsum((quantity[active,]*price[active,]) %*% rep(1, V), firm_city[active], reorder = FALSE)
  inc[setdiff(1:V,firm_city[active])] <- 0 
  
  price_prec[active,] <- price[active,]
  
  #update prices for incumbents
  price[active,] <- pmin(firm_var_cost[active,] + (price[active,] - firm_var_cost[active,]) * (1+m)^(2*atan(pi*var_index[active,]/2)/pi), 10^50)
  
  
  #firms with negative profit exit
  active[active] <- (rowSums(profit_prec[active,]) - fixed ) >= 0 
  
  #update technologies
  
  tech[active] <- tech[active] + prop_rd * pmax(rowSums(profit_prec[active,]) - fixed, 0)
                  
                  
  #update cost
  firm_var_cost[active,] <-  matrix(rep(w[firm_city[active]], V), nrow= length(firm_city[active]), ncol = V) * 
    ((1+geo_distance[firm_city[active], 1:V])^g - 1 + matrix(rep(1/tech[active]^t_e, V), nrow= length(firm_city[active]), ncol = V))
  
  #change the income share of non active industries to 0 and update other to 1/#industries
  ind_inc_share <- rep(0, N)
  ind_inc_share[unique(firm_ind[active])] <- 1/length(unique(firm_ind[active]))
  
  ##Upate the industries present in each city
  ind_city_current <- aggregate.data.frame(active[active], by = list(ind = firm_ind[active], city = firm_city[active]), FUN = any)
  ind_city_current <- mapply(function(indus, cit) ind_city_current[ind_city_current$city == cit & ind_city_current$ind == indus, 3],
                      rep(1:N, V), unlist(lapply(1:V, rep, N)), SIMPLIFY = TRUE)
  ind_city_current[sapply(ind_city_current, length) == 0] <- 0
  ind_city <- matrix(unlist(ind_city_current), nrow = N, ncol = V, byrow = FALSE)
  
  ###Clean memory every n steps to speed up the processing
  #if (j %% 20 == 0 ) {
  
  #(j %in% as.integer(T_max*seq(0.1, 1, 40/T_max)^3)) {
    
    price <- price[entry_time >= j | active, ]
    quantity <- quantity[ entry_time >= j | active, ]
    firm_var_cost <- firm_var_cost[entry_time >= j | active, ]
    firm_var_cost_prec <- firm_var_cost_prec[entry_time >= j | active, ]
    var_index <- var_index[entry_time >= j | active,]
    profit_prec <- profit_prec[entry_time >= j | active, ]
    price_prec <- price_prec[entry_time >= j | active, ]
    
    firm_city <- firm_city[entry_time >= j | active]
    firm_ind <- firm_ind[entry_time >= j | active]
    tech <- tech[entry_time >= j | active]
    
    active_prec <- active
    active <- active[entry_time >= j | active]
    entry_time <- entry_time[entry_time >= j | active_prec]
    
    lines(1:sum(inc > 0), sort(inc[inc>0], decreasing = TRUE), col = rainbow(20)[15 - j/20])
    
  #}
  
 ##extend the firms vectors in case there is no more space
  if(sum(entry_time == Inf) <  N*V) {
    
    ###State vectors
    n_row_to_add <- N*V + as.integer(N*V*100)
    
    active <- c(active, rep(FALSE, n_row_to_add))
    
    firm_ind <- c(firm_ind, rep(0, n_row_to_add)) ## as.vector collapses by columns by default
    
    firm_city <- c(firm_city, rep(0, n_row_to_add))
    
    tech <- c(tech, rep(0, n_row_to_add))
    
    entry_time <- c(entry_time, rep(Inf, n_row_to_add))
    
    ##State matrices
    var_index <- rbind(var_index, matrix(0, nrow = n_row_to_add, ncol = V))
    
    quantity <- rbind(quantity, matrix(0, nrow = n_row_to_add, ncol = V))
    
    price <- rbind(price, matrix(10^50, nrow = n_row_to_add, ncol = V))
    
    firm_var_cost <- rbind(firm_var_cost, matrix(0, nrow = n_row_to_add, ncol = V))
    
    firm_var_cost_prec <- rbind(firm_var_cost_prec, matrix(0, nrow = n_row_to_add, ncol = V))
    
    profit_prec <- rbind(profit_prec, matrix(0, nrow = n_row_to_add, ncol = V))
    
    price_prec <- rbind(price_prec, matrix(0, nrow = n_row_to_add, ncol = V))
      }
  j <- j + 1
}

ptm <- proc.time() - ptm
ptm

###Results verctors for analysis
##Size

size_regression <- summary(lm(log(1:(length(inc[inc > 0]))) ~ log(sort(inc[inc > 0], decreasing = TRUE))))

size_cofficient <- size_regression$coefficients[2,"Estimate"]

size_explained_var <- size_regression$r.squared

###Specialization
##We use the chi^2 index
ind_city_size <- rowSums(profit_prec[active,])
ind_city_size <- aggregate(ind_city_size - fixed,
                           by = list(city = firm_city[active], ind = firm_ind[active]),
                           FUN = sum, drop = FALSE)
ind_city_size <- as.numeric(mapply(function (x,y) ind_city_size[ind_city_size$city == x & ind_city_size$ind == y, 3],
                        rep(1:V, N), unlist(lapply(1:N, rep, times = V))))
ind_city_size[is.na(ind_city_size)] <- 0
ind_city_size <- matrix(ind_city_size, nrow = N, ncol = V, byrow = TRUE)

ind_size <- apply(ind_city_size, 1, sum)
city_size <- apply(ind_city_size, 2, sum)

#Removing empty cities and industries
ind_city_size <- ind_city_size[ind_size != 0,]
ind_city_size <- ind_city_size[,city_size != 0]
ind_size <- ind_size[ind_size != 0]
city_size <- city_size[city_size != 0]

##Transforming the sizes into suitable ratios
ind_city_size <- ind_city_size %*% diag(1/city_size)

##Herfindahl index
herfindal <- colSums(ind_city_size^2)
plot(city_size, herfindal, log = "xy")
plot( (j - entry_time[active])/j, rowSums(profit_prec[active,]) - fixed)
plot( (j - entry_time[active])/j, tech[active])


specialization_regression <- summary(lm(log(herfindal) ~ log(city_size)))
specialization_cofficient <- specialization_regression$coefficients[2,"Estimate"]
specialization_explained_var <- specialization_regression$r.squared


#Plotting cities sizes

library(ggplot2)

#pdf("/Volumes/Donnees/switchdrive/Mehdi-Celine/Chapitre - Theories and models of urbanization/city_final.pdf",
#    width = 10, height = 10)
par(mar=c(1,1,1,1))
plot(cos(2*pi/V * 1:V), sin(2*pi/V * 1:V), cex = sqrt(inc)/12, xlim = c(-1.5,1.5),
     ylim = c(-1.5,1.5), xaxt='n', yaxt='n', xlab = "", ylab = "", pch = 16)

lines(cos(2*pi*(1:5000) / 5000), sin(2*pi*(1:5000) / 5000), lwd = 0.1)
#ev.off()

##Any points to remove??
n_remove <- 0

size_regress <- log(sort(inc[inc > 0], decreasing = TRUE))

#For the regression
reg_plot <- summary(lm(size_regress[1:(length(size_regress)-n_remove)]  ~ log(1:(length(size_regress)-n_remove)) ))
x_line_to_plot <- seq(1, (V-n_remove)*2, by = 0.1)
y_line_to_plot <- exp(reg_plot$coefficients[1]) * x_line_to_plot^reg_plot$coefficients[2]

##plot the points
leg_lab <- paste0("log(y) = ", round(reg_plot$coefficients[1],2), " - ", -round(reg_plot$coefficients[2],2), " log(x) ")

hierarchy <- ggplot(data.frame(Rank = 1:(V-n_remove), Size = exp(size_regress[1:(V-n_remove)])), aes(x=Rank, y=Size)) + 
  theme_minimal() +
  geom_line(aes(x = x, y = y, 
                color = leg_lab), 
            data = data.frame(x = x_line_to_plot, y = y_line_to_plot),
            size = 1) +
  scale_x_continuous(trans='log10', limits = c(1, 1.2*length(size_regress))) +
  scale_y_continuous(trans='log10', limits = c(0.5*min(exp(size_regress)), 2*max(exp(size_regress)))) +
  geom_point() +
  labs(x = "Rank", y = "City size", color = "") +
  theme(axis.text.x = element_text(color="black", 
                                   size=16, angle=0),
        axis.text.y = element_text(color="black", 
                                   size=16, angle=0),
        axis.title=element_text(size= 20),
        legend.position = c(0.7, 0.8),
        legend.text = element_text(size=18))
  

hierarchy

#ggsave(paste0("/Volumes/Donnees/switchdrive/Mehdi-Celine/Chapitre - Theories and models of urbanization/simulation_rank_size.pdf"), hierarchy,
#       device = "pdf", width = 20, height = 20, units = "cm")


#############################Specialization
##To plot the regression line
reg_plot <- summary(lm(log(herfindal)  ~ log(city_size)))
x_line_to_plot <- seq(0.1, 1.2*max(city_size), by = 0.1)
y_line_to_plot <- exp(reg_plot$coefficients[1]) * x_line_to_plot^reg_plot$coefficients[2]


leg_lab <- paste0("log(y) = ", round(reg_plot$coefficients[1],2), " - ", -round(reg_plot$coefficients[2],2), " log(x) ")

spec <- ggplot(data.frame(Size = city_size, Spec = herfindal), aes(x=Size, y=Spec)) + 
  theme_minimal() +
  geom_line(aes(x = x, y = y, 
                color = leg_lab), 
            data = data.frame(x = x_line_to_plot, y = y_line_to_plot),
            size = 1) +
  scale_x_continuous(trans='log10', limits = c(min(city_size)/1.2, max(city_size)*1.2)) +
  scale_y_continuous(trans='log10', limits = c(min(herfindal),1)) +
  geom_point() +
  labs(x = "City size", y = "Hirschmann-Herfindahl index", color = "") +
  theme(axis.text.x = element_text(color="black", 
                                   size=16, angle=0),
        axis.text.y = element_text(color="black", 
                                   size=16, angle=0),
        axis.title=element_text(size= 20),
        legend.position = c(0.7, 0.8),
        legend.text = element_text(size=18)) +
        scale_color_manual(values=c("dodgerblue2", "black"))

#spec

#ggsave(paste0("/Volumes/Donnees/switchdrive/Mehdi-Celine/Chapitre - Theories and models of urbanization/simulation_HHI_index.pdf"), spec,
#       device = "pdf", width = 20, height = 20, units = "cm")


plot(cos(2*pi/V * 1:V), sin(2*pi/V * 1:V), cex = sqrt(inc)/12, xlim = c(-1.5,1.5),
     ylim = c(-1.5,1.5), xaxt='n', yaxt='n', xlab = "", ylab = "", pch = 16)
inc
sum(discovered)
unique(firm_ind[active])
sum(inc > 0)
length(unique(firm_ind[active]))

plot(entry_time[active], apply(profit_prec, 1, sum)[active], col = firm_ind[active])

plot(entry_time[active], tech[active], col = firm_city[active])

plot(1:sum(active), sort(apply(profit_prec, 1, sum)[active], decreasing = TRUE), log= "xy")

