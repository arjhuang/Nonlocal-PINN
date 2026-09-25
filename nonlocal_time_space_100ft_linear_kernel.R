### LAX–Friedrichs Scheme (R version)
### Import Density Data
density_data <- as.matrix( read.csv("data/NGSIM_US101_Density_Data.csv", header = FALSE) )

################# 
dx <- 20 
dt_in <- 5
ttot <- 2695
N <- as.integer(ttot / dt_in) # number of coarse time steps 
x_left <- 0
x_right <- 2060

# right_end <- 40
right_end <- 100

gamma <- 0.02
xtot <- x_right - x_left
K <- as.integer(xtot / dx) # 2060 / 20 = 103 
etatot <- as.integer(right_end / dx) # 60 / 20 = 3

x <- seq(x_left, x_right, length.out = K + 1) 
ti <- seq(0, ttot, length.out = N + 1) 
############# Kernel definition
d <- right_end
dx <- 20
s <- seq(0, d, by = dx)

# unnormalized linear kernel: K(s) = d - s
K_raw <- d - s

# unnormalized exponen kernel1 K(s) = exp(-s/d)
# K_raw <- exp(-s/d)

# unnormalized exponen kernel2 K(s) = exp(-s/d) - exp(-1)
# K_raw <- (exp(-s / d)-exp(-1))

## unnormalized exponen kernel 3 K(s) = exp(- 1/(s-d)^2)
# K_raw <- (exp(-1/(s-d)^2))


# trapezoidal rule function 
trapz <- function(x, y) 
{ sum((y[-1] + y[-length(y)]) * diff(x) / 2) }

# compute area of unnormalized kernel 
area_raw <- trapz(s, K_raw) 

# normalized kernel so integral = 1 
eta_n <- K_raw / area_raw 
# show results 
x_collar <- x_right - right_end
xtot_coll <- x_collar - x_left 
xk <- xtot_coll / dx 
### Max velocity and density 
vm <- 80 
um <- 0.12
### Time refinement (robust)
time_scale <- 100
dt <- dt_in / time_scale # refined dt
N_b <- as.integer(N * time_scale + 1) # refined time steps

t_preserve <- gamma * right_end # seconds
t_ind <- as.integer(round(t_preserve / dt)) 
lag_sec <- gamma * dx 
lag_step <- as.integer(round(lag_sec / dt))
u1 <- matrix(0, nrow = K + 1, ncol = N_b)

### Initial condition that is thick

i_preserve <- t_ind + 1

for (i in 1:i_preserve)
{
  k1 <- (i - 1) %% time_scale
  k2 <- (i - 1 - k1) / time_scale 
  if (k1 == 0) 
  { 
    interp_col <- density_data[, k2 + 1]
  }
  else 
  {
    interp_col <- (density_data[, k2 + 1] * (time_scale - k1) + k1 * density_data[, k2 + 2]) / time_scale
  } 
  u1[2:xk, i] <- interp_col[2:xk]
}





### Left boundary condition
for (i in 1:(N_b - 1))
{ 
  k1 <- (i - 1) %% time_scale
  k2 <- (i - 1 - k1) / time_scale 
  if (k1 == 0) 
  { 
    u1[1, i] <- density_data[1, k2 + 1] 
  }
  else 
  {
    u1[1, i] <- (density_data[1, k2 + 1] * (time_scale - k1) + k1 * density_data[1, k2 + 2]) / time_scale 
  } 
}
u1[1, N_b] <- density_data[1, k2 + 2]


### Right boundary condition
for (j in (xk + 1):(K + 1)) 
{
  for (i in 1:(N_b - 1)) 
  {
    k1 <- (i - 1) %% time_scale
    k2 <- (i - 1 - k1) / time_scale
    if (k1 == 0) 
    {
      u1[j, i] <- density_data[j, k2 + 1] 
    }
    else 
    {
      u1[j, i] <- (density_data[j, k2 + 1] * (time_scale - k1) + k1 * density_data[j, k2 + 2]) / time_scale 
    } 
  }
  u1[j, N_b] <- density_data[j, k2 + 2] 
}

### Main time loop
#lag_sec <- gamma * dx
#lag_step <- lag_sec/dt



for (i in i_preserve:(N_b - 1))
{
  #### doing the convolution 
  u_d <- u1[, i]
  for (ud_index in 1:xk+1)
  {
    acc <- 0 
    for (k in 0:(etatot - 1))
    { 
      space_idx <- ud_index + k
      time_idx <- i - k * lag_step # t, t-gamma*dx, t-2*gamma*dx...
      acc <- acc + u1[space_idx, time_idx] * eta_n[k + 1]  
    } 
    u_d[ud_index] <- acc * dx
  }
  ### Flux 
  f <- vm * u1[, i] * (1 - u_d / um)
  ### Space loop 
  for (j in 2:xk)
  {
    u1[j, i + 1] <- 0.5 * (u1[j - 1, i] + u1[j + 1, i]) - (dt / (2 * dx)) * (f[j + 1] - f[j - 1]) 
  }
}




### Downsampling back to original time grid
index <- 1:540
scaled_up_index <- (index - 1) * time_scale + 1
compressed_u <- matrix(1, nrow = K + 1, ncol = N + 1)
for (j in 1:(K + 1)) 
{ 
  compressed_u[j, ] <- u1[j, scaled_up_index]
}


### Error computation
U <- compressed_u[1:(xk + 1), 1:(N + 1)]
D <- density_data[1:(xk + 1), 1:(N + 1)] 
S <- sum((U - D)^2) 
S_e <- sum(D^2) 
error <- 100 * S / S_e
print(error)

compressed_u <- t(compressed_u)
density_data <- t(density_data)

# open a new device with larger left margin
par(mar = c(5, 6.5, 4, 2))  # bottom, left, top, right

### Plotting Density
filled.contour(ti, x, compressed_u, color.palette = terrain.colors,  xlab = "",
               ylab = "",
               key.axes = axis(4, cex.axis = 1.5),
               
               plot.axes = {
                 axis(1, cex.axis = 1.5)
                 axis(2, cex.axis = 1.5)
                 
                 mtext("time (s)", side = 1, line = 2.5, cex = 1.5)
                 mtext("x (feet)", side = 2, line = 3.8, cex = 1.5, las = 3)
               }
               
               
)