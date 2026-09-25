### LAX-FRIEDRICHS SCHEME: LOCAL VERSION

### Import Density Data
density_data <- as.matrix(read.csv("data/NGSIM_US101_Density_Data.csv", header = FALSE))

#################
dx <- 20
dt_in <- 5
ttot <- 2695

N <- as.integer(ttot / dt_in)
x_left <- 0
x_right <- 2060
xtot <- x_right - x_left
K <- as.integer(xtot / dx)

x <- seq(x_left, x_right, length.out = K + 1)
ti <- seq(0, ttot, length.out = N + 1)

### Max velocity and density
vm <- 80
um <- 0.12

### Time refinement
time_scale <- 100
dt <- dt_in / time_scale
N_b <- as.integer(N * time_scale + 1)

u1 <- matrix(0, nrow = K + 1, ncol = N_b)

### Initial condition
u1[, 1] <- density_data[, 1]

### Left boundary condition
for (i in 1:(N_b - 1)) {
  k1 <- (i - 1) %% time_scale
  k2 <- (i - 1 - k1) / time_scale
  
  if (k1 == 0) {
    u1[1, i] <- density_data[1, k2 + 1]
  } else {
    u1[1, i] <- (density_data[1, k2 + 1] * (time_scale - k1) +
                   k1 * density_data[1, k2 + 2]) / time_scale
  }
}
u1[1, N_b] <- density_data[1, k2 + 2]

### Right boundary condition
for (i in 1:(N_b - 1)) {
  k1 <- (i - 1) %% time_scale
  k2 <- (i - 1 - k1) / time_scale
  
  if (k1 == 0) {
    u1[K + 1, i] <- density_data[K + 1, k2 + 1]
  } else {
    u1[K + 1, i] <- (density_data[K + 1, k2 + 1] * (time_scale - k1) +
                       k1 * density_data[K + 1, k2 + 2]) / time_scale
  }
}
u1[K + 1, N_b] <- density_data[K + 1, k2 + 2]

### Main time loop: LOCAL flux only
for (i in 1:(N_b - 1)) {
  
  ### Local flux: f(u) = vm * u * (1 - u/um)
  f <- vm * u1[, i] * (1 - u1[, i] / um)
  
  ### Space loop
  for (j in 2:K) {
    u1[j, i + 1] <- 0.5 * (u1[j - 1, i] + u1[j + 1, i]) -
      (dt / (2 * dx)) * (f[j + 1] - f[j - 1])
  }
}

### Downsampling back to original time grid
scaled_up_index <- (0:N) * time_scale + 1
compressed_u <- matrix(0, nrow = K + 1, ncol = N + 1)

for (j in 1:(K + 1)) {
  compressed_u[j, ] <- u1[j, scaled_up_index]
}

### Error computation
U <- compressed_u[1:(K + 1), 1:(N + 1)]
D <- density_data[1:(K + 1), 1:(N + 1)]
S <- sum((U - D)^2)
S_e <- sum(D^2)
error <- 100 * S / S_e
print(error)

### Transpose for plotting
compressed_u_plot <- t(compressed_u)
density_data_plot <- t(density_data)

# open a new device with larger left margin
par(mar = c(5, 6.5, 4, 2))  # bottom, left, top, right

filled.contour(
  ti, x, compressed_u_plot,
  color.palette = terrain.colors,
  xlab = "",
  ylab = "",
  key.axes = axis(4, cex.axis = 1.5),
  
  plot.axes = {
    axis(1, cex.axis = 1.5)
    axis(2, cex.axis = 1.5)
    
    mtext("time (s)", side = 1, line = 2.5, cex = 1.5)
    mtext("x (feet)", side = 2, line = 3.8, cex = 1.5, las = 3)
  }
  
  
)

