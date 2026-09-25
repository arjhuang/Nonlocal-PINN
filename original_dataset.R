### Import Density Data
density_data <- as.matrix(read.csv("data/NGSIM_US101_Density_Data.csv", header = FALSE))

### Parameters
dx <- 20
dt <- 5
ttot <- 2695
x_left <- 0
x_right <- 2060

### Grid sizes
K <- as.integer((x_right - x_left) / dx)
N <- as.integer(ttot / dt)

### Spatial and time axes
x <- seq(x_left, x_right, length.out = K + 1)
ti <- seq(0, ttot, length.out = N + 1)

### Transpose data for plotting (time rows, space columns)
density_data <- t(density_data)


# open a new device with larger left margin
par(mar = c(5, 6.5, 4, 2))  # bottom, left, top, right

filled.contour(
  ti,
  x,
  density_data,
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