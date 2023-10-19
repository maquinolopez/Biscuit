BISCUIT_Logo_GitHub_Hex <- function() {
  par(mar=c(2,2,2,2))
  
  # Create an empty square plot
  plot(1, 1, xlim=c(0, 100), ylim=c(0, 100),
       xlab="", ylab="", type="n", xaxt='n', yaxt='n', asp=1)
  
  # Adjusted Hexagon coordinates for a bigger hexagon
  hex_x <- c(5, 95, 100, 95, 5, 0)
  hex_y <- c(100, 100, 50, 0, 0, 50)
  
  # Set the clipping region to the hexagon
  clip(min(hex_x), max(hex_x), min(hex_y), max(hex_y))
  polygon(hex_x, hex_y, col=rgb(0.9, 0.9, 0.9, 0.2), lwd=2, border=NA)
  
  # Drawing DNA double helices as sequences
  x_vals <- seq(0, 100, length.out = 200)
  y1_vals <- 12 * sin(0.2 * x_vals) + 57
  y2_vals <- 12 * cos(0.2 * x_vals) + 57
  
  mat_gradient <- colorRampPalette(c("darkorchid4", "mediumpurple1"))(length(x_vals))
  for (i in 1:(length(x_vals)-1)) {
    polygon(c(x_vals[i], x_vals[i+1], x_vals[i+1], x_vals[i]), 
            c(y1_vals[i], y1_vals[i+1], y2_vals[i+1], y2_vals[i]), 
            col=mat_gradient[i], border=NA)
  }
  
  lines(x_vals, y1_vals, col="darkorchid4", lwd=3)
  lines(x_vals, y2_vals, col="darkorchid4", lwd=3)
  
  # Representing tie-points at peaks and troughs
  tie_points_x <- seq(10, 90, by=20)
  points(tie_points_x, 12 * sin(0.2 * tie_points_x) + 57, col="darkred", pch=19, cex=2.5)
  points(tie_points_x, 12 * cos(0.2 * tie_points_x) + 57, col="darkred", pch=19, cex=2.5)
  
  # Drawing "BISCUIT" above the sequences
  text(50, 82, "BISCUIT", cex = 2.5, font = 2, col = "black")
  
  # Adjusted font size for better visibility of the bottom legend
  text(50, 18, "Bayesian Integration of Sequences", cex = 1.5, font = 3, col = "grey30", adj = c(0.5, 0.5))
  text(50, 10, "Using Chronology Intersection Tie-points", cex = 1.5, font = 3, col = "grey30", adj = c(0.5, 0.5))
}

# Run the function
BISCUIT_Logo_GitHub_Hex()

# Save the plot as a PNG
png(filename = "~/GitHub/BISCUIT_Logo.png", width = 400, height = 400)
BISCUIT_Logo_GitHub_Hex_v4()
dev.off()
