library(ggplot2)
library(dplyr)
library(tidyr)
setwd("D:/postgraduate.research/my subject/Mendelian randomization/figure/circular pattern")
# Load your data
data <- read.csv('circular pattern91.csv')

# Ensure data has both columns Abbreviation and FullName
if (!("Abbreviation" %in% colnames(data)) || !("FullName" %in% colnames(data))) {
  stop("数据表格必须包含Abbreviation和FullName列")
}

# Add row numbers as id column
data$id <- seq_len(nrow(data))

# Calculate angles
num_points <- nrow(data)
data <- data %>%
  mutate(angle = 2 * pi * (id - 1) / num_points,
         angle_degrees = angle * 180 / pi)  # Calculate angles in degrees

# Calculate coordinates
data <- data %>%
  mutate(x = 0.225 * cos(angle),  # Adjust this value to change the radius of the blue circle
         y = 0.225 * sin(angle),
         text_angle = ifelse(angle_degrees > 90 & angle_degrees < 270, angle_degrees + 180, angle_degrees))  # Adjust text angles

# Plot circular pattern with reduced radius
p <- ggplot(data, aes(x = x, y = y)) +
  geom_point(color = "#1E90FF", size = 3) +
  geom_text(aes(label = Abbreviation, 
                x = 0.2 * cos(angle),  # Reduced radius for inside circle
                y = 0.2 * sin(angle), 
                angle = text_angle),  # Adjusted angle for inward radiating
            hjust = ifelse(data$angle_degrees > 90 & data$angle_degrees < 270, 0, 1), 
            vjust = 0.5, 
            size = 3) +  # Abbreviation inside circle, radiating inward
  geom_text(aes(label = FullName, 
                x = 0.25 * cos(angle),  # Reduced radius for outside circle
                y = 0.25 * sin(angle), 
                angle = text_angle), 
            hjust = ifelse(data$angle_degrees > 90 & data$angle_degrees < 270, 1, 0), 
            vjust = 0.5, 
            size = 3) +  # FullName outside circle, radiating outward
  coord_fixed() +
  theme_void() +
  labs(title = "环形文本图", x = "", y = "") +
  theme(legend.position = "none") +
  theme(plot.title = element_text(size = 20),  # Adjust title size
        text = element_text(size = 12)) +  # Adjust text size
  xlim(-1, 1) +  # Adjust the limits of x-axis
  ylim(-1, 1)  # Adjust the limits of y-axis

# Save as PDF
ggsave("circular_pattern.pdf", p, width = 25, height = 25, units = "in")
