#In Desktop/manuscript folder 
# Assuming your data is stored in a data frame called 'my_data'
# Filter out missing values
file_path <- "RFDist.xlsx"  # Update with the correct path
data <- read_excel(file_path)
#my_data <- my_data[complete.cases(my_data),]
library(readxl)
library(ggplot2)
library(tidyverse)
library(reshape2)
# Transform the data to long format using melt
data_long <- melt(data, id.vars = "Completeness", variable.name = "Method", value.name = "Value")
# my_colors <- c("#66C2A5", "#99D8C9", "#E78AC3", "#DADAEB", "#d3d3d3")
my_colors <- c("#99D8C9", "#66C2A5", "#A6CEE3", "#1F78B4", "#E78AC3")
dark_colors <- c("#238B45", "#238B45", "#004C8A", "#004C8A", "#AA3377")  # Dark Green, Dark Green, Dark Purple, Dark Purple, Grey

p <- ggplot(data_long, aes(x = as.factor(Completeness), y = Value, fill = Method)) +
  geom_boxplot(width = 0.95,  # Increase width of the boxes
               position = position_dodge2(preserve = "single", padding = 0.15),  # Adjust dodge width and spacing
               aes(color = Method)) +  # Remove outliers and color by Method
  scale_y_reverse(limits = c(0.9, 0), breaks = seq(0, 0.9, by = 0.1), labels = seq(0, 0.9, by = 0.1)) +  # Reverse y-axis and set breaks and labels
  scale_x_discrete(expand = expansion(add = c(0.2, 0.2)), labels = c("20%", "40%", "60%", "80%", "99%")) +  # Adjust space between categories and set custom labels
  theme_minimal() +
  scale_fill_manual(values = my_colors, name = "Sampling Type and Software Tool") +  # Apply custom colors and set legend title
  scale_color_manual(values = dark_colors, name = "Sampling Type and Software Tool") +  # Set manual colors for boxplot outlines and legend
  labs(x = "Backbone Tree Completeness", y = "Robinson-Foulds Distance") +  # Adjust axis labels
  scale_x_discrete(labels = c("20%", "40%", "60%", "80%", "99%", expression(atop("~20%", "(biased sampling)")))) +  # Set custom labels with line break
  # Adjust theme elements
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Set axis line color
    axis.ticks = element_line(color = "black"),  # Set tick color
    axis.ticks.length = unit(5, "pt"),  # Set tick length
    panel.border = element_blank(),  # Remove panel border
    plot.margin = margin(10, 10, 10, 10),  # Reduce margin around the plot
    legend.key = element_rect(color = NA, fill = "white"),  # Remove border color around legend keys, keep white fill
    legend.title = element_text(color = "black"),  # Set legend title color
    axis.text = element_text(size = 18),  # Increase axis labels text size
    axis.title = element_text(size = 22),  # Increase axis titles text size
    legend.position = "none"  
  )


p 


#---------------------------------------------------------
#Correct Placement Percentage
file_path <- "CorrectPlacementPercentageNew.xlsx"  # Update with the correct path
data <- read_excel(file_path)
#my_data <- my_data[complete.cases(my_data),]
library(readxl)
library(ggplot2)
library(tidyverse)
library(reshape2)
# Transform the data to long format using melt
data_long <- melt(data, id.vars = "Completeness", variable.name = "Method", value.name = "Value")

# Define custom color palettes
# my_colors <- c("#99cfe0", "#ffe680", "#b3e6b3", "#ffb366", "#d3d3d3")  # Lighter shades of blue, yellow, green, orange, grey
# dark_colors <- c("#00688B", "#B8860B", "#006400", "#8B4513", "#696969")  # Darker shades of blue, yellow, green, orange, grey

my_colors <- c("#99D8C9", "#66C2A5", "#A6CEE3", "#1F78B4", "#E78AC3")
dark_colors <- c("#238B45", "#238B45", "#004C8A", "#004C8A", "#AA3377")  # Dark Green, Dark Green, Dark Purple, Dark Purple, Grey

p <- ggplot(data_long, aes(x = as.factor(Completeness), y = Value, fill = Method)) +
  geom_boxplot(width = 0.95,size=0.6, position = position_dodge2(preserve = "single", padding = 0.15), aes(color = Method)) +  # Adjust width and dodge position, set outline color by Method
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +  # Set y-axis limits and breaks
  scale_x_discrete(expand = expansion(add = c(1, 0.05)), labels = c("20%", "40%", "60%", "80%", "99%")) +  # Adjust space between categories and set custom labels
  theme_minimal() +
  scale_fill_manual(values = my_colors, name = "Sampling Type and Software Tool") +  # Apply custom colors and set legend title
  scale_color_manual(values = dark_colors, name = "Sampling Type and Software Tool") +  # Set manual colors for boxplot outlines and legend
  labs(x = "Backbone Completeness", y = "Correct Placement Percentage (%)") +  # Adjust axis labels
  scale_x_discrete(labels = c("20%", "40%", "60%", "80%", "99%", expression(atop("~20%", "(biased sampling)")))) +  # Set custom labels with line break
  coord_cartesian(ylim = c(0, 100), expand = TRUE) +  # Set y-axis limits and prevent expansion for y-axis only
  # Adjust theme elements
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Set axis line color
    axis.ticks = element_line(color = "black"),  # Set tick color
    axis.ticks.length = unit(5, "pt"),  # Set tick length
    panel.border = element_blank(),  # Remove panel border
    plot.margin = margin(20, 20, 20, 20),  # Add margin around the plot
    legend.key = element_rect(color = NA, fill = "white"),  # Remove border color around legend keys, keep white fill
    legend.title = element_text(color = "black"),  # Set legend title color
    axis.text = element_text(size = 18),  # Increase axis labels text size
    axis.title = element_text(size = 22),  # Increase axis titles text size
    legend.position = "none" 
  ) +
  guides(fill = guide_legend(override.aes = list(color = dark_colors)))  # Adjust legend keys and boxplot outlines to match dark_colors
p

#------------

#Correct Genus Percentage
file_path <- "CGP.xlsx"  # Update with the correct path
data <- read_excel(file_path)
#my_data <- my_data[complete.cases(my_data),]
library(readxl)
library(ggplot2)
library(tidyverse)
library(reshape2)
# Transform the data to long format using melt
data_long <- melt(data, id.vars = "Completeness", variable.name = "Method", value.name = "Value")
# Define custom color palettes
my_colors <- c("#99D8C9","#66C2A5", "#A6CEE3", "#1F78B4")
dark_colors <- c("#238B45", "#238B45", "#004C8A", "#004C8A")  # Dark Green, Dark Green,blue, blue


p <- ggplot(data_long, aes(x = as.factor(Completeness), y = Value, fill = Method)) +
  geom_boxplot(width = 0.9,size=0.6, position = position_dodge2(width = 1.2,preserve = "single", padding = 0.19), aes(color = Method)) +  # Adjust width and dodge position, set outline color by Method
  scale_y_continuous(limits = c(0, 100), breaks = seq(0, 100, by = 20)) +  # Set y-axis limits and breaks
  scale_x_discrete(expand = expansion(add = c(0.5, 1)), labels = c("20%", "40%", "60%", "80%", "99%")) +  # Adjust space between categories and set custom labels
  theme_minimal() +
  scale_fill_manual(values = my_colors, name = "Sampling Type and Software Tool") +  # Apply custom colors and set legend title
  scale_color_manual(values = dark_colors, name = "Sampling Type and Software Tool") +  # Set manual colors for boxplot outlines and legend
  labs(x = "Backbone Completeness", y = "Correct Genus Placement Percentage (%)") +  # Adjust axis labels
  coord_cartesian(ylim = c(20, 100), expand = TRUE) +  # Set y-axis limits and prevent expansion for y-axis only
  # Adjust theme elements
  theme(
    panel.grid.major = element_blank(),  # Remove major grid lines
    panel.grid.minor = element_blank(),  # Remove minor grid lines
    axis.line = element_line(color = "black"),  # Set axis line color
    axis.ticks = element_line(color = "black"),  # Set tick color
    axis.ticks.length = unit(5, "pt"),  # Set tick length
    panel.border = element_blank(),  # Remove panel border
    plot.margin = margin(20, 20, 20, 20),  # Add margin around the plot
    legend.key = element_rect(color = NA, fill = "white"),  # Remove border color around legend keys, keep white fill
    legend.title = element_text(color = "black"),  # Set legend title color
    axis.text = element_text(size = 18),  # Increase axis labels text size
    axis.title = element_text(size = 22),  # Increase axis titles text size
    legend.position = "none" 
  ) 
p

