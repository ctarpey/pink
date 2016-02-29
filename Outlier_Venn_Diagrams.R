### Making Venn Diagrams for Outliers
###    
### Carolyn Tarpey | February 2016
### ---------------------------------------


#install.packages("VennDiagram")
library(VennDiagram)


draw.quad.venn()

venn.plot <- draw.quad.venn(
  area1 = 272,
  area2 = 502,
  area3 = 293,
  area4 = 755,
  n12 = 18,
  n13 = 25,
  n14 = 22,
  n23 = 35,
  n24 = 53,
  n34 = 22,
  n123 = 11,
  n124 = 2,
  n134 = 6,
  n234 = 6,
  n1234 = 2,
  category = c("PT1", "PT2", "Latitude", "Longitude"),
  fill = c("orange", "red", "green", "blue"),
  cex = 2,
  cat.cex = 2,
  cat.col = c("orange", "red", "green", "blue")
)