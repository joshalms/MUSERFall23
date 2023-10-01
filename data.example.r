#### Create Data for Example

library(tidyverse)

wide <- data.frame (Team  = c("A", "B", "C", "D"),
                  Points = c(88,91,99,94),
                  Assists = c(12,17,24,28),
                  Rebounds = c(22,28,30,31.)
)

long <- data.frame(Team = c("A", "A" , "A" , "B" , "B" , "B" , "C" , "C" , "C", "D", "D", "D"),
                  Variable = rep(c("Points", "Assists", "Rebounds"), 4),
                  Value = c(88,12,22,91,17,28,99,24,30,94,28,31))



ggplot(wide,
       aes(x = Points , y = Assists, color = Team)) +
       geom_point()


ggplot(long,
       aes(x = Value , y = Variable, color = Team)) +
  geom_point()
