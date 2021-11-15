
library(titanic)
library(corrplot)
library(magrittr)
library(tidyverse) 

df = titanic_train %>%
  as_tibble %>%
  select(Survived, Class = Pclass, Sex, Age, 
         TicketID = Ticket, TicketPrice = Fare, Port = Embarked) %>%
  mutate(Class = as.character(Class)) %>%
  drop_na %>%
  filter(Port!="") %>%
  mutate(across(where(function(col) is.character(col) & length(unique(col))=="2"),
                ~ .x %>% factor(labels = c(0, 1)) %>% as.character %>% as.numeric ) )
  
k = ncol(df)

n = ncol(df)

cardinality_threshold = 0.00

coef_deter = function(x) {
  
  x %>%
    summary %>%
    extract2("adj.r.squared") %>%
    sqrt %>% 
    replace_na(0) %>%
    return()
  
}

score_mat = matrix(data = NA, nrow = k, ncol = k)

for (i in 1:k) { # loop over targets
  
  rowv = pull(df, i)
  
  for (j in 1:k) { # loop over features
    
    colv = pull(df, j)
    
    if (i==j) {
      score_mat[i, j] = 1
      next
      
    } else if (is.character(rowv)) {
      
      flev = (table(rowv)/n) 
      
      flev = flev[flev > cardinality_threshold]
      
      if (is_empty(flev)) {
        
        score_mat[i, j] = 0
        next
        
      }
      
      vecrsq = vector(length = length(flev))
      
      for (l in 1:length(flev)) {
        
        rowv_level = if_else(rowv==names(flev)[l], 1, 0)
        
        if (is.character(colv)) {
          
          vecrsq[l] =
            lm(formula = rowv_level ~ colv) %>%
            coef_deter
          
        } else {
          
          vecrsq[l] =
            lm(formula = rowv_level ~ colv + I(colv^2)) %>%
            coef_deter
          
        }
        
      }
      
      score_mat[i, j] = weighted.mean(x = vecrsq,
                                      w = flev)
      
    } else {
      
      if (is.character(colv)) {
        
        score_mat[i, j] =
          lm(formula = rowv ~ colv) %>%
          coef_deter
        
      } else {
        
        score_mat[i, j] =
          lm(formula = rowv ~ colv + I(colv^2)) %>%
          coef_deter
        
      }
      
    }
    
  }
  
}

windows(); corrplot(score_mat %>%
                set_colnames(names(df)) %>%
                set_rownames(names(df)),
                    is.corr = FALSE,
              method = "color",
              addCoef.col = "red",
              tl.col = "black",
              tl.pos = "lt")
gridGraphics::grid.echo() # Convert a scene that was drawn using the graphics package to an identical scene drawn with the grid package.
corplot = grid::grid.grab() 
grid::grid.text("Correlation matrix plot", 
          x = unit(0.5, "npc"), 
          y = unit(1, "npc"), 
          gp = grid::gpar(fontface = "bold",
                    fontsize = 10))
ggsave(corplot, file="mult_corr_plot.png")

windows(); list(c(1, 0.36, 0.59, 0.19, 0.33, 0.39, 0.24),
     c(0.2, 1, 0, 0.23, 0.51, 0.9, 0.17),
     c(0.57, 0, 1, 0.085, 0.24, 0.27, 0.042),
     c(0,0.069, 0, 1, 0, 0,0),
     c(0.0069, 0.019, 0.0084, 0.013, 1, 0.67, 0.029),
     c(0, 0.2, 0, 0, 0.64, 1, 0),
     c(0, 0, 0, 0.054, 0.45, 0.74, 1)) %>%
  reduce(rbind) %>%
  set_colnames(names(df)) %>%
  set_rownames(names(df)) %>%
  corrplot(is.corr = FALSE,
           method = "color",
           addCoef.col = "red",
           tl.col = "black",
           tl.pos = "lt")
gridGraphics::grid.echo() # Convert a scene that was drawn using the graphics package to an identical scene drawn with the grid package.
corplot = grid::grid.grab() 
grid::grid.text("Correlation matrix plot", 
                x = unit(0.5, "npc"), 
                y = unit(1, "npc"), 
                gp = grid::gpar(fontface = "bold",
                                fontsize = 10))
ggsave(corplot, file="PPS_plot.png")
