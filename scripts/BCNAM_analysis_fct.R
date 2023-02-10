############################
# BCNAM analysis functions #
############################


#### DD function ####

DD <- function(tmin, tmax, Tbase = 8, Tmax = 45){
  
  n_days <- length(tmin)
  DD <- rep(NA, n_days)
  for(d in 1:n_days){
    tav <- (tmin[d] + tmax[d])/2
    if(tav < Tmax){
      DD[d] <- max(tav - Tbase, 0) 
    } else DD[d] <- 0
    
  }
  
  DD
  
}

#### function to extract EC values from list of points ####

# grid_list: list of data frame that contain a number of environmental characteristics
# for each point.

# EC: character that specify which EC (column) of the individual matrix should
# be used.

# stat: the type of statistic: mean or sum. Other could be added if needed

# year: the year that need to be selected

# start: the starting point of the selected period in DOY

# end: the end point of the selected period in DOY

# grid_list = met_list_ext
# EC = 'T2M'
# stat = 'mean'
# year = '2012'
# start = 172
# end = start + 120
# m_i <- met_list_ext[[1]]

get_EC_val <- function(EC, stat = 'mean', year = 2012, grid_list, start = 172, end = 292){
  
  n_pt <- length(grid_list)
  EC_val <- rep(NA, n_pt)
  
  for(i in 1:n_pt){
    
    m_i <- grid_list[[i]]
    m_i <- m_i[m_i$YEAR == year, ]
    EC_vect <- m_i[start:end, EC]
    
    if(stat == 'sum'){
      m_EC <- mean(EC_vect, na.rm = TRUE)
      EC_vect[is.na(EC_vect)] <- m_EC
      EC_val[i] <- sum(EC_vect)
    } else {
      EC_val[i] <- mean(EC_vect, na.rm = TRUE)
    }
    
  }
  
  return(EC_val)
  
}


#### function to process data.frame before plotting ####

add_xy_and_dim_id <- function(d_p){
  
  d_p <- data.frame(d_p)
  
  # pop_env lk
  pop_env_un <- unique(as.character(d_p[, 1]))
  pop_env_lk <- length(pop_env_un):1
  names(pop_env_lk) <- pop_env_un
  
  # par lk
  p_un <- unique(as.character(d_p[, 2]))
  p_lk <- 1:length(p_un)
  names(p_lk) <- p_un
  
  d_p$y <- pop_env_lk[as.character(d_p[, 1])]
  d_p$x <- p_lk[as.character(d_p[, 2])]
  
  
  # Reduce the parent names to plot
  p_un <- as.character(p_un)
  p_un[p_un == 'IS15401'] <- 'IS15'
  p_un[p_un == 'IS23540'] <- 'IS23'
  p_un[p_un == 'SC56614'] <- 'SC56'
  
  return(list(d_p = d_p, pop_env_un = pop_env_un,  p_un = p_un))
  
}

#### plot Qmain QxE ####

plot_Qmain_QxE <- function(d_pl, pop_env_un, par_un, var_unit, main, text.size = 12){
  
  # define the bolded text.y [could be removed later if problem with ggplot2]
  face.y <- rep('plain', length(pop_env_un))
  face.y[grep(pattern = '_main', x = pop_env_un)] <- 'bold'
  
  p <- ggplot(d_pl) +
    geom_point(aes(x=x, y=y, fill = effect, size = log10pval), shape = 21) +
    scale_fill_gradient2(low = "red", mid = "white", high = "blue",
                         guide = guide_colorbar(order = 1), name = var_unit) +
    scale_y_continuous(breaks = 1:length(pop_env_un), labels = rev(pop_env_un)) +
    scale_x_continuous(breaks = 1:length(par_un), labels = as.character(par_un),
                       position = "top") +
    xlab("Parents") + ylab("Population x Environment") +
    ggtitle(main) + 
    theme(axis.title.x = element_text(size=text.size),
          axis.title.y = element_text(size=text.size),
          axis.text.y = element_text(size = text.size-2, face = rev(face.y)),
          axis.text.x = element_text(size = text.size-4, angle = 45),
          plot.title = element_text(size=(text.size+2)),
          strip.text.x =  element_text(size=text.size),
          legend.title = element_text(size=(text.size)),
          legend.text = element_text(size=(text.size)))
  
  return(p)
  
}

#### Reduce QpxEC matrix given significance ####

red_tab <- function(dt, d_cond, thre){
  
  dt[d_cond < thre] <- NA
  dt <- dt[, which(!colSums(is.na(dt)) == nrow(dt)), drop = FALSE]
  dt <- dt[which(!rowSums(is.na(dt)) == ncol(dt)), , drop = FALSE]
  
  return(dt)
}

#### Format EC matrix for projection  ####

pop_EC_par_mat_form <- function(d_B, d_main, d_sign_pt, v_sel, meta_inf){
  
  pop_vect <- sapply(X = meta_inf[v_sel], FUN = `[[`, 1)
  par_vect <- sapply(X = meta_inf[v_sel], FUN = `[[`, 2)
  
  d_m <- data.table(pop = pop_vect, par = par_vect, d_main[v_sel, , drop = FALSE])
  d_B <- data.table(pop = pop_vect, par = par_vect, d_B[v_sel, , drop = FALSE])
  d_s <- data.table(pop = pop_vect, par = par_vect, d_sign[v_sel, , drop = FALSE])
  
  d_m <- melt(data = d_m, id.vars = c('pop', 'par'))
  d_B <- melt(data = d_B, id.vars = c('pop', 'par'))
  d_s <- melt(data = d_s, id.vars = c('pop', 'par'))
  
  d_m <- d_m[!is.na(value)]
  d_B <- d_B[!is.na(value)]
  d_s <- d_s[!is.na(value)]
  
  d_EC <- data.table(d_m, B = d_B$value, sign = d_s$value)
  colnames(d_EC)[3:4] <- c('EC', 'int')
  
  return(d_EC) 
  
}

#### calculate y_pred = int + B_EC * EC ####

EC_proj <- function(d_EC, EC_mat, lat_lon){
  
  d <- c()
  
  for(i in 1:nrow(d_EC)){
    
    EC_i <- as.character(d_EC$EC[i])
    y_pred <- d_EC$int[i] + d_EC$B[i] * EC_mat[, EC_i]
    d_i <- data.frame(pop = d_EC$pop[i], par = d_EC$par[i], EC = EC_i,
                      int = d_EC$int[i], B = d_EC$B[i], lat = lat_lon$lat,
                      long = lat_lon$long, y_pred, logpval = d_EC$sign[i])
    
    d <- rbind(d, d_i)
    
  }
  
  # order according to most sign EC and parent
  EC_freq <- sort(tapply(d$EC, d$EC, length), decreasing = TRUE)
  EC_sign <- d %>% group_by(EC) %>% summarise(sign = mean(logpval))
  EC_tab <- data.frame(EC_sign, freq = EC_freq[EC_sign$EC])
  EC_tab <- EC_tab %>% arrange(desc(freq), desc(sign))
  
  EC_lev <- EC_tab$EC
  par_lev <- c("EC", names(sort(table(d$par), decreasing = TRUE)))
  
  y_range <- range(d$y_pred)
  av_sign <- mean(d$logpval, na.rm = TRUE)
  
  # Add the raw EC values
  
  for(j in 1:ncol(EC_mat)){
    
    # scale the value
    y <- EC_mat[, j]
    EC_range <- range(y)
    Beta <- (y_range[2] - y_range[1])/(EC_range[2] - EC_range[1])
    y_bar <- y_range[1] + (y - EC_range[1]) * Beta
    
    d_j <- data.frame(pop = d_EC$pop[1], par = 'EC', EC = colnames(EC_mat)[j],
                      int = NA, B = NA, lat = lat_lon$lat, long = lat_lon$long,
                      y_pred = y_bar, logpval = av_sign)
    
    d <- rbind(d, d_j)
    
  }
  
  d$EC <- factor(d$EC, levels = EC_lev)
  d$par <- factor(d$par, levels = par_lev)
  
  d <- data.table(d)
  
  return(d)
  
}

#### plot QTLxEC projection Mali ####

plot_QpxEC <- function(d_p, EC_range, Mali_layer, size_txt = 3, main, var_id){
  
  # text information
  d_text <- d_p[, .(B = round(mean(B), 3)), by = .(par, EC)]
  d_text <- d_text[!(par == 'EC')]
  d_text$B <- paste0('B = ', d_text$B)
  
  EC_range$EC <- factor(rownames(EC_range), levels = levels(d_p$EC))
  EC_range$par = factor('EC', levels = levels(d_p$par))
  EC_range$rg <- paste0('DAS: ', EC_range$start, '-', EC_range$end)
  
  # plot QpxEC
  QEC_plot <- Mali_layer +
    geom_point(data=d_p, shape = 20, aes(long, lat, color = y_pred, size = logpval)) +
    coord_equal() +
    scale_colour_gradient(low = "#2166AC", high = "#B2182B", name = var_id) +
    facet_grid(rows = vars(EC), cols = vars(par), switch = "y") +
    geom_text(data = EC_range, size = size_txt,
              mapping = aes(x = -7.5, y = 14.7, label = rg)) +
    geom_text(data = d_text, size = size_txt,
              mapping = aes(x = -7.5, y = 14.7, label = B)) +
    ggtitle(main)
  
  return(QEC_plot)
  
}