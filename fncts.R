create_circle <- function(DET_RAD){
  if(is.integer(DET_RAD/2)){
    ss <- 2*pi/(DET_RAD*12)
  } else {
    ss <- 2*pi/(DET_RAD*2*(5+(DET_RAD-2)))
  }
  
  lapply(seq(ss, 2*pi, ss), function(a){
    return(c(x=cos(a)*DET_RAD,
             y=sin(a)*DET_RAD) %>% round())
  }) %>% unique() %>% bind_rows() %>%
    return()
}

get_single_index <- function(x,y,nr){
  
  nr*(x-1)+y
  
}

get_xy_index <- function(i, nr){
  
  y <- i%%nr
  if(y==0){
    y_ret <- nr
    x <- (i-y)/nr 
  } else {
    y_ret <- y
    x <- (i-y)/nr+1
  }
  
  return(c(x=x, y=y_ret))
  
}

dist_pts <- function(x1, y1, x2, y2){
  
  
  sqrt((x1-x2)^2+(y1-y2)^2)
  
}

measure_mn <- function(sel, radius, pixel_shift, channel_dir, output_dir){
  print(sel)
  df_pos_raw <- bioimagetools::readTIF(sel)[,,1,1] %>%
    as_tibble()
  nrow_img <- nrow(df_pos_raw)
  df_pos <- df_pos_raw %>%
    rownames_to_column("y") %>%
    gather(x,i, -y) %>%
    
    mutate(x=str_replace_all(x, "V", "") %>% as.numeric(),
           y=as.numeric(y)) 
  
  whites <- df_pos %>%
    .[which(.$i==1),] %>%
    #filter(i==1) %>%
    select(x,y) %>%
    mutate(c=fpc::dbscan(., eps = sqrt(2), MinPts = 1)$cluster)
  
  center_group <- whites %>%
    group_by(c) %>%
    tally() %>%
    filter(n==min(.$n)) %>%
    pull(c)
  
  center <- whites %>%
    filter(c==center_group) %>%
    select(x,y) %>%
    kmeans(centers = 1) %>%
    .$centers
  
  print("center detected")
  
  border <- whites %>%
    filter(!c==center_group) %>%
    select(x,y) %>%
    as.matrix()
  
  
  mask <- apply(border, 1, function(px){
    #print(px)
    full_circle <- lapply(seq_len(radius), create_circle) %>%
      #Reduce(function(x,y)rbind(x,y),.)
      
      bind_rows() %>%
      mutate(x=x+px[["x"]],
             y=y+px[["y"]])
    
    
    # return(cbind(full_circle$x+px[["x"]],
    #             full_circle$y+px[["y"]]))
  }) %>% 
    #Reduce(function(x,y)rbind(x,y),.) %>%
    bind_rows() %>%
    unique()
  
  
  mask_ind <- mask  %>%
    rowwise() %>%
    mutate(dist_to_center=dist_pts(x,y,center[1], center[2])) %>%
    filter(dist_to_center<=cutoff_px) %>%
    mutate(ind=get_single_index(x,y,nrow_img)) %>%
    pull(ind)
  
  
  print("mask created")
  
  border_length <- length(border)
  
  ## load channels 
  ## red channel
  ch00 <- readTIF(file.path(channel_dir, str_replace(basename(sel), "mn.+.tif", "ch00.tif"))) #%>%
  
  intensity_sum <- ch00 %>%
    .[,,1,1] %>%
    .[mask_ind] %>%
    sum()
  
  show <- ch00 %>%
    .[,,1,1]%>%
    as_tibble() %>%
    rownames_to_column("y") %>%
    gather(x,i, -y) %>%
    
    mutate(x=str_replace_all(x, "V", "") %>% as.numeric(),
           y=as.numeric(y))  %>%
    rowwise() %>%
    mutate(ind=get_single_index(x,y,nrow_img))
  
  print("channel loaded")
  
  rel <- show[which(show$ind %in% mask_ind),]
  
  cutoff <- EBImage::otsu(matrix(rel$i, nrow=1))
  
  rel2 <- rel %>%
    mutate(i=ifelse(i>cutoff, 1, 0)) %>%
    filter(i==1)
  
  rel3 <- rel2 %>%
    select(x,y) %>% 
    ungroup() %>%
    mutate(c=fpc::dbscan(., eps = sqrt(2), MinPts = 9)$cluster) %>%
    .[which(.$c!=0),]
  
  cluster_size= rel3 %>% group_by(c) %>% tally() %>%
    rowwise() %>%
    mutate(
      #n_synapse_raw=n/32,
      n_synapse_raw=n/(2*median(.$n)),
      n_synapse=ifelse(n_synapse_raw<1, 1, floor(n_synapse_raw)))
  
  counted <- sum(cluster_size$n_synapse)    
  
  print("synapses counted")
  ### plots
  
  p_raw <- ggplot(rel, aes(x=x, y=-y, fill=i))+
    geom_tile(show.legend = F)+
    scale_fill_gradient2(mid="black", high="red")+
    geom_point(data=as_tibble(center), inherit.aes=F, 
               aes(x=x, y=-y),shape=16, size=5)+
    geom_point(data=as_tibble(center), inherit.aes=F, 
               aes(x=x, y=-y),shape=4, size=10)+
    theme_minimal()+
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())  
  p_otsu <- ggplot(rel2, aes(x=x, y=-y))+
    geom_tile(show.legend = F, fill="red")+
    geom_point(data=as_tibble(center), inherit.aes=F, 
               aes(x=x, y=-y),shape=16, size=5)+
    geom_point(data=as_tibble(center), inherit.aes=F, 
               aes(x=x, y=-y),shape=4, size=10)+
    theme_minimal()+
    theme(axis.text = element_blank(),
          axis.title = element_blank(),
          panel.grid = element_blank())  
  
  if(!sum(cluster_size$n_synapse)==nrow(cluster_size)){
    sec_data <- lapply(cluster_size %>% filter(n_synapse>1) %>% pull(c), function(C){
      
      relevant <- rel3 %>%
        .[which(.$c==C),] %>%
        kmeans(centers = 1) %>%
        .$centers %>%
        as_tibble() %>%
        return()
      
    }) %>% bind_rows() %>%
      left_join(cluster_size %>% select(c,n_synapse))
    
    
    
    p_dc <- ggplot(rel3 %>%
                     left_join(cluster_size %>% select(c,n_synapse)), 
                   aes(x=x, y=-y, fill=n_synapse))+
      scale_fill_gradientn(colours=c("grey", "yellow", "red2"))+
      geom_tile(show.legend = F, 
                #fill="red"
      )+
      
      geom_text(data=sec_data, inherit.aes=F,
                aes(x=x, y=-y, label=as.character(n_synapse)))+
      
      geom_point(data=as_tibble(center), inherit.aes=F, 
                 aes(x=x, y=-y),shape=16, size=5)+
      geom_point(data=as_tibble(center), inherit.aes=F, 
                 aes(x=x, y=-y),shape=4, size=10)+
      theme_minimal()+
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank())
    
    
    
  } else {
    p_dc <- ggplot(rel3 %>%
                     left_join(cluster_size %>% select(c,n_synapse)), 
                   aes(x=x, y=-y, fill=n_synapse))+
      scale_fill_gradientn(colours=c("grey", "yellow", "red2"))+
      geom_tile(show.legend = F, 
                #fill="red"
      )+
      
      
      geom_point(data=as_tibble(center), inherit.aes=F, 
                 aes(x=x, y=-y),shape=16, size=5)+
      geom_point(data=as_tibble(center), inherit.aes=F, 
                 aes(x=x, y=-y),shape=4, size=10)+
      theme_minimal()+
      theme(axis.text = element_blank(),
            axis.title = element_blank(),
            panel.grid = element_blank())
  }
  
  plot_name <- file.path(output_dir,
                         basename(sel) %>% str_replace(".tif$", ".pdf"))
  pdf(file=plot_name,
      width=0.03*(3*abs(min(mask$x)-max(mask$x))),
      height=0.03*(abs(min(mask$y)-max(mask$y))))
  gridExtra::grid.arrange(p_raw, p_otsu, p_dc, nrow=1)
  dev.off()
  
  
  plot_name <- file.path(output_dir,
                         basename(sel) %>% str_replace(".tif$", ".png"))
  png(file=plot_name,
      width=2*(3*abs(min(mask$x)-max(mask$x))),
      height=2*(abs(min(mask$y)-max(mask$y))))
  gridExtra::grid.arrange(p_raw, p_otsu, p_dc, nrow=1)
  dev.off()  
  
  
  return(tibble(name=basename(sel),
                n_px=length(mask_ind),
                n_px_over_otsu=nrow(rel2),
                n_synapse=counted,
                int_sum=intensity_sum,
                estimated_border_length=border_length))  
  
}

