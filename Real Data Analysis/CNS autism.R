#Componentwise Nodewise Selection
#Dataset ASD
data.path<-""
setwd(data.path)

# abide dataset: autism group
load("time.series.ABIDE.RData")
n1<-dim(h.1)[1]
p1<-dim(h.1)[2]
t1<-dim(h.1)[3]

n2<-dim(h.2)[1]
p2<-dim(h.2)[2]
t2<-dim(h.2)[3]

autism<- h.1
control<-h.2
library(fda)
library(mvtnorm)
library(glmnet)
library(igraph)
library(pROC)
set.seed(2025)
#autism group
###########################################
# Functional representation with b-spline #
###########################################
#L=6, lambda=lambda.1se, lasso

L <- 6 #numero di basi
t_seq <- seq(0, 1, length.out = t1)

basis <- create.bspline.basis(rangeval = c(0,1), nbasis = L)

# b-spline coefficients
coef_array <- array(0, dim = c(n1, p1, L))
for(i in 1:n1){
  for(j in 1:p1){
    fd_obj <- smooth.basis(t_seq, autism[i,j,], basis)$fd
    coef_array[i,j,] <- fd_obj$coefs
  }
}

# scaling
dim(coef_array)
coef_array_scaled <- array(0, dim = dim(coef_array))
for (j in 1:dim(coef_array)[2]) { 
  for (k in 1:dim(coef_array)[3]) { 
    coef_array_scaled[, j, k] <- scale(coef_array[, j, k])
  }
}

adj_est <- matrix(0, p1, p1)

for(j in 1:p1){
  Y_j <- coef_array_scaled[, j, ]
  
  X_j <- matrix(0, nrow = n1, ncol = (p1-1)*L)
  idx <- 1
  for(k in 1:p1){
    if(k != j){
      X_j[, idx:(idx + L - 1)] <- coef_array_scaled[, k, ]
      idx <- idx + L
    }
  }
  
  fit <- cv.glmnet(X_j, Y_j, family = "mgaussian", alpha = 1)
  
  lambda_stronger <- fit$lambda.1se * 3
  coef_list <- coef(fit, s = lambda_stronger)
  
  idx <- 1
  for(k in 1:p1){
    if(k != j){
      start_col <- idx
      end_col <- idx + L - 1
      beta_kj <- sapply(coef_list, function(mat) mat[start_col:end_col, 1])  # L x L
      if(any(beta_kj != 0)){
        adj_est[j, k] <- 1
      }
      idx <- idx + L
    }
  }
}

adj_est_symm_autism <- (adj_est & t(adj_est)) * 1  # AND

layout <- function(list, option = c("xy","xz","yz")) {
  p <- length(list)
  result <- matrix(NA, nrow = p, ncol = 2)
  coord.vec <- unlist(list)
  x.vec <- coord.vec[seq(1, length(coord.vec), 3)]
  y.vec <- coord.vec[seq(2, length(coord.vec), 3)]
  z.vec <- coord.vec[seq(3, length(coord.vec), 3)]
  
  if(option == "xy") {
    result[,1] <- x.vec
    result[,2] <- y.vec
  } else if(option == "xz") {
    result[,1] <- x.vec
    result[,2] <- z.vec
  } else if(option == "yz") {
    result[,1] <- y.vec
    result[,2] <- z.vec
  }
  
  return(result)
}

#graph.plotter <- function(adj.mat, coordinate.list, name.tag.vec, 
#                          option = "yz", title = "Graph", color = "black") {

#  p <- ncol(adj.mat)
#  nbr.count <- colSums(adj.mat) 

#  net <- graph_from_adjacency_matrix(adj.mat, mode = "undirected") %>%
#    set_vertex_attr("label", value = name.tag.vec)

#  layout_orig <- layout(coordinate.list, option)


#  plot(net, vertex.label = name.tag.vec, vertex.shape = "circle",
#       vertex.color = "red", vertex.frame.color = "red",
#       vertex.label.font = 1, vertex.label.family = "Helvetica",
#       vertex.label.cex = 0.3, vertex.label.dist = 0.4, vertex.label.degree = pi/4,
#       edge.color = color, edge.width = 1.2, vertex.size = 0.5 * (nbr.count + 1),
#       layout = layout_orig, edge.curved = 0.8)
#}
#thesis layout
graph.plotter <- function(adj.mat, coordinate.list, name.tag.vec, 
                          option = "yz", title = "Graph", color = "black") {
  
  p <- ncol(adj.mat)
  nbr.count <- colSums(adj.mat) 
  
  net <- graph_from_adjacency_matrix(adj.mat, mode = "undirected") %>%
    set_vertex_attr("label", value = name.tag.vec)
  
  layout_orig <- layout(coordinate.list, option)
  
  plot(net,
       vertex.label = name.tag.vec,
       vertex.shape = "circle",
       vertex.color = "#FF1493",
       vertex.frame.color = "black", 
       vertex.label.font = 2,
       vertex.label.family = "Helvetica",
       vertex.label.cex = 1,  
       vertex.label.dist = 0.8, 
       vertex.label.degree = 0,   
       vertex.label.color = "black",
       edge.color = "black",
       edge.width = 1.5,
       vertex.size = 1.2 * (nbr.count + 1),  
       layout = layout_orig,
       edge.curved = FALSE)
  
}




adj.mat <- adj_est_symm_autism
file.path <- ""
graph.path <- paste(file.path, "Graph", sep="/")

func.path <- ""
save.path <- ""
runtime.path <- ""

source(paste(file.path,"AAL.info.R",sep=""))
coordinate.list <- coordinate.list 

name.tag.vec <- name.tag.vec

graph.plotter(adj.mat, coordinate.list, name.tag.vec, "yz", "My Graph")
num_edges_s <- sum(adj_est_symm_autism) / 2
cat("Numero di archi nel grafo stimato:", num_edges_s, "\n")


# 2% sparsity - Zhao
setwd(file.path)
load("G.spa.ctrl.autism.RData")
graph.plotter(G.autism.sym, coordinate.list, name.tag.vec, "yz", "My Graph")

num_edges_s <- sum(G.autism.sym) / 2
cat("n. edges zhao:", num_edges_s, "\n")


#Dataset CONTROL

############################################
# Functional representation with b-spline #
############################################
#L=6, lambda=lambda.1se, lasso

L <- 6
t_seq <- seq(0, 1, length.out = t2)

basis <- create.bspline.basis(rangeval = c(0,1), nbasis = L)

# B-spline coefficients
coef_array <- array(0, dim = c(n2, p2, L))
for(i in 1:n2){
  for(j in 1:p2){
    fd_obj <- smooth.basis(t_seq, control[i,j,], basis)$fd
    coef_array[i,j,] <- fd_obj$coefs
  }
}

# scaling
dim(coef_array)
coef_array_scaled <- array(0, dim = dim(coef_array))
for (j in 1:dim(coef_array)[2]) { 
  for (k in 1:dim(coef_array)[3]) { 
    coef_array_scaled[, j, k] <- scale(coef_array[, j, k])
  }
}

adj_est <- matrix(0, p2, p2)

for(j in 1:p1){
  Y_j <- coef_array_scaled[, j, ]
  
  X_j <- matrix(0, nrow = n2, ncol = (p2-1)*L)
  idx <- 1
  for(k in 1:p2){
    if(k != j){
      X_j[, idx:(idx + L - 1)] <- coef_array_scaled[, k, ]
      idx <- idx + L
    }
  }
  
  fit <- cv.glmnet(X_j, Y_j, family = "mgaussian", alpha = 1)
  
  lambda_stronger <- fit$lambda.1se * 3.7
  coef_list <- coef(fit, s = lambda_stronger)
  
  idx <- 1
  for(k in 1:p2){
    if(k != j){
      start_col <- idx
      end_col <- idx + L - 1
      beta_kj <- sapply(coef_list, function(mat) mat[start_col:end_col, 1])  # L x L
      if(any(beta_kj != 0)){
        adj_est[j, k] <- 1
      }
      idx <- idx + L
    }
  }
}

adj_est_symm_control <- (adj_est & t(adj_est)) * 1  # AND


adj.mat <- adj_est_symm_control
file.path <- ""
graph.path <- paste(file.path, "Graph", sep="/")

func.path <- ""
save.path <- ""
runtime.path <- ""

source(paste(file.path,"AAL.info.R",sep=""))
coordinate.list <- coordinate.list 

name.tag.vec <- name.tag.vec

graph.plotter(adj.mat, coordinate.list, name.tag.vec, "yz", "My Graph")
num_edges_s <- sum(adj_est_symm_control) / 2
cat("n. edges estimated graph:", num_edges_s, "\n")


#2% sparsity - Zhao
setwd(file.path)
load("G.spa.ctrl.control.RData")
graph.plotter(G.control.sym, coordinate.list, name.tag.vec, "yz", "My Graph")

num_edges_s <- sum(G.control.sym) / 2
cat("n. edges Zhao:", num_edges_s, "\n")





##############################
# AAL codes
codes <- c(2001, 2002, 2101, 2102, 2111, 2112, 2201, 2202, 2211, 2212,
           2301, 2302, 2311, 2312, 2321, 2322, 2331, 2332, 2401, 2402,
           2501, 2502, 2601, 2602, 2611, 2612, 2701, 2702, 3001, 3002,
           4001, 4002, 4011, 4012, 4021, 4022, 4101, 4102, 4111, 4112,
           4201, 4202, 5001, 5002, 5011, 5012, 5021, 5022, 5101, 5102,
           5201, 5202, 5301, 5302, 5401, 5402, 6001, 6002, 6101, 6102,
           6201, 6202, 6211, 6212, 6221, 6222, 6301, 6302, 6401, 6402,
           7001, 7002, 7011, 7012, 7021, 7022, 7101, 7102, 8101, 8102,
           8111, 8112, 8121, 8122, 8201, 8202, 8211, 8212, 8301, 8302,
           9001, 9002, 9011, 9012, 9021, 9022, 9031, 9032, 9041, 9042,
           9051, 9052, 9061, 9062, 9071, 9072, 9081, 9082,
           9100, 9110, 9120, 9130, 9140, 9150, 9160, 9170)

names <- c("Precentral_L", "Precentral_R", "Frontal_Sup_L", "Frontal_Sup_R",
           "Frontal_Sup_Orb_L", "Frontal_Sup_Orb_R", "Frontal_Mid_L", "Frontal_Mid_R",
           "Frontal_Mid_Orb_L", "Frontal_Mid_Orb_R", "Frontal_Inf_Oper_L", "Frontal_Inf_Oper_R",
           "Frontal_Inf_Tri_L", "Frontal_Inf_Tri_R", "Frontal_Inf_Orb_L", "Frontal_Inf_Orb_R",
           "Rolandic_Oper_L", "Rolandic_Oper_R", "Supp_Motor_Area_L", "Supp_Motor_Area_R",
           "Olfactory_L", "Olfactory_R", "Frontal_Sup_Medial_L", "Frontal_Sup_Medial_R",
           "Frontal_Med_Orb_L", "Frontal_Med_Orb_R", "Rectus_L", "Rectus_R",
           "Insula_L", "Insula_R", "Cingulum_Ant_L", "Cingulum_Ant_R",
           "Cingulum_Mid_L", "Cingulum_Mid_R", "Cingulum_Post_L", "Cingulum_Post_R",
           "Hippocampus_L", "Hippocampus_R", "ParaHippocampal_L", "ParaHippocampal_R",
           "Amygdala_L", "Amygdala_R", "Calcarine_L", "Calcarine_R",
           "Cuneus_L", "Cuneus_R", "Lingual_L", "Lingual_R",
           "Occipital_Sup_L", "Occipital_Sup_R", "Occipital_Mid_L", "Occipital_Mid_R",
           "Occipital_Inf_L", "Occipital_Inf_R", "Fusiform_L", "Fusiform_R",
           "Postcentral_L", "Postcentral_R", "Parietal_Sup_L", "Parietal_Sup_R",
           "Parietal_Inf_L", "Parietal_Inf_R", "SupraMarginal_L", "SupraMarginal_R",
           "Angular_L", "Angular_R", "Precuneus_L", "Precuneus_R",
           "Paracentral_Lobule_L", "Paracentral_Lobule_R", "Caudate_L", "Caudate_R",
           "Putamen_L", "Putamen_R", "Pallidum_L", "Pallidum_R",
           "Thalamus_L", "Thalamus_R", "Heschl_L", "Heschl_R",
           "Temporal_Sup_L", "Temporal_Sup_R", "Temporal_Pole_Sup_L", "Temporal_Pole_Sup_R",
           "Temporal_Mid_L", "Temporal_Mid_R", "Temporal_Pole_Mid_L", "Temporal_Pole_Mid_R",
           "Temporal_Inf_L", "Temporal_Inf_R", "Cerebelum_Crus1_L", "Cerebelum_Crus1_R",
           "Cerebelum_Crus2_L", "Cerebelum_Crus2_R", "Cerebelum_3_L", "Cerebelum_3_R",
           "Cerebelum_4_5_L", "Cerebelum_4_5_R", "Cerebelum_6_L", "Cerebelum_6_R",
           "Cerebelum_7b_L", "Cerebelum_7b_R", "Cerebelum_8_L", "Cerebelum_8_R",
           "Cerebelum_9_L", "Cerebelum_9_R", "Cerebelum_10_L", "Cerebelum_10_R",
           "Vermis_1_2", "Vermis_3", "Vermis_4_5", "Vermis_6",
           "Vermis_7", "Vermis_8", "Vermis_9", "Vermis_10")

get_area <- function(nm) {
  if (grepl("^Frontal|^Precentral|^Rectus|^Olfactory|Rolandic|Supp_Motor", nm)) return("frontal lobe")
  if (grepl("^Parietal|Postcentral|Precuneus|Angular|SupraMarginal|Paracentral", nm)) return("parietal")
  if (grepl("^Temporal|Heschl", nm)) return("temporal")
  if (grepl("^Occipital|Calcarine|Cuneus|Lingual|Fusiform", nm)) return("occipital")
  return("others")
}

# Dictionary: area<- color
node_area_dict <- setNames(sapply(names, get_area), codes)

area_colors <- c(
  "frontal lobe" = "#00C853", 
  "parietal"     = "#2962FF",  
  "temporal"     = "#FF1493",
  "occipital"    = "#FF6D00",  
  "others"       = "#9E9E9E"  
)



# dictionary: code<-color
node_color_dict <- setNames(area_colors[ node_area_dict ], names(node_area_dict))
graph.plotter <- function(adj.mat,
                          coordinate.list,
                          name.tag.vec,
                          node_color_dict,
                          option = "yz") {
  
  net <- graph_from_adjacency_matrix(adj.mat, mode = "undirected") %>%
    set_vertex_attr("label", value = name.tag.vec)
  
  layout_orig <- layout(coordinate.list, option)
  vertex_cols <- node_color_dict[ as.character(name.tag.vec) ]
  nbr.count   <- colSums(adj.mat)
  
  plot(net,
       main                = title,
       vertex.label        = name.tag.vec,
       vertex.shape        = "circle",
       vertex.color        = vertex_cols,
       vertex.frame.color  = "black",
       vertex.label.font   = 2,
       vertex.label.family = "Helvetica",
       vertex.label.cex    = 1,
       vertex.label.dist   = 0.8,
       vertex.label.degree = 45 * pi/180,
       vertex.label.color  = "black",
       edge.color          = "black",
       edge.width          = 1.5,
       vertex.size         = 1.2 * (nbr.count + 1),
       layout              = layout_orig,
       edge.curved         = FALSE)
  
}

graph.plotter(adj.mat         = adj.mat,
              coordinate.list = coordinate.list,
              name.tag.vec    = codes,
              node_color_dict = node_color_dict)

# 2% sparaity-zhao
setwd(file.path)
load("G.spa.ctrl.autism.RData")
graph.plotter(adj.mat         = G.autism.sym, 
              coordinate.list = coordinate.list,
              name.tag.vec    = codes,
              node_color_dict = node_color_dict)


load("G.spa.ctrl.control.RData")
graph.plotter(adj.mat         = G.control.sym,
              coordinate.list = coordinate.list,
              name.tag.vec    = codes,
              node_color_dict = node_color_dict)



plot_difference_graph <- function(G1, G2, coordinate.list, name.tag.vec, node_color_dict, option = "yz") {
  only_G1 <- (G1 == 1) & (G2 == 0)
  only_G2 <- (G1 == 0) & (G2 == 1)
  both    <- (G1 == 1) & (G2 == 1)

  adj.diff <- matrix(0, nrow = nrow(G1), ncol = ncol(G1))
  adj.diff[only_G1] <- 1
  adj.diff[only_G2] <- 2
  adj.diff[both]    <- 3

  net <- graph_from_adjacency_matrix((adj.diff > 0)*1, mode = "undirected")
  layout_orig <- layout(coordinate.list, option)
  vertex_cols <- node_color_dict[ as.character(name.tag.vec) ]
  
  edge_types <- E(net)$weight <- apply(as_edgelist(net), 1, function(x) adj.diff[as.numeric(x[1]), as.numeric(x[2])])
  edge_cols <- ifelse(edge_types == 1, "#2962FF", 
                      ifelse(edge_types == 2, "#FF1744",     
                             "#000000"))                     
  
  plot(net,
       layout              = layout_orig,
       vertex.label        = name.tag.vec,
       vertex.color        = vertex_cols,
       vertex.label.cex    = 1,
       vertex.size         = 10, 
       vertex.label.color  = "black",
       edge.color          = edge_cols,
       edge.width          = 3,
       edge.curved         = FALSE,
       vertex.frame.color  = "black",
       vertex.label.font   = 2)
}

plot_difference_graph(adj.mat, G.autism.sym, coordinate.list, codes, node_color_dict)
plot_difference_graph(adj.mat, G.control.sym, coordinate.list, codes, node_color_dict)
