#' MultiKParallel main algorithm
#'
#' MultiKParallel main algorithm: takes a preprocessed gene expression matrix as input. Then subsamples 80\% of the cells and applies standard Seurat pipeline on the subsampled data matrix 100 times over multiple resolution parameters.
#' @param seu A Seurat object with normalized count
#' @param resolution A vector Seurat resolution parameters. Default is from 0.05 to 2 with step size of 0.05
#' @param nPC Number of principal components to use in clustering
#' @param reps Number of subsampling runs. Integer value. Default is 100
#' @param pSample Proportion of cells to sample. Numerical value. Default is 0.8
#' @param seed Optional numerical value. This sets a random seed for generating reproducible results
#' @return A list with components: k is a vector of number of runs for each K. clusters is a list containing the clustering labels for each subsampling run at each resolution parameter. consensus is a list containing a consensus matrix for each K.
#' @export
MultiKParallel <- function(seu, resolution = seq(0.05, 2, 0.05), nPC = 30, reps = 100, pSample = 0.8, seed = NULL, numCores = NULL) {
        # setting seed for reproducibility
        if (is.null(seed) == TRUE) {
                seed <- timeSeed <- as.numeric(Sys.time())
        }
        set.seed(seed)
        
        # setting up doParallel
        if (is.null(numCores)) {
                numCores <- parallel::detectCores() / 2
        }
        cl <- parallel::makeCluster(numCores)
        doSNOW::registerDoSNOW(cl)
        #doParallel::registerDoParallel(numCores)
        print(paste0("doSNOW: numCores = ", numCores, "."))
        
        # step 1: subsampling
        subcol <- list()
        for (i in 1: reps) {
                subcol[[i]] <- sample(x=ncol(seu), size=round(ncol(seu) * pSample), replace=FALSE)
        }
        
        # step 2: loop over subsampling runs, with each run subsampling 80% of cells, reselect genes for clustering
        print("Starting subsampling and clustering...")
        
        # setting up progress bar
        progress <- function(n) cat(sprintf("-- rep %d is complete\n", n))
        opts <- list(progress=progress)
        
        results <- foreach::foreach (i = 1:reps, .options.snow=opts) %dopar% {

                print(paste("Rep: ", i, sep=""))
                # subsample the columns (the cells) from the full matrix
                subX <- seu[, subcol[[i]] ]
                
                # normalizing the data
                #subX <- NormalizeData(object = subX, normalization.method = "LogNormalize", scale.factor = 10000, verbose=F)
                
                # Find HVG genes ~ 2000 genes
                subX <- Seurat::FindVariableFeatures(object = subX, selection.method = "vst", 
                                             nfeatures = 2000, loess.span = 0.3, 
                                             clip.max = "auto", num.bin = 20, 
                                             binning.method = "equal_width", verbose = F)
                
                # Scaling unwanted variation
                all.genes <- rownames(x = subX)
                subX <- Seurat::ScaleData(object = subX, features = all.genes, verbose = F)
                # Run PCA to reduce dimensions
                subX <- Seurat::RunPCA(object = subX, features = Seurat::VariableFeatures(object = subX), npcs = nPC, verbose=F)
                # Run Clustering
                subX <- Seurat::FindNeighbors(object = subX, k.param = 20, reduction = "pca", dims = 1:nPC, verbose = F)
                
                clusters <- list()
                messages <- c()
                ks <- c()
                
                for (res in resolution) {
                        print(paste("Rep", i, "Res", res, sep=" "))
                        subX <- Seurat::FindClusters(subX, resolution = res, verbose = F)
                        subX.clusters <- Seurat::Idents(subX)
                        message <- paste("Rep_", i, "_res_", res, sep = "")
                        messages <- c(messages, message)
                        clusters[[message]] <- subX.clusters
                        ks <- c(ks, length(unique(subX.clusters)))
                }

                # send results
                results_temp <- list("messages" = messages, "clusters" = clusters, "ks" = ks)
                results_temp
                
        }
        
        results_combined <- rep_combine(results)
        messages <- results_combined$messages
        clusters <- results_combined$clusters
        ks <- results_combined$ks
        
        print("DONE!")
        print("Starting consensus...")
        
        # step 3: calculate consensus matrix across subsampling runs for each unique K
        mInit <- matrix(0, ncol = ncol(seu), nrow = ncol(seu))
        unique.ks <- unique(ks)[order(unique(ks))]
        
        # setting up progress bar
        progress <- function(n) cat(sprintf("-- k-value %d is complete\n", n))
        opts <- list(progress=progress)
        
        res <- foreach::foreach (count.k = 1:length(unique.ks), .options.snow=opts) %dopar% {
                
                k = unique.ks[count.k]
                print(paste("k =", k, sep=" "))
                idx <- which(ks == k)
                cluster.k <- clusters[idx]

                for (s in 1: length(cluster.k) ) {
                        print(paste("run", s, sep = ""))
                        currClusters = names(cluster.k[[s]])
                        sampleKey <- sapply(currClusters, function(x, srtObj){which(colnames(srtObj) == x)}, srtObj=seu)
                        sampleKey <- as.numeric(sampleKey)
                        if (s == 1){
                                ml <- connectivityMatrix(cluster.k[[s]], mInit, sampleKey)
                                m.count <- connectivityMatrix(rep(1, length(sampleKey)), mInit, sampleKey)
                        } else {
                                ml <- connectivityMatrix(cluster.k[[s]], ml, sampleKey)
                                m.count <- connectivityMatrix(rep(1, length(sampleKey)), m.count, sampleKey)
                        }
                }
        
                res_k <- triangle(ml, mode = 3) / triangle(m.count, mode = 3)
                res_k[which(triangle(m.count, mode = 3) == 0)] = 0
                print(paste(k, " finished", sep = ""))
                
                # return results
                res_k
        }
        
        print("ALL DONE!")
        
        names(res) = 1:length(res)
        return(list("consensus" = res, "k" = ks))
        
}



rep_combine <- function(ret_list) {
        
        messages <- c()
        clusters <- list()
        ks <- c()
        for (i in 1:length(ret_list)) {
                rep_results <- ret_list[[i]]
                messages <- c(messages, rep_results$messages)
                clusters <- append(clusters, rep_results$clusters)
                ks <- c(ks, rep_results$ks)
          
        }
        
        combined_ret_list = list("messages" = messages, "clusters" = clusters, "ks" = ks)
        
        return(combined_ret_list)
        
}