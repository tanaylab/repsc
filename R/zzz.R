.onLoad <- function(libname, pkgname) {
  doMC::registerDoMC(60)
  doFuture:: registerDoFuture()
  future::plan(future::multiprocess(workers = future::availableCores()))
  
  ggplot2::theme_set(ggplot2::theme_light() %+replace% 
                     theme(panel.background = element_blank(), 
                           panel.grid.minor = element_blank(),
                           panel.spacing = unit(0, "lines"),
                           axis.text.x = element_text(angle = 90, hjust = 1),
                           strip.text = element_text(size = 6, colour = "white", face = 'bold')))
  options(dplyr.width = 200)
  options(tgs_use.blas = TRUE)
  options(future.globals.maxSize=100E9)
  options(width=200)
}

# plan(batchtools_sge(resources = list(queue = "all.q", threads = 1, memory = 1), workers = Inf) )
custom.col <- c("#FFDB6D", "#C4961A", "#F4EDCA", 
                "#D16103", "#C3D7A4", "#52854C", "#4E84C4", "#293352")