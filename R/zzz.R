.onLoad <- function(libname, pkgname) {
  devtools::load_all('/net/mraid14/export/tgdata/users/davidbr/src/dbutils/')
  doMC::registerDoMC(60)
  doFuture:: registerDoFuture()
  future::plan(future::multiprocess)
  
  ggplot2::theme_set(ggplot2::theme_light() %+replace% 
                     theme(panel.background = element_blank(), 
                           panel.grid.minor = element_blank(),
                           axis.text.x = element_text(angle = 90, hjust = 1),
                           strip.text = element_text(size = 6, colour = "white", face = 'bold')))
  options(dplyr.width = 200)
  options(future.globals.maxSize=100E9)
  options(width=200)
}

# plan(batchtools_sge(resources = list(queue = "all.q", threads = 1, memory = 1), workers = Inf) )