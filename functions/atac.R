subset_bed_file = function(x, intervals=NULL) {

  d = 
    rtracklayer::import(
      con = x,
      which = intervals, 
      format="bed"
    )
  
  strand(d) = d$name
  start(d) = start(d) - 1
  d$name = NULL
  
  return(d)
}
