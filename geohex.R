gh_normalize <- function(lon) {
  if (lon < -180) {
    return(lon + 360) 
  } else {
      if lon > 180 {
      return(lon - 360)
    }
  }
  return(lon)
}