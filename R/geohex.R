library("hash")
library("cwhmisc")

pkg.env <- new.env()

pkg.env$H_chars=c("A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P","Q","R","S","T","U","V","W","X","Y","Z","a","b","c","d","e","f","g","h","i","j","k","l","m","n","o","p","q","r","s","t","u","v","w","x","y","z")
pkg.env$H_base= 20037508.34
pkg.env$H_d2r   <- pi / 180
pkg.env$H_k <- tan(pi/6)
pkg.env$H_er    <- 6371007.2
pkg.env$H_index <- list()
pkg.env$H_units <- vector("list", 13)

for (i in 1:length(pkg.env$H_chars)) {
  val = pkg.env$H_chars[i]
  pkg.env$H_index[[val]] = i-1
}


gh_point_to_ne <- function(point) {
  northing = (pkg.env$H_k * point$x * point$unit$width + point$y * point$unit$height) / 2
  easting = (northing - point$y * point$unit$height) / pkg.env$H_k
  return(list('northing'=northing, 'easting'=easting))
}

gh_normalize <- function(lon) {
  if (lon < -180) {
    return(lon + 360) 
  } else {
    if (lon > 180) {
      return(lon - 360)
    }
  }
  return(lon)
}

gh_digitize <- function(code) {
  ind1 = substr(code, 1, 1)
  ind2 = substr(code, 2, 2)
  n1 = pkg.env$H_index[[ind1]]
  n2 = pkg.env$H_index[[ind2]]
  
  digits = toString(n1 * 30 + n2)

  if (nchar(code) >= 3) {
    digits = paste(digits, substring(code, 3), sep="")
  }
  
  return(digits) 
}

gh_easting <- function(lon) {
  return( gh_normalize(lon) * pkg.env$H_base / 180.0 )
}

gh_northing <- function(lat) {
  val = tan((90+lat) * pkg.env$H_d2r / 2)
  return( log(tan((90+lat) * pkg.env$H_d2r / 2)) / pi * pkg.env$H_base ) 
}

gh_point <- function(lat, lon, level) {
  u = gh_unit(level)
  e = gh_easting(lon)
  n = gh_northing(lat)
  x = (e + ( n / pkg.env$H_k ) ) / u$width
  y = (n - ( pkg.env$H_k * e ) ) / u$height
  
  x0 = floor(x)
  y0 = floor(y)
  
  xd = x - x0
  yd = y - y0
  
  if( (yd > -xd+1) && (yd < 2*xd) && (yd > 0.5*xd) ) {
    xn = x0 + 1
    yn = y0 + 1
  } else if ( (yd < -xd+1) && (yd > (2*xd-1)) && (yd < 0.5*xd + 0.5) ) {
    xn = x0
    yn = y0
  } else {
    xn = floor(x + 0.499999)
    yn = floor(y + 0.499999)
  }
  
  return( list( x=xn , y=yn, unit=u ))
}

gh_unit <- function(level) {
  
  if ( is.null(pkg.env$H_units[[level+1]]) ) {
    size   = pkg.env$H_base / 3^(level+3)
    scale  = size / pkg.env$H_er
    width  = 6 * size
    height = width * pkg.env$H_k
    h = list()
    h$level=level
    h$size=size
    h$width=width
    h$height=height
    h$scale=scale
    pkg.env$H_units[[level+1]] = h
  } 
  val = pkg.env$H_units[[level+1]]
  return(val)
}

gh_parse <- function(code) {
  x = 0
  y = 0
  len = nchar(code)
  digits = gh_digitize(code)
  
  pos = 0
  for (i in 1:nchar(digits)) {
    n10 = as.numeric(substr(digits, i, i))
    pow = 3^(len-pos)
    pos = pos + 1
    
    n3 = ""
    repeat {
      if (n10 <= 0) {
        break
      }
      n3 = paste( (n10 %% 3), n3, sep="")
      n10 = floor(n10 / 3)
    }
    
    n3 = as.numeric(n3)
    
    xd = floor(n3 / 10)
    if (xd == 0) {
      x = x - pow
    } else if (xd == 2) {
      x = x + pow
    }

    yd = floor(n3 %% 10)
    if (yd == 0) {
      y = y - pow
    } else if (yd == 2) {
      y = y + pow
    }
    
  }
  return( list(x=x, y=y, unit=gh_unit(len-2) ))
}

#' GeoHex V3 encode function
#'
#' This function takes a lat,lon, level and returns GeoHex V3 code
#' @param 
#' lat Latitude
#' lon Longitude
#' level GeoHex level
#' @keywords geohex
#' @export
#' @examples 
#' gh_encode(85.05112507763846, 89.37952995300293, 15)
#' 
gh_encode <- function(lat, lon, level){
  point = gh_point(lat, lon, level)
  ne = gh_point_to_ne(point)
  code = ""
  mod_x = 0
  mod_y = 0
  
  if ( (pkg.env$H_base - ne$easting) < point$unit$size) {
    mod_x = point$y
    mod_y = point$x
  } else {
    mod_x = point$x
    mod_y = point$y
  }
  
  level = point$unit$level+2
  for( i in level:0) {
    pow = 3^i
    p2c = ceiling(pow / 2)
    
    if (mod_x >= p2c) {
      mod_x = mod_x - pow
      c3_x = 2
    } else if (mod_x <= -p2c) {
      mod_x = mod_x + pow
      c3_x = 0
    } else {
      c3_x = 1
    }
    
    if (mod_y >= p2c) {
      mod_y = mod_y - pow
      c3_y = "2"
    } else if (mod_y <= -p2c) {
      mod_y = mod_y + pow
      c3_y = 0
    } else {
      c3_y = 1
    }
      
    c3xy = paste(c3_x, c3_y, sep="")
    code = paste(code, strtoi(c3xy, 3), sep="" )
    
  }
  
  num = as.numeric( substr(code, 1, 3) )
  first = floor(num/30) + 1
  second = floor(num%%30) + 1
  final = paste(pkg.env$H_chars[[first]], pkg.env$H_chars[[second]], substring(code, 4), sep="")
  return(final)
}

#' GeoHex V3 decode function
#'
#' This function takes a GH V3 code, then return lat,lon and level
#' @param code Valid GeoHex V3 string
#' @keywords geohex
#' @export
#' @examples 
#' gh_decode('bb3371844188')
#' 
gh_decode <- function(code) {
  point = gh_parse(code)
  
  ne = gh_point_to_ne(point)
  
  if ( pkg.env$H_base - ne$easting < point$unit$size ) {
    lon = 180
  } else {
    lon = gh_normalize(ne$easting / pkg.env$H_base * 180) 
  }
  
  lat = 180 / pi * ( 2 * atan( exp(ne$northing / pkg.env$H_base * 180 * pkg.env$H_d2r)) - pi / 2) 

  return(list(lat=lat, lon=lon, level=point$unit$level))
}