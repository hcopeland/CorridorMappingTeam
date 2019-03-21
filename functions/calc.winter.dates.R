# This function takes the metadata file and calculates the start and end of the winter period
# based on quantiles of the end of fall migration (for start of winter)
# and the starts of spring migration (for end of winter).
# written by Jerod Merkle, 7 Feb 2019 

calc.winter.dates <- function(mig.metadata.file="C:/Users/jmerkle/Desktop/Mapp2/tab6output/metadata.csv",
                              qtl.end.fall.mig=0.95, qtl.start.spring.mig=0.05){
  
  #manage packages
  if(all("stringr" %in% installed.packages()[,1])==FALSE)
    stop("You must install the following packages: stringr")
  require(stringr)
  
  fl <- read.csv(mig.metadata.file)
  fl$seas <- substr(str_split_fixed(fl$input.file, "_",2)[,2],1,2)
  fl$Start.Date <- as.POSIXct(strptime(paste0("2018", substr(fl$Start.Date,5,19)),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  fl$End.Date <- as.POSIXct(strptime(paste0("2018", substr(fl$End.Date,5,19)),format = "%Y-%m-%d %H:%M:%S"), tz="GMT")
  
  
  sp <- fl[fl$seas == "sp",]
  fa <- fl[fl$seas == "fa",]

  wint.start <- round(quantile(fa$End.Date, probs = qtl.end.fall.mig),0)
  wint.end <- round(quantile(sp$Start.Date, probs = qtl.start.spring.mig),0)
  
  return(data.frame(winter_start=paste(strftime(wint.start, format = "%m", tz = "GMT"),strftime(wint.start, format = "%d", tz = "GMT"),sep="-"),
              winter_end=paste(strftime(wint.end, format = "%m", tz = "GMT"),strftime(wint.end, format = "%d", tz = "GMT"),sep="-")))

}#end function

