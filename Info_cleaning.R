Info <- read.csv("./ordered_strain_info.csv", header = T)
str(Info)
unique(Info$Substrate_of_isolation)
unique(Info$Geographic_location)

#Regroup substrate of isolation into categories
Info$Substrate_of_isolation <- gsub(".*(Wine|wine).*", "Wine", Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub("(Quercus|Quecus|Oak).*", "Oak", Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub(".*(Fermentation|fermentation|Fermented|Brewery|brewing|
                                    brewerie|beer|sake|Bioethanol|bioethanol|Champagne|ethanol).*", 
                                    "Other fermentation", Info$Substrate_of_isolation)
Info$Substrate_of_isolation <- gsub(".*(Fruit|fruit|Grape|fig|Apple|apple|Vineyard).*", "Fruit", Info$Substrate_of_isolation)



#Regroup geographic location by country 
Info$Geographic_location <- gsub(".*,\\s(\\w+$)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(Greece)", "\\1", Info$Geographic_location)
Info$Geographic_location <- gsub(".*(Portugal).*", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(Greece)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(South Africa)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(UK)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub(".*(France)", "\\1", Info$Geographic_location) 
Info$Geographic_location <- gsub("(Asia).*", "\\1", Info$Geographic_location)