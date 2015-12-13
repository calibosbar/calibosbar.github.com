



library(ggplot2)
color_mag_sersic[c("U.r","Mr")]
ggplot(color_mag_sersic, aes(x=Mr, y=U.r, colour=sersic))+geom_point(size=4)+ scale_colour_gradientn(colours=rainbow(4))



ggplot(color_mag_sersic, aes(x=Log.Mass, y=U.r, colour=morph))+geom_point(size=5)


lista_agn[c("U.r","Mr")]
ggplot(lista_agn, aes(x=Mr, y=U.r, colour=morph))+geom_point()
ggplot(lista_agn, aes(x=Log.Mass, y=U.r, colour=morph))+geom_point()

