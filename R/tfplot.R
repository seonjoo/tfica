#' Draw 3-panel figure including raw time series, time-frequency power, and average spectral density.
#'
#' @param obj : morlet wavelet object (from dplR morlet function)
#' @param interpolate
#' @param figstr
#' @param fig.title
#' @param fig.width
#' @param fig.height
#' @param fig.res
#'
#' @return
#' @import ggplot2
#' @import gridExtra
#' @export
#'
#' @examples
tfplot<-function(obj,interpolate=TRUE,
                       figstr='fig.png',fig.title='',
                       fig.width=1800, fig.height=1800,fig.res=300){
  nt = length(obj$y)
  y.c = scale(obj$y, center=TRUE, scale=FALSE)
  mat = as.matrix(Mod(obj$Power))
  colnames(mat)<-(obj$period)
  rownames(mat)<-obj$x
  perid.x = obj$period
  perid.x0=1/((1:(nt/2-1))/nt)
  y.c.fft = Mod(fft(y.c)[2:(nt/2)])^2/(2*pi*nt)
  lspec.fit = lspec(period = y.c.fft)
  pred = dlspec( 1/perid.x *2*pi,lspec.fit)
  lspec.est = pred$d + pred$m

  mat.long = melt(mat) %>%
    rename(Period=Var2, Time=Var1)
  fig1<-ggplot(data.frame(x=obj$x, y=y.c),aes(x=x,y=y))+
    geom_line() + theme_bw() +
    theme(plot.margin = margin(0.1,0, 0,0.7, "cm"),
          axis.title.x=element_blank(),
          axis.title.y=element_blank(),
          axis.text.x=element_blank(),
          title =element_text(size=20, face='bold')) +
    ggtitle(fig.title)
  fig3<-ggplot(rbind(data.frame(x=perid.x0, y=y.c.fft, g='periodogram'),
                     data.frame(x=perid.x, y=lspec.est,g='lspec')),
               aes(x=x,y=y,colour=g)) + geom_line(alpha=0.5) +
    theme_bw() + coord_flip() + scale_x_continuous(trans = 'log2',limits=range(obj$period)) +
    theme(legend.position='bottom',legend.direction = 'vertical',
          legend.key.height = unit(0.1, 'cm'),
          legend.key.width = unit(0.2, 'cm'),
          legend.margin = margin(0,0,0,0,'cm'),
          legend.box.margin=margin(0,0,0,0, "cm"),
          plot.margin = margin(0,0.2, 0.2,0.1, "cm"),
          axis.title.y=element_blank(),
          axis.text.y=element_blank(),
          legend.title =element_blank())+
    #ggtitle('Power spectrum') +
    ylab('Power')

  fig2<-ggplot(mat.long, aes(x=Time, y=Period, fill=value)) +
    geom_raster(interpolate=interpolate) +
    scale_y_continuous(trans = 'log2') +
    scale_fill_distiller(palette = "Spectral", direction = -1) +
    theme_bw()+
    #ggtitle('Time-frequecy Power (Mod)') +
    theme(legend.position='bottom',plot.margin = margin(0,0, -0.5,0, "cm"))

  lay <- rbind(c(1,1,1,NA),
               c(2,2,2,3),
               c(2,2,2,3),
               c(2,2,2,3))

  png(figstr, width=fig.width, height=fig.height, res=fig.res)
  a = grid.arrange(fig1,fig2,fig3,
                   layout_matrix = lay)
  return(list(a,fig1,fig2,fig3))
}
