library(getopt)
options(bitmapType = "cairo")
###传参信息
spec <- matrix(c(
  'file', 'f', 0, 'character',
  'name', 'n', 0, 'character',
  'insert_name', 'i', 0, 'character'), byrow = TRUE, ncol = 4)
opt <- getopt(spec)
print_usage <- function(spec=NULL){
  q(status = 1);
}
if ( !is.null(opt$help))   { print_usage(spec) }
if ( is.null(opt$file))   { print_usage(spec)}
if ( is.null(opt$name))   { print_usage(spec)}
if ( is.null(opt$insert_name))   { print_usage(spec)}
library(ggplot2)

a<-read.table(opt$file,sep="\t",header=FALSE)
colnames(a) = c("Chr","Position","Depth")
b<-readLines(opt$insert_name)
#b <- ggplot(a, aes(x = a$Position, y = a$Depth)) +
#  geom_bar(stat = "identity", position = "dodge") +
#  labs(title = "My Bar Plot", x = "Position", y = "Depth")
#ggsave(paste0(opt$name,".png"), plot = b)
for (i in b){
  p<- ggplot(a[a$Chr== i, ])+geom_col(aes(x=Position,y=Depth),col=NA,fill="royalblue4")+geom_hline(yintercept=mean(a$Depth),linetype="dashed")+labs(title=paste0("Tdna sequence ",i," coverage in ",opt$name))+theme_bw()+theme(plot.title.position = "plot")
  ggsave(paste0(opt$name,"_",i,".png"), plot = p,width = 16, height = 9)
}
