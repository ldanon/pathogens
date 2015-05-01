setwd('/Users/ldanon/tmp/DengueLatest/Model/')

system("./a.out 2")


inf<-read.table('Infecteds.dat',)
not_inf<-read.table('not_inf.dat')
not_imm<-read.table('not_imm.dat')

matplot(inf,type="l")
