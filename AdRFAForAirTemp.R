# Function Add_AFR
Add_AFR=function(temperaturas){
  if (! require (simstudy)) {install.packages ("simstudy")}
  if (! require (lmomRFA)) {install.packages ("lmomRFA")}
  if (! require (rrcov)) {install.packages ("rrcov")}
  library(lmomRFA)
  library(simstudy)
  library(rrcov)
question=menu(c("1 (mais rigoroso), digite 1 ", "1.28, digite 2","1.64 , digite 3", "2 (menos rigoroso), digite 4"),
            title="Escolha o valor crítico para a medidade H de heterogeniedade?")
if(question==1){H.crit=1}
if(question==2){H.crit=1.28}
if(question==3){H.crit=1.64}
if(question==4){H.crit=2}
question=menu(c("Totalmente independentes (sigma=0), digite 1 ", "sigma=0.1, digite 2","sigma=0.2, digite 3",
             "sigma=0.3, digite 4","sigma=0.4, digite 5","sigma=0.5, digite 6","sigma=0.6, digite 7","sigma=0.7, digite 8",
             "sigma=0.8, digite 9","sigma=0.9 (elevada dependencia), digite 10")
           , title="Escolha o valor medio (sigma) para dependencia espacial entre as localidades?")
if(question==1){sigma=0}
if(question==2){sigma=0.1}
if(question==3){sigma=0.2}
if(question==4){sigma=0.3}
if(question==5){sigma=0.4}
if(question==6){sigma=0.5}
if(question==7){sigma=0.6}
if(question==8){sigma=0.7}
if(question==9){sigma=0.8}
if(question==10){sigma=0.9}
data.temp=temperaturas 
Ns=1000
message("Calculando. Por favor aguarde alguns instantes")
n.years=length(temperaturas[,1])
n.locals=length(temperaturas[1,])

  x1.atoutset=regsamlmu(temperaturas, lcv = FALSE)
  med=mean(x1.atoutset$l_1)
  mcd=CovMcd(x1.atoutset[,4:6])
  New.D=as.matrix(sqrt(getDistance(mcd)))
  lmom.atsite=matrix(NA,n.locals,4)
  lmom.atsite.sim=matrix(NA,n.locals,4)
  vetor.numerador=matrix(NA,n.locals,1)
  at.site.par=matrix(NA,n.locals,3)
  V.sim=matrix(NA,Ns,1)
  numerador.bias=matrix(NA,Ns,1)
  numerador.bias.quad=matrix(NA,Ns,1)
  sknumerador.bias=matrix(NA,Ns,1)
  sknumerador.bias.quad=matrix(NA,Ns,1)
  kurt.sKew.numerador.bias=matrix(NA,Ns,1)
  for (site in 1:n.locals){
    data.temp[,site]=temperaturas[,site]-mean(temperaturas[,site],na.rm=TRUE)
    lmom.atsite[site,]=samlmu(data.temp[,site])}
  rmom=colMeans(lmom.atsite)
  para=try(pelkap(rmom),TRUE)
  for (v in 1:n.locals){
    vetor.numerador[v]=x1.atoutset[v,2]*(((lmom.atsite[v,2])-rmom[2])^2)}
  V=sqrt(sum(vetor.numerador)/sum(x1.atoutset[,2]))
  for (ns in 1:Ns){
    z=genCorGen(n.years, nvars = n.locals, params1 = 0, params2 = 1, dist = "normal", 
                corstr = "ar1", rho=sigma, wide = TRUE)
    z1=as.matrix(z[,2:(n.locals+1)])
    u <- pnorm(z1)
    u.sim <- pnorm(z1)
    
    bloco.kap.sim=(quakap(u.sim,c(para[1],para[2],para[3],para[4])))
    ff=1; for (ff in 1:n.locals){
      bloco.kap.sim[,ff]=bloco.kap.sim[,ff]
      lmom.atsite.sim[ff,]=samlmu(bloco.kap.sim[,ff])}
    rmom.sim=colMeans(lmom.atsite.sim)
    
    for (v in 1:n.locals){
      vetor.numerador[v]=x1.atoutset[v,2]*(((lmom.atsite.sim[v,2])-rmom.sim[2])^2)}
    
    V.sim[ns,1]=sqrt((sum(vetor.numerador)/sum(x1.atoutset[,2])))
    
    numerador.bias[ns]=rmom.sim[4]-rmom[4]
    numerador.bias.quad[ns]=(rmom.sim[4]-rmom[4])^2
    sknumerador.bias[ns]=rmom[3]-rmom.sim[3]
    sknumerador.bias.quad[ns]=(rmom[3]-rmom.sim[3])^2
    kurt.sKew.numerador.bias[ns]=(rmom[3]-rmom.sim[3])*(rmom[4]-rmom.sim[4])
  }
  H=as.matrix((V-mean(V.sim))/sd(V.sim))
  colnames(H)=c("Medida de heterogeneidade")
  rownames(H)=c("H")
  print(H)
  if (H<=H.crit){
    
  message("A regiao pode ser considerada homogenea")
  x1=regsamlmu(temperaturas)
  regional.lmom=regavlmom(x1)
  para.gev=pelgev(regional.lmom)
  para.glo=pelglo(regional.lmom)
  para.gum=pelgum(regional.lmom)
  para.pe3=pelpe3(regional.lmom)
  para.gno=pelgno(regional.lmom)
  para.gp=pelgpa(regional.lmom)
  B4=sum(numerador.bias)/Ns
  SD=sqrt((sum(numerador.bias.quad)-(Ns*(B4^2)))/(Ns-1))
  lmomgev=lmrgev(para = para.gev, nmom = 4)
  lmomglo=lmrglo(para = para.glo, nmom = 4)
  lmomgum=lmrgum(para = para.gum, nmom = 4)
  lmompe3=lmrpe3(para = para.pe3, nmom = 4)
  lmomgno=lmrgno(para = para.gno, nmom = 4)
  lmomgp=lmrgpa(para = para.gp, nmom = 4)
  Z.GEV=(lmomgev[4]-rmom[4]+B4)/SD
  Z.GLO=(lmomglo[4]-rmom[4]+B4)/SD
  Z.gum=(lmomgum[4]-rmom[4]+B4)/SD
  Z.PE3=(lmompe3[4]-rmom[4]+B4)/SD
  Z.GNO=(lmomgno[4]-rmom[4]+B4)/SD
  Z.gp=(lmomgp[4]-rmom[4]+B4)/SD
  Z=cbind(Z.GEV,Z.GLO,Z.gum,Z.PE3,Z.GNO,Z.gp)
  ###########
  B3=sum(sknumerador.bias)/Ns
  SD.Skew=sqrt(((sum(sknumerador.bias.quad))-Ns*(B3^2))/(Ns-1))
  SD.kurt.Skew=(sum(kurt.sKew.numerador.bias)-(Ns*B3*B4))/(Ns-1)
  StR<-matrix(c(SD.Skew^2,SD.kurt.Skew,SD.kurt.Skew,SD^2),nrow=2,ncol=2)
  tb=c(rmom[3]-B3,rmom[4]-B4)
  Tdist=c(lmomgev[3],lmomgev[4])
  Zbivar.GEV=mahalanobis(Tdist, tb, StR)
  Tdist=c(lmomglo[3],lmomglo[4])
  Zbivar.GLO=mahalanobis(Tdist, tb, StR)
  Tdist=c(lmomgum[3],lmomgum[4])
  Zbivar.gum=mahalanobis(Tdist, tb, StR)
  Tdist=c(lmompe3[3],lmompe3[4])
  Zbivar.PE3=mahalanobis(Tdist, tb, StR)
  Tdist=c(lmomgno[3],lmomgno[4])
  Zbivar.gno=mahalanobis(Tdist, tb, StR)
  Tdist=c(lmomgp[3],lmomgp[4])
  Zbivar.gp=mahalanobis(Tdist, tb, StR)
  Zbivar=cbind(Zbivar.GEV,Zbivar.GLO,Zbivar.gum,Zbivar.PE3,Zbivar.gno,Zbivar.gp)
  Goodness=rbind(Z,Zbivar)
  colnames(Goodness)=c("GEV","GLO","GUM","PE3","GNO","GPA")
  if (n.locals<=15){
    best=names(which.min(Goodness[1,]))
  }else{best=names(which.min(Goodness[2,]))}
  
  means=as.matrix(x1$l_1)
  lmom.atsite=matrix(NA,n.locals,4)
  lmom.atsite.weighted=matrix(NA,n.locals,4)
  for (j in 1:n.locals){
  lmom.atsite[j,]=samlmu(data.temp[,j])
  lmom.atsite.weighted[j,]=lmom.atsite[j,]*x1$n[j]}
  lmom.atsite.weighted.final=colSums(lmom.atsite.weighted)/sum(x1$n)
  temp=as.numeric(readline(prompt="Forneça o valor da temperatura do ar para calculo da frequencia esperada(%) "))
  means=temp-means
  if (best=="GEV"){
    message("Distribuição selecionada:Generalizada dos valores extremos")
    para.reg=pelgev(lmom.atsite.weighted.final)
    prob=as.matrix(100*cdfgev(means,para.reg))
  }
  if (best=="GLO"){
    message("Distribuição selecionada:Generalizada logistica")
    para.reg=pelglo(lmom.atsite.weighted.final)
    prob=as.matrix(100*cdfglo(means,para.reg))
  }
  if (best=="GUM"){
    message("Distribuição selecionada:Gumbel")
    para.reg=pelgum(lmom.atsite.weighted.final)
    prob=as.matrix(100*cdfgum(means,para.reg))
  }
  if (best=="PE3"){
    message("Distribuição selecionada:Pearson tipo III")
    para.reg=pelpe3(lmom.atsite.weighted.final)
    prob=as.matrix(100*cdfpe3(means,para.reg))
  }
  if (best=="GNO"){
    message("Distribuição selecionada:Nornal generalizada (3 parametros)")
    para.reg=pelgno(lmom.atsite.weighted.final)
    prob=as.matrix(100*cdfgno(means,para.reg))
  }
  if (best=="GPA"){
    message("Distribuição selecionada:Generalizada Pareto")
    para.reg=pelgpa(lmom.atsite.weighted.final)
    prob=as.matrix(100*cdfgpa(means,para.reg))
  }
  nomes=as.matrix(colnames(temperaturas))
  message("Probabilidade de ocorrencia de valores iguais ou inferiores a temperatura fornecida")
  rownames(prob)=nomes
  colnames(prob)=("Freq(%)")
  print(prob)
  } else{
    message("A regiao nao pode ser considerada homogenea")
    message("A AFR nao pode ser aplicada")
    message("recomenda-se remorver as localidades com os valores mais elevados da medida de discordancia e re-calcular H")
    nomes=colnames(temperaturas)
    rownames(New.D)=nomes
    colnames(New.D)=c("discordancia")
    print(New.D)
  }}

# O arquivos meusdados.csv deve conter dados de temperaturas do ar máxima ou minimas de pelo menos sete pontos de coleta.
# Cada Columa deve conter os dados de cada ponto.
# Aviso para softwares como o R e R-studio, o separador decinal é o ".".
# Veja o arquivo meusdados.csv como modelo.

temperaturas=read.csv("meusdados.csv",sep=",", dec=".", header = TRUE)
Add_AFR(temperaturas)
  