betaparametab <-
    function(xn,w,P,V){        
n<-length(xn)
W2<-0
for(i in seq(n)){
W2<-W2+w[i]^2	
}
W2X<-0
for(i in seq(n)){
W2X<-W2X+w[i]^2/xn[i]	
}  	
b<-(P*(1-P)*W2-V)/(
(V/(1-P))-P*W2X)
a<-1
a<-b*P/(1-P)
c<-c(a,b)	
return(c)
}
