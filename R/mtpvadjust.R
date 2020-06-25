mtpvadjust <-
function(pv,C){
	# pv from largest to smallest
	N<-length(pv)
 
	q<-rep(1,N)
	for(i in seq(N)){
		q[i]<-i/N
		}
	S<-1	
	SS<-rep(1,N)
	for(k in 2:N){
		if(q[k]==0&q[k-1]==0){
		S<-S+1	
			}else{
				S<-S+(q[k]-q[k-1])/(q[k]+q[k-1])
				
				}
				for(i in 2:k){
					if(q[i]==0&q[i-1]==0){
						SS[k]<-SS[k-1]+1
						}else{
							SS[k]<-SS[k-1]+(q[i]-q[i-1])/(q[i]+q[i-1])
							}
					}
		       
		}
	 R<-rep(1,N)	
	 paj<-rep(1,N) # adjusted pvalue
		for(k in seq(N)){
			R[k]<-(SS[k]/S)^C
			
			}
	R<-sort(R, decreasing=TRUE)	
	pv<-sort(pv,decreasing=TRUE)
	paj[1]<-R[1]
	for (i in 2:N){
	paj[i]<-min(paj[i-1],(pv[i]*(N/(N^R[N-i+1]))))
		}
		
	return(paj)	
		
	}
