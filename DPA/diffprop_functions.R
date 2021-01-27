
# Differential proportion analysis of single-cell populations

###
# Simulation of a data set with different number of cells and their cell types
###
simCells<-function(n=10^7, prop=c(0.6, 0.3, 0.05, 0.03, 0.01, 0.005, 0.005)){
	ct<-paste("Cell_Type_",1:length(prop),sep="");
	nct<-round(n*prop);
	lab<-rep(ct,nct);	
	return(lab);
}

###
# Simulate a data-set with some provided error rate to simulate experimental/biological noise
###
simCellsNoisy<-function(n=10^7, prop, error = 0.05){
  # Create a 'noisy' distribution for sampling cells
  prop.noisy = unlist(lapply(prop, function(x) 
                     {ifelse(sample(c(1:2), size=1) == 1, x+x*error, x-x*error)}))
  prop.noisy = prop.noisy/sum(prop.noisy)
  
  # Sample as before but with a modified distribution
  ct<-paste("Cell_Type_",1:length(prop.noisy),sep="");
  nct<-round(n*prop.noisy);
  lab<-rep(ct,nct);	
  return(lab);
}

###
# Subsampling of cells
###
subCells<-function(cells, n=round(length(cells)/10)){
	sample(cells,n);
}



###
# Create the condition and cluster variables from the observed contingency table
###
makeVariables<-function(obs.table){
   cond_names<-dimnames(obs.table)[[1]];
   clus_names<-dimnames(obs.table)[[2]];
   
   cond<-rep(cond_names, apply(obs.table,1,sum) );
   clus_list<-apply(obs.table,1,function(x){
     rep(clus_names, x );
   })
   clus<-base::do.call(c,clus_list);
   #table(cond,clus);
   return(list(cond=cond, clus=clus));
}

###
# Generate null distributions based on a permutation procedure
###
generateNull<-function(obs.table, n=10000, p=0.2){
    obs<-makeVariables(obs.table);
    all.exp<-array(NA,dim=c(dim(obs.table)[1],dim(obs.table)[2],n))
    dimnames(all.exp)[[1]]<-dimnames(obs.table)[[1]];
    dimnames(all.exp)[[2]]<-dimnames(obs.table)[[2]];
    dimnames(all.exp)[[3]]<-1:n;
    clus_names<-dimnames(obs.table)[[2]];
    cond_names<-dimnames(obs.table)[[1]];
    
    v<-makeVariables(obs.table);
    
    ### Perform permutation n times
    for(i in 1:n){
      pv<-v;
      
      if (i %% 1000 == 0){
        print(i);}
      
      ### Permutation: Randomly select p*100% of points to be re-sampled from the background 
      randInd<-sample(1:length(pv$clus),round(length(pv$clus)*p));
      pv$clus[randInd]<-sample(v$clus,length(randInd),replace=F);  
      
      this.exp<-table(pv$cond,pv$clus);
      exp.table<-obs.table;
      exp.table[1:dim(exp.table)[1],1:dim(exp.table)[2]]<-NA
      exp.table[dimnames(this.exp)[[1]],dimnames(this.exp)[[2]]]<-this.exp;
      exp.table
      all.exp[,,i]<-exp.table;
    }
    return(all.exp);
}

###
# Perform sum by ignoring NA
###
sumNoNA<-function(x){
  sum(x[which(!is.na(x))])
}

###
# Perform a two-class test of significance
###
two.class.test<-function(obs.table, all.exp, cond.control="C", cond.treatment="PA",to.plot=T){
  clus_names<-dimnames(obs.table)[[2]]
  pp<-array(-1,length(clus_names));
  names(pp)<-clus_names;

  if(to.plot){
     par(mfrow=c(3,ceiling(length(clus_names)/3)));
  }
  for(this_clus in clus_names){
    obs.diff<-obs.table[cond.treatment,this_clus]/sum(obs.table[cond.treatment,]) - 
              obs.table[cond.control,this_clus]/sum(obs.table[cond.control,]);
    all.diff<-all.exp[cond.treatment,this_clus,]/apply(all.exp[cond.treatment,,],2,sumNoNA) -
              all.exp[cond.control,this_clus,]/apply(all.exp[cond.control,,],2,sumNoNA);
    if(to.plot){
       hist(all.diff,breaks=50,col="grey",border="grey",main=this_clus)
       abline(v=obs.diff,col="red")
    }
    p.r<-length(which(obs.diff>all.diff))/length(which(!is.na(all.diff)));
    p.l<-length(which(obs.diff<all.diff))/length(which(!is.na(all.diff)));
    pp[this_clus]<-min(p.r,p.l);
  }
  return(pp);
}


