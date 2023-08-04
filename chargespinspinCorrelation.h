#include "itensor/all.h"
#include "itensor/util/print_macro.h"

using namespace itensor;
//~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
////computing charge spin spin correlation
////~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
void chargespinspinCorrelator(MPS psi, const SiteSet sites, int c){
	auto N = length(psi);
	printfln("i j c <Nhc Si^z Sj^z> <Nhc Si^+ Sj^-> <Nhc Si^- Sj^+>");
        for (int i=1;i<=N;i++)
	{ 
         for (int j=i;j<=N;j++)
	    {
	     if (i<c)
	     {   psi.position(i);
		auto Cz = psi(i)*sites.op("Sz",i);
		auto Cpm = psi(i)*sites.op("S+",i);
	        auto Cmp = psi(i)*sites.op("S-",i);
		if (j==i)
		 {
	          Cz *= prime(sites.op("Sz",i));
	          Cz *= dag(prime(prime(psi(i),"Site"),"Site"));
	          Cpm *= prime(sites.op("S-",i));
		  Cpm *= dag(prime(prime(psi(i),"Site"),"Site"));
		  Cmp *= prime(sites.op("S+",i));
		  Cmp *= dag(prime(prime(psi(i),"Site"),"Site"));
		  auto ir = commonIndex(psi(i),psi(i+1),"Link");
		  Cz *= dag(prime(prime(psi(i),"Site"),ir));
		  Cpm *= dag(prime(prime(psi(i),"Site"),ir));
		  Cmp *= dag(prime(prime(psi(i),"Site"),ir));
		  for (int k=i+1;k<c;k++){
			  Cz *= psi(k);
			  Cz *= dag(prime(psi(k),"Link"));
			  Cpm *= psi(k);
			  Cpm *= dag(prime(psi(k),"Link"));
			  Cmp *= psi(k);
			  Cmp *= dag(prime(psi(k),"Link"));}
		  Cz *= psi(c)*sites.op("Nh",c);
		  Cpm *= psi(c)*sites.op("Nh",c);
		  Cmp *= psi(c)*sites.op("Nh",c);
		  auto il = commonIndex(psi(c),psi(c-1),"Link");
		  Cz *= dag(prime(prime(psi(c),"Site"),il));
		  Cpm *= dag(prime(prime(psi(c),"Site"),il));
		  Cmp *= dag(prime(prime(psi(c),"Site"),il));
                  } // closes if j==i 
		 else if ((j>i) && (j>c))
		 { auto ir = commonIndex(psi(i),psi(i+1),"Link");
                   Cz *= dag(prime(prime(psi(i),"Site"),ir));
		   Cpm *= dag(prime(prime(psi(i),"Site"),ir));
		   Cmp *= dag(prime(prime(psi(i),"Site"),ir));
                   for (int k=i+1;k<c;k++){
			    Cz *= psi(k);
	 		    Cz *= dag(prime(psi(k),"Link"));
			    Cpm *= psi(k);
			    Cpm *= dag(prime(psi(k),"Link"));
			    Cmp *= psi(k);
			    Cmp *= dag(prime(psi(k),"Link"));}		    
		Cz *= psi(c)*sites.op("Nh",c);
	        Cpm *= psi(c)*sites.op("Nh",c);
		Cmp *= psi(c)*sites.op("Nh",c);
		auto il = commonIndex(psi(c),psi(c-1),"Link");
	        Cz *= dag(prime(prime(psi(c),"Site"),il));
	        Cpm *= dag(prime(prime(psi(c),"Site"),il));
	        Cmp *= dag(prime(prime(psi(c),"Site"),il));
		auto irr = commonIndex(psi(c),psi(c+1),"Link");
		Cz *= dag(prime(prime(psi(c),"Site"),irr));
		Cpm *= dag(prime(prime(psi(c),"Site"),irr));
		Cmp *= dag(prime(prime(psi(c),"Site"),irr));
         	for (int h=i+1;h<j;h++){
		     Cz *= psi(h);
		     Cz *= dag(prime(psi(h),"Link"));
		     Cpm *= psi(h);
		     Cpm *= dag(prime(psi(h),"Link"));
		     Cmp *= psi(h);
		     Cmp *= dag(prime(psi(h),"Link"));}
		Cz *= psi(j)*sites.op("Sz",j);
		Cpm *= psi(j)*sites.op("S-",j);
	        Cmp *= psi(j)*sites.op("S+",j);
	        auto ill = commonIndex(psi(j),psi(j-1),"Link");
		Cz *= dag(prime(prime(psi(j),"Site"),ill));
		Cpm *= dag(prime(prime(psi(j),"Site"),ill));
	        Cmp *= dag(prime(prime(psi(j),"Site"),ill));
		} //closes if j>i and j>c
		 else if ((j>i) && (j<c))
		 { auto ir = commonIndex(psi(i),psi(i+1),"Link"); 
			Cz *= dag(prime(prime(psi(i),"Site"),ir)); 
			Cpm *= dag(prime(prime(psi(i),"Site"),ir));  
			Cmp *= dag(prime(prime(psi(i),"Site"),ir));
			for (int k=i+1;k<j;k++){
				Cz *= psi(k);
				Cz *= dag(prime(psi(k),"Link"));
				Cpm *= psi(k);
				Cpm *= dag(prime(psi(k),"Link"));
				Cmp *= psi(k);
				Cmp *= dag(prime(psi(k),"Link"));}
			Cz *= psi(j)*sites.op("Sz",j);
			Cpm *= psi(j)*sites.op("S-",j);
			Cmp *= psi(j)*sites.op("S+",j);
			auto il = commonIndex(psi(j),psi(j-1),"Link");
			Cz *= dag(prime(prime(psi(j),"Site"),il));
			Cpm *= dag(prime(prime(psi(j),"Site"),il));  
			Cmp *= dag(prime(prime(psi(j),"Site"),il));
			auto irr = commonIndex(psi(j),psi(j+1),"Link"); 
			Cz *= dag(prime(prime(psi(j),"Site"),irr));  
			Cpm *= dag(prime(prime(psi(j),"Site"),irr)); 
			Cmp *= dag(prime(prime(psi(j),"Site"),irr));
			for (int h=i+1;h<c;h++){ 
				Cz *= psi(h);
				Cz *= dag(prime(psi(h),"Link"));  
				Cpm *= psi(h);
				Cpm *= dag(prime(psi(h),"Link"));
				Cmp *= psi(h); 
				Cmp *= dag(prime(psi(h),"Link"));}
			Cz *= psi(c)*sites.op("Nh",c);
			Cpm *= psi(c)*sites.op("Nh",c);  
			Cmp *= psi(c)*sites.op("Nh",c);
			auto ill = commonIndex(psi(c),psi(c-1),"Link");  
			Cz *= dag(prime(prime(psi(c),"Site"),ill));  
			Cpm *= dag(prime(prime(psi(c),"Site"),ill));
			Cmp *= dag(prime(prime(psi(c),"Site"),ill));
                 } // closes if j>i and j>c
	printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp));	
	     } // closes if i<c
	          else if (i>c)
		  { psi.position(c);           
			 auto Cz = psi(c)*sites.op("Nh",c);
			 auto Cpm = psi(c)*sites.op("Nh",c);
			 auto Cmp = psi(c)*sites.op("Nh",c);
			  auto ir = commonIndex(psi(c),psi(c+1),"Link");
			  Cz *= dag(prime(prime(psi(c),"Site"),ir));
			  Cpm *= dag(prime(prime(psi(c),"Site"),ir));
			  Cmp *= dag(prime(prime(psi(c),"Site"),ir)); 
			  if (j==i)
			  { for (int k=c+1;k<i;k++){
				  Cz *= psi(k);
				  Cz *= dag(prime(psi(k),"Link"));
				  Cpm *= psi(k);
				  Cpm *= dag(prime(psi(k),"Link"));
				  Cmp *= psi(k);
				  Cmp *= dag(prime(psi(k),"Link"));}
			  Cz *= psi(i)*sites.op("Sz",i);
			  Cpm *= psi(i)*sites.op("S+",i);
			  Cmp *= psi(i)*sites.op("S-",i);
                          Cz *=  prime(sites.op("Sz",i));
			  Cz *= dag(prime(prime(psi(i),"Site"),"Site"));
			  Cpm *= prime(sites.op("S-",i));
			  Cpm *= dag(prime(prime(psi(i),"Site"),"Site"));
			  Cmp *= prime(sites.op("S+",i));
			  Cmp *= dag(prime(prime(psi(i),"Site"),"Site"));
			  auto il = commonIndex(psi(j),psi(j-1),"Link");
			  Cz *= dag(prime(prime(psi(j),"Site"),il));
			  Cpm *= dag(prime(prime(psi(j),"Site"),il));
			  Cmp *= dag(prime(prime(psi(j),"Site"),il));
			  } // closes if j==i  
			  else if (j>i)
			  {
				      for (int k=c+1;k<i;k++){
					  Cz *= psi(k);
					  Cz *= dag(prime(psi(k),"Link"));
					  Cpm *= psi(k);
					  Cpm *= dag(prime(psi(k),"Link"));
					  Cmp *= psi(k);
					  Cmp *= dag(prime(psi(k),"Link"));}
				     Cz *= psi(i)*sites.op("Sz",i);
				     Cpm *= psi(i)*sites.op("S+",i);
				     Cmp *= psi(i)*sites.op("S-",i);
				     auto il = commonIndex(psi(i),psi(i-1),"Link");
				     Cz *= dag(prime(prime(psi(i),"Site"),il));
				     Cpm *= dag(prime(prime(psi(i),"Site"),il));
				     Cmp *= dag(prime(prime(psi(i),"Site"),il));
				     auto irr = commonIndex(psi(i),psi(i+1),"Link");
				     Cz *= dag(prime(prime(psi(i),"Site"),irr));
				     Cpm *= dag(prime(prime(psi(i),"Site"),irr));
				     Cmp *= dag(prime(prime(psi(i),"Site"),irr));
				     for (int h=i+1;h<j;h++){
					     Cz *= psi(h);
					     Cz *= dag(prime(psi(h),"Link"));
					     Cpm *= psi(h);
					     Cpm *= dag(prime(psi(h),"Link"));
					     Cmp *= psi(h); 
					     Cmp *= dag(prime(psi(h),"Link"));}
				     Cz *= psi(j)*sites.op("Sz",j);
				     Cpm *= psi(j)*sites.op("S-",j);
				     Cmp *= psi(j)*sites.op("S+",j); 
				     auto ill = commonIndex(psi(j),psi(j-1),"Link"); 
				     Cz *= dag(prime(prime(psi(j),"Site"),ill)); 
				     Cpm *= dag(prime(prime(psi(j),"Site"),ill));
				     Cmp *= dag(prime(prime(psi(j),"Site"),ill));
			       } // close j>i in i>c

		printfln("",i," ",j," ",elt(Cz)," ",elt(Cpm)," ",elt(Cmp));
		  }//closes if i>c
	         } // closes for (int j=i;j<=N;j++)
	         } // closes for (int i=1;i<=N;i++)
	         } // closes the function

