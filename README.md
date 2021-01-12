# Network_Cores
C++ Algorithms for network core extraction/k-core/M-core/GK-core/Weak core and related observables

Author: Bernat Corominas-Murtra (2021)    

Based on the following Bibliography

*K-core:

	  
	K-core organization of complex networks
   	S.N. Dorogovtsev, A.V. Goltsev, J.F.F. Mendes
	Phys. Rev. Lett. 96, 040601 (2006)

*Generalized K-core:

	Detection of the Elite Structure in a Virtual Multiplex Social System by Means of a Generalised K-Core
   	B. Corominas-Murtra, B. Fuchs, S. Thurner PloS one 9 (12), e112606 (2014)

*M-core:

	Deciphering the global organization of clustering in real complex networks
	Pol Colomer-de-Simón, M. Ángeles Serrano, Mariano G. Beiró, J. Ignacio Alvarez-Hamelin & Marián Boguñá
	Scientific Reports 3, Article number: 2517 (2013)

*Weak Core & Core_Cluster

	The weak core and the structure of elites in social multiplex networks
	B. Corominas-Murtra and S. Thurner
	Springer Complexity: Interconnected Networks, Chapter 10, (2015)
	ArXiv preprint: https://arxiv.org/abs/1412.6913

****Contents

Header Cores_subgraph.h includes:										 

    Graph connected Component structure		
    
    K-core										 
    
    Generalized K-core					
    
    M-core										 
    
    Weak Core									 
		
    Core Clusters (MCore of GKcore)				 
		
    KCoreness									 
 		
    GKCoreeness									 
		
    Mcoreness									 
													 
main.cpp includes an example of computation of the above core subgraph-related observables. Some example networks are provided in order to test the algorithms.

***WARNING 

	    GKCoreeness	may be time expensive in dense and massive graphs with strong hubs								 



