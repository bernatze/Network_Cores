#ifndef CORES_SUBGRAPH
#define CORES_SUBGRAPH

/******************************************************
**   CORE EXTRACTING Functions                       **
**   Bernat Corominas-Murtra                         **
**   Last revision: 16.11.2017                       **
**   Vienna Complexity Science Hub					 **
**   Medical University of Vienna					 **
**   Section for the Science of Complex Systems		 **
**	 Includes:										 **
**		Graph Component structure					 **
**		K-core										 **	
**		Generalized K-core							 **
**		M-core										 **
**		Weak Core									 **
**		Core Clusters (MCore of GKcore)				 **
**		KCoreness									 **
**		GKCoreeness									 **
**		Mcoreness									 **
**													 **
******************************************************/

/******************************************************
* Bibliography
*
*	K-core:
*   K-core organization of complex networks
*   S.N. Dorogovtsev, A.V. Goltsev, J.F.F. Mendes
*	Phys. Rev. Lett. 96, 040601 (2006)
*
*   Generalized K-core:
*   Detection of the Elite Structure in a Virtual Multiplex Social System by Means of a Generalised K-Core
*   B. Corominas-Murtra, B. Fuchs, S. Thurner PloS one 9 (12), e112606 (2014)
*
*   M-core:
*	Deciphering the global organization of clustering in real complex networks
*	Pol Colomer-de-Simón, M. Ángeles Serrano, Mariano G. Beiró, J. Ignacio Alvarez-Hamelin & Marián Boguñá
*	Scientific Reports 3, Article number: 2517 (2013)
*
*	Weak Core & Core_Cluster
*	The weak core and the structure of elites in social multiplex networks
*	B. Corominas-Murtra and S. Thurner
*	Springer Complexity: Interconnected Networks, Chapter 10, (2015)
*
******************************************************/

#include <iostream>
#include <vector>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>
#include <fstream>
#include <cstdio>
#include <iterator> 
#include "Cores_subgraph.h"

using namespace std;

namespace Cores_BCM
{

class Cores {

    public:

        vector<vector<int> > graph; //Input graph
	
	//Core extracting Functions
        vector<int> Nodes_Graph();
        vector<vector<int> >  Copy_graph () {return graph;};
        void Component_Structure (const vector<vector<int> >, vector<vector<int> > &);
        void Kcore(vector<vector<int> >, int, vector<vector<int> > &);
        void GKcore(vector<vector<int> > , int, vector<vector<int> > &);
        void M_core (vector<vector<int> >, vector<vector<int> >, vector<vector<int> > &, int);
        void Core_clusters (vector<vector<int> >, int, int, vector<vector<int> > &C);
        void Weak_core(const vector<vector<int> >, int, int, vector<vector<int> > &);
        void KCoreness(const vector<vector<int> >, vector<vector<int> > &);
        void GKCoreness(const vector<vector<int> >, vector<vector<int> > &);
        void MCoreness(const vector<vector<int> >, vector<vector<int> > &);
};

/***************** Elementary functions *************************/
/*void printVector(vector<int> v, ostream &output_file)*/
/*void printVectorString(vector<string> v, ostream &output_file)*/
/*void printVectorScreen(vector<int> v)*/
/*void getGraph(istream &input_file, vector<vector<int> > &v)*/
/*void printGraph(ostream &output_file, vector<vector<int> > a)*/
/*void printGraph_Screen(vector<vector<int> > a)*/

vector<int> Cores::Nodes_Graph()
{
    vector<int>  A;
    A.clear();
    if(graph.size()>0)
    {
        vector< vector<int> >::const_iterator it;
        vector<int>::const_iterator edg;
        int k=0;
        for(it=graph.begin(); it != graph.end(); ++it)
        {
            for (edg = it->begin(); edg != it->end(); edg++)
            {
                int h=0;
                int i;
                for(i=0; i<k; ++i)
                {
                    if(A[i]==*edg)
                    {
                        h=1;
                        break;
                    }
                }
                if(h==0)
                {
                    A.push_back(*edg);
                    k=k+1;
                }
            }
        }
    }
    return A;
}

void printVector(vector<int> v, ostream &output_file)
{
	for(int i=0; i< v.size(); i++)
		output_file << v[i] <<endl;
}

void printVectorString(vector<string> v, ostream &output_file)
{
	for(int i=0; i< v.size(); i++)
		output_file << v[i] <<endl;
}

void printVectorScreen(vector<int> v)
{
    for(int i=0; i< v.size(); i++)
        cout<< v[i] <<endl;
}

/*****Function getGraph(dataFile,  2-vector graph  v)
creates a 2-vector representing the edge  list of the 
graph ***WARNING*** To work with, function  CheckDupl
must be applied over v to avoid duplicated edges ***/

void getGraph(istream &input_file, vector<vector<int> > &v)
{
	int TEMP;
	int k=0;
	int i;
	while(!input_file.eof())
	{
        if(!input_file)
			break;
		else
		{	
			v.push_back(vector <int>());
			for(i=0; i<2; ++i)
			{
				input_file >> TEMP;
				v[k].push_back(TEMP);
				//TEMP=-1;	
			}
			k=k+1;
		}
	}	
}

/**Prints the graph as a two column adjacency list into the file
specified in &output_file *************************************/

void printGraph(ostream &output_file, vector<vector<int> > a)
{
	if(a.size()>0)
	{
		vector< vector<int> >::iterator it;
		vector<int>::iterator edg;
		for(it=a.begin(); it != a.end(); ++it)	
		{
			for (edg = it->begin(); edg != it->end(); edg++)
			{
				output_file << *edg << " ";
			}
			output_file << endl;
		}
	}
}

void printGraph_Screen(vector<vector<int> > a)
{
	if(a.size()>0)
	{
		vector< vector<int> >::iterator it;
		vector<int>::iterator edg;
		for(it=a.begin(); it != a.end(); ++it)	
		{
			for (edg = it->begin(); edg != it->end(); edg++)
			{
				cout << *edg << " ";
			}
			cout << endl;
		}
	}
}

/***************** Elementary graph functions *************************/
/*void Nodes(const vector<vector<int> > a, vector<int>  &A)*/
/*void Connectivity(const vector<vector<int> > a, const vector<int>  A, vector<int> &B)*/
/*void First_Neighbors_List (vector<vector<int> > g, vector<vector<int> > &b)*/
/*void Connected_Component (const vector<vector<int> > g, vector<vector<int> > &b)*/
/*void Component_Structure (const vector<vector<int> > g, vector<vector<int> > &b)*/
/*void Pairs_merged_components(const vector<vector<int> > g,  const vector<vector<int> > Comp, vector<vector<int> > &gCompMerged)*/

/****Function Nodes(Original graph, vector of nodes) 
Creates a vector of all the nodes of the graph.***/

void Nodes(const vector<vector<int> > a, vector<int>  &A)
{
	A.clear();
	if(a.size()>0)
	{
		vector< vector<int> >::const_iterator it;
		vector<int>::const_iterator edg;
		int k=0;
		for(it=a.begin(); it != a.end(); ++it)	
		{
			for (edg = it->begin(); edg != it->end(); edg++)
			{
				int h=0;
				int i;
				for(i=0; i<k; ++i)
				{
					if(A[i]==*edg)
					{
						h=1;
						break;
					}
				}
				if(h==0)
				{
					A.push_back(*edg);
					k=k+1;
				}	
			}				
		}	
	}
}

/******** Vector A: vector whose position i contains the integer id of node "i"
 Computed from "Nodes" Function from the original edge list vector<vector<int>>
 ATTENTION: The forthcoming computations must be performed AFTER the removal of
 duplicated edges*************************************************************/

void Connectivity(const vector<vector<int> > a, const vector<int>  A, vector<int> &B)
{
	vector< vector<int> >::const_iterator it;
	vector<int>::const_iterator edg;
	vector<int>::const_iterator Nod;
	B.clear();
	for(Nod=A.begin(); Nod !=A.end(); ++Nod )
	{
		int k=0;
		for(it=a.begin(); it != a.end(); ++it)
		{
			for (edg = it->begin(); edg != it->end(); edg++)
			{
				if(*edg==*Nod)
					k=k+1;				
			}						
		}	
		B.push_back(k);
	}
}

/* First_negighbors_List creates a list where the first element of the
rows is a node and the remaining elements are those nodes connected to 
it. The argument is a non-redundant list of edges g.*/

void First_Neighbors_List (vector<vector<int> > g, vector<vector<int> > &b)
{
	vector<int> A;
	b.clear();
	Nodes(g, A);
	int i,j;

	for(i=0; i<A.size(); ++i)
	{
		b.push_back(vector <int>());
		
		b[i].push_back(A[i]);
		for(j=0; j<g.size(); ++j)
		{
            if(g[j][0]==A[i])
				b[i].push_back(g[j][1]);
			if(g[j][1]==A[i])
				b[i].push_back(g[j][0]);			
		}
	}
}

/* Connected_Component creates a list of nodes paired to a -1 if 
they belong to the connected component the first node belongs to 
and to a 0 if it doesn't.*/

void Connected_Component (const vector<vector<int> > g, vector<vector<int> > &b)
{
	vector<vector<int> > A; //Vector where we store the visits of the queue
	vector<vector<int> > B; //List of first neighbors
	vector<int> Queue; //Queue where we introduce neighbors and sistematically burn them
	
	b.clear();
	Queue.clear();
	A.clear();
	B.clear();
	int k;

	First_Neighbors_List (g, B); //Computing the first neighbors list
	if(B.size()==0) //If the list is empty, there is no component
	{
		b.clear();
	}
	else
	{
		for(int i=0; i<B.size(); i++)
		{
			A.push_back(vector <int>());
			A[i].push_back(B[i][0]);
			A[i].push_back(0);
		}
		Queue.push_back(B[0][0]); //Root of the traverse process
		int idx=0;
		while (Queue.size()>0)  
		{
			for(int j=1; j<B[idx].size(); j++)
			{	
				for(k=0; k<A.size(); k++)
				{
					if((A[k][0]==B[idx][j]) && (A[k][1]==0))
					{
						Queue.push_back(B[idx][j]); //Listing the neighbors of Queue[0] 
						A[k][1]=-1;
					}	
				}
			}
			Queue.erase(Queue.begin());//Removing the first element of the queue when all its neighbors are listed
			for(k=0;  k<B.size(); k++)
			{
				if(B[k][0]==Queue[0]) //Finding the new first element of the queue
				{
					idx=k;
					break;
				}
			}
		}
		b=A;//Output: node paired with -1 if it belongs to the component of B[0][0] and 0 otherwise
	}
}

/*Component_Structure creates a list of nodes (b) in which each node is 
paired with the  label  (0,1...,h)  of  the component  it belongs to */

void Cores::Component_Structure (const vector<vector<int> > g, vector<vector<int> > &b)
{
	vector<vector<int> > a, B;
	a.clear();
	b.clear();
	B.clear();
	a=g;
	int h=0; //Label of the component.
	int u=0;

	//We remove components using Connected_Component until there are no remaining nodes 
	//Each erased component is labeled with h
	while(a.size()>0) 
	{
		Connected_Component(a, B);	
		if(B.size()==0)
		{
			a.clear();
		}
		else
		{
			for(int k=0; k<B.size(); k++)
			{
				if(B[k][1]==-1)
				{
					b.push_back(vector<int>());
					b[u].push_back(B[k][0]);
					b[u].push_back(h);
					u=u+1;
					int j;
					j=0;
					while(j<a.size())
					{
						if((a[j][0]==B[k][0]) || (a[j][1]==B[k][0]))
						{
							a.erase(a.begin()+j);
						}
						else
							j=j+1;
					}
				}			
			}
			h=h+1;
			B.clear();
		}	
	}
}

/*****************CORE EXTRACTING FUNCTIONS*****************************/
/*void KThreshold(const vector<vector<int> > a, int K, vector<vector<int> > &b)*/
/*void Kcore(vector<vector<int> > a, int K, vector<vector<int> > &b)*/
/*void GKThreshold(const vector<vector<int> > a, int K, vector<vector<int> > &b)*/
/*void GKcore(vector<vector<int> > a, int K, vector<vector<int> > &b)*/
/*void M_TGraph (const vector<vector<int> > A, const vector<vector<int> > B, vector<vector<int> > &M_Threshold, int M)*/
/*void M_core (vector<vector<int> > A, vector<vector<int> > B, vector<vector<int> > &MCSubgraph, int M)*/
/*void Cores::Core_clusters (vector<vector<int> > g, int K, int M, vector<vector<int> > &CoreCluster)*/

/***Function Ktreshold: f(Original graph a, Threshold k, New graph b)
New graph b is the subgraph of all nodes of a having k>K-1**********/

void KThreshold(const vector<vector<int> > a, int K, vector<vector<int> > &b)
{
	if(a.size()>0)
	{
		vector<int> B, C, KT;
		KT.clear();
		b.clear();
		B.clear();
		C.clear();

		Nodes(a, B);
		Connectivity(a, B, C);
		vector< vector<int> >::const_iterator it;
		vector<int>::iterator  edg, Nod;
	
		int	k=0;
		for(Nod=B.begin(); Nod !=B.end(); ++Nod)
		{	
			if(C[k] > K-1)
				KT.push_back(*Nod);
			k=k+1;
		}
	
		k=0;
		int h=0;
		int i;
		for(it=a.begin(); it !=a.end(); ++it)
		{
			int a1=0;
			int a2=0;
			for(Nod=KT.begin(); Nod !=KT.end(); ++Nod)
			{	
				if(a[k][0]==*Nod)
					a1=1;
				if(a[k][1]==*Nod)
					a2=1;
				if(a1==1 && a2==1)
					break;				
			}
			if(a1==1 && a2==1)
			{
				b.push_back(vector <int>());			
				for(i=0; i<2; ++i)
				{				
					b[h].push_back(a[k][i]);						
				}
				h=h+1;
			}
			k=k+1;
		}
	}
}	

/***Function Kcore:  f(Original graph a, Threshold k, New graph b)
New graph b is the  subgraph  whose  minimal connectivity is k>K-1. 
Is obtained through iteration  of KTreshold***********************/

void Cores::Kcore(vector<vector<int> > a, int K, vector<vector<int> > &b)
{
	KThreshold(a,K,b);	
	while(a.size()> b.size())
	{
		a=b;
		KThreshold(a,K,b);
	}
}

/***Function Ktreshold: f(Original graph a, Threshold k, New graph b)
New  graph  b  is  the  subgraph of all nodes of a  having  k>K-1  OR 
those nodes connecting two or more nodes having k>K-1***************/

void GKThreshold(const vector<vector<int> > a, int K, vector<vector<int> > &b)
{
	if(a.size()>0)
	{
		vector<int> B, C, KT;
		KT.clear();
		b.clear();
		B.clear();
		C.clear();

		Nodes(a, B);
		Connectivity(a, B, C);
		vector< vector<int> >::const_iterator it;
		vector<int>::iterator  edg, Nod, con;
		int k;

		for(k=0; k< B.size(); k++)
		{
			if(C[k] > K-1)
				KT.push_back(B[k]);
			else
			{
				int y=0;
				int i=0;
				for(i=0; i< a.size(); ++i)
				{
					if(a[i][0]==B[k])
					{	int j=0;
						for(j=0; j<B.size(); ++j)
						{
							if(a[i][1]==B[j])
							{
								if(C[j]>K-1)
									y=y+1;
								if(y>1)
									break;
							}							
						}
					}
					if(a[i][1]==B[k])
					{	int j=0;
						for(j=0; j<B.size(); ++j)
						{
							if(a[i][0]==B[j])
							{
								if(C[j]>K-1)
									y=y+1;
								if(y>1)
									break;
							}							
						}
					}										
					if(y>1)
					{
						KT.push_back(B[k]);
						break;
					}
				}
			}
		}

		k=0;
		int h=0;
		int i;
		for(it=a.begin(); it !=a.end(); ++it)
		{
			int a1=0;
			int a2=0;
			for(Nod=KT.begin(); Nod !=KT.end(); ++Nod)
			{	
				if(a[k][0]==*Nod)
					a1=1;
				if(a[k][1]==*Nod)
					a2=1;
				if(a1==1 && a2==1)
					break;				
			}
			if(a1==1 && a2==1)
			{
				b.push_back(vector <int>());			
				for(i=0; i<2; ++i)
				{				
					b[h].push_back(a[k][i]);						
				}
				h=h+1;
			}
			k=k+1;
		}
	}
}	

/***Function GKcore:  f(Original graph a, Threshold k, New graph b)
New  graph b is the  subgraph  whose  nodes have connectivity k>K-1
OR, having k<K, connect  two  nodes  having  k>K-1.  It is obtained 
through iteration  of  GKTreshold*********************************/

void Cores::GKcore(vector<vector<int> > a, int K, vector<vector<int> > &b)
{
	GKThreshold(a,K,b);	
	while(a.size()> b.size())
	{
		a=b;
		GKThreshold(a,K,b);
	}
}

/* M_Tgraph generates the subgraph from graph A (with list of first neighbors B) 
containing  all  edges  participating in,  at  least,  M  triangles.  Arguments: 
A: edgelist, B, list of  first  neighbors, results stored in an graph, edge list 
form: M_Threshold subragph */

void M_TGraph (const vector<vector<int> > A, const vector<vector<int> > B, vector<vector<int> > &M_Threshold, int M)
{
	int v,w	,mv, T;
	v=-1;
	w=-1;
	M_Threshold.clear();
	T=0;
	for(int i=0; i<A.size(); i++)
	{
		for(int k=0; k<B.size(); k++)
		{
			if(A[i][0]==B[k][0])
			{
				v=k;
				break;
			}
		}
		for(int j=0; j<B.size(); j++)
		{
			if(A[i][1]==B[j][0])
			{
				w=j;
				break;
			}				
		}
		mv=0;
		for(int n=1; n<B[v].size(); n++)
		{
			for(int m=1; m<B[w].size(); m++)
			{
				if((B[v][n]==B[w][m]) && ((B[v][n]!=B[w][0]) && (B[w][m]!=B[v][0])))
					mv=mv+1;
			}
		}
				
		if(mv>M-1)
		{
			M_Threshold.push_back(vector <int>());
			M_Threshold[T].push_back(A[i][0]);
			M_Threshold[T].push_back(A[i][1]);
			T=T+1;
		}
	}
}

/*Function M_core: From a graph A (with list of first neighbors B) generates the maximal 
subgraph (MCSubgraph)  in  which every link participates,  at  least,  in  M  triangles.  
Arguments: A  (edge  list  of the  original  graph)  B (list of first neighbors of the original graph),
MCSubgrah: obtained subgraph, M: Threshold*/

void Cores::M_core (vector<vector<int> > A, vector<vector<int> > B, vector<vector<int> > &MCSubgraph, int M)
{
	vector<vector<int> > C;
	M_TGraph (A, B, MCSubgraph, M);
	while(A.size()> MCSubgraph.size())
	{	
		A=MCSubgraph;
		First_Neighbors_List (A, B);
		M_TGraph (A, B, MCSubgraph, M);
	}
}

/*Function Core_clusters: from a graph g performs the gKcore (K parameter) and then the M-core (M parameter).
Returns a list of edges CoreCluster of the derived subgrpah*/

void Cores::Core_clusters (vector<vector<int> > g, int K, int M, vector<vector<int> > &CoreCluster)
{
	vector<vector<int> > B, gK;
	Cores Core;
	
	CoreCluster.clear();
	B.clear();
	gK.clear();
	
	Core.GKcore(g,K, gK);
	First_Neighbors_List (gK, B);
	Core.M_core (gK, B, CoreCluster, M);	
}

/* MC_Core_Raw_Comp Computes the potential list of nodes that could belong to the Weak core after 
computing the GK-core (parameter K) and,
over the GKcore, the M-core (parameter M): Those who are not at the 
M-core or those of the M-core that are connected to some node outside the M-core. It outputs a 4dim array
(v,u, c(v),c(u)). If the node is not in the M-core, c(v)=-1*/

void MC_Core_Raw_Comp (const vector<vector<int> > g, int K, int M, vector<vector<int> > &gKcomp)
{
	vector<vector<int> > B, gK, MC, Comp;
	vector<int> NodgK, NodMC, diff;
	Cores Core;
	MC.clear();
	B.clear();
	Comp.clear();
	gKcomp.clear();

	Core.GKcore(g,K, gK);
	First_Neighbors_List (gK, B);
	Core.M_core (gK, B, MC, M);
	Core.Component_Structure (MC, Comp);	

	NodgK.clear();
	NodMC.clear();
	diff.clear();
	Nodes(gK, NodgK);
	Nodes(MC, NodMC);
	
	sort(NodgK.begin(), NodgK.end());
	sort(NodMC.begin(), NodMC.end());
	set_difference(NodgK.begin(), NodgK.end(), NodMC.begin(), NodMC.end(), inserter(diff, diff.begin()));
		
	int s=Comp.size();
	for(int i=0; i<diff.size(); i++)
	{
		Comp.push_back(vector<int>());
		Comp[i+s].push_back(diff[i]);
		Comp[i+s].push_back(-1);		
	}
	
	int h=0;
	for(int i=0; i<gK.size(); i++)
	{
		int c,d;
		c=-2;
		d=-2;
		for(int k=0; k<Comp.size(); k++)
		{
			if(gK[i][0]==Comp[k][0])
			{
				c=Comp[k][1];
				break;
			}
		}

		for(int k=0; k<Comp.size(); k++)
		{
			if(gK[i][1]==Comp[k][0])
			{
				d=Comp[k][1];
				break;
			}
		}

		if((c!=d) || (c==-1 || d==-1))
		{
			gKcomp.push_back(vector<int>());	
			gKcomp[h].push_back(gK[i][0]);
			gKcomp[h].push_back(gK[i][1]);
			gKcomp[h].push_back(c);
			gKcomp[h].push_back(d);
			h=h+1;
		}		
	}
}

/*List communities: From the intput 4dim array (u,v,c(u),c(v)) creates a list of first neighbours FN and a list 
with the same dimensions in which the node v in the position i,j of the FN is occupied by c(v)*/

void Lists_communities (const vector<vector<int> > Comps, vector<vector<int> > &FN, vector<vector<int> > &FNComp)
{
	vector<vector<int> > edges, NodComp;
	vector<int> Nod;
	edges.clear();
	
	FN.clear();
	FNComp.clear();
	
	for(int i=0; i< Comps.size(); i++)
	{
		edges.push_back(vector<int>());
		edges[i].push_back(Comps[i][0]);
		edges[i].push_back(Comps[i][1]);		
	}
	
	First_Neighbors_List (edges, FN);	
	Nodes(edges,Nod);
	
	//Array NodComp with pairs (v, c(v))
	
	for(int i=0; i<Nod.size(); i++)
	{
		NodComp.push_back(vector<int>());
		NodComp[i].push_back(Nod[i]);
		for(int k=0; k<Comps.size(); k++)
		{
			//cout<<Nod[i]<<" "<<Comps[k][0]<<endl;

			if(Nod[i]==Comps[k][0])
			{
				NodComp[i].push_back(Comps[k][2]);
				break;
			}
			if(Nod[i]==Comps[k][1])
			{
				NodComp[i].push_back(Comps[k][3]);
				break;
			}			
		}
	}
	
	//Array FNComp: for each position of the nodes in FN, we print the component they are part of
	for(int i=0; i<FN.size(); i++)
	{
		FNComp.push_back(vector<int>());
		for(int k=0; k<FN[i].size(); k++)
		{
			for(int h=0; h<NodComp.size(); h++)
			{
				if(FN[i][k]==NodComp[h][0])
				{
					FNComp[i].push_back(NodComp[h][1]);
					break;
				}
			}
		}
	}
}

/*Intra community: Finds the nodes outside the  K-M core cluster which connect clusters of the K-M core. Input:
FN and FNComp, computed in Lists_communities. Output: The list of links connecting different k-M core clusters*/

void IntraCommunity (const vector<vector<int> >  FN, const vector<vector<int> >  FNComp, vector<vector<int> > &Intra)
{
	vector<vector<int> > A; //Vector where we store the visits of the queue
	vector<int> Nod, QueueV, QueueC, NodIntra;
	Intra.clear();
	Nod.clear();
	NodIntra.clear();

	for(int i=0; i<FN.size(); i++)
	{
		QueueV.clear();
		QueueC.clear();
		QueueV.push_back(FN[i][0]); //Root of the traverse process
		QueueC.push_back(FNComp[i][0]);
		int idx=i;
		A.clear();

		for(int m=0; m<FN.size(); m++)
		{
			A.push_back(vector <int>());
			A[m].push_back(FN[m][0]);
			A[m].push_back(0);
		}
		
		for(int k=0; k<A.size(); k++)
		{
			if(A[k][0]==FN[idx][0])
				A[k][1]=-1;
		}

        while (QueueV.size()>0)
		{
			for(int j=1; j<FN[idx].size(); j++)
			{	
				for(int k=0; k<A.size(); k++)
				{
					if((A[k][0]==FN[idx][j]) && (A[k][1]==0))
					{
						QueueC.push_back(FNComp[idx][j]);
						A[k][1]=-1;
				
						if(FNComp[idx][j]==-1)
							QueueV.push_back(FN[idx][j]); //Listing the neighbors of Queue[0] 
					}
				}
			}
			
			QueueV.erase(QueueV.begin());
						
			//Removing the first element of the queue when all its neighbors are listed
			if(QueueV.size()>0)
			{
				for(int k=0;  k<FN.size(); k++)
				{
					if(FN[k][0]==QueueV[0]) //Finding the new first element of the queue
					{
						idx=k;
						break;
					}
				}
			}
		}

		int h=0;
		int c=-1;
		for(int k=0; k<QueueC.size(); k++)
		{
			if(QueueC[k]>-1 && QueueC[k]!=c)
			{	
				h=h+1;
				c=QueueC[k];				
			}
			if(h>1)
			{	
				NodIntra.push_back(FN[i][0]);
				break;
			}
		}
	}
		
	/*Work in progress*/
	//FirstNeighbours list to edge list
	vector<vector<int> > edges;
	edges.clear();
	int w=0;
	for(int i=0; i<FN.size(); i++)
	{	
		for(int k=0; k<FN[i].size();k++)
		{
			if(FN[i][k]>FN[i][0])
			{
				edges.push_back(vector<int>());
				edges[w].push_back(FN[i][0]);
				edges[w].push_back(FN[i][k]);
				w=w+1;
			}			
		}
	}
	
	int l=0;
	for(int i=0; i<edges.size(); i++)
	{
		int h=0;
		for(int k=0; k<NodIntra.size(); k++)
		{
			if(NodIntra[k]==edges[i][0])
				h=h+1;
		}
		
		for(int k=0; k<NodIntra.size(); k++)
		{	
			if(NodIntra[k]==edges[i][1])
				h=h+1;
		}
		if(h>1)
		{
			Intra.push_back(vector<int>());
			Intra[l].push_back(edges[i][0]);
			Intra[l].push_back(edges[i][1]);
			l=l+1;			
		}
	}
}

/* Function Weak_core: Computes the weak core of a graph. Input: Graph g, K, M. Output: List of edges corresponding to the weak core WeakC*/

void Cores::Weak_core(const vector<vector<int> > g, int K, int M, vector<vector<int> > &WeakC)
{
	vector<vector<int> > Weak_Raw, gKcomp, FN, FNComp;
	gKcomp.clear();
	Weak_Raw.clear();
	FN.clear();
	FNComp.clear();
	WeakC.clear();

	MC_Core_Raw_Comp (g, K,  M,gKcomp);
	Lists_communities (gKcomp, FN, FNComp);
	IntraCommunity (FN, FNComp, WeakC);
}

/*Function KCoreness: Computes the deepest K-core a node belongs to. Input: graph g, output, list of pairs node--coreness of the node*/

void Cores::KCoreness(const vector<vector<int> > g, vector<vector<int> > &KCor)
{
	vector<vector<int> > KC1,KC2;
	vector<int> Nod1,Nod2, diff;
	Cores Core;
	KCor.clear();
	KC1.clear();
	
	KC1=g;
	Nodes(KC1,Nod1);
	int K=2;
	int h=0;
	while(KC1.size()>0)
	{
		diff.clear();
		KC2.clear();
		Nod2.clear();
		Core.Kcore(KC1,K,KC2);
		Nodes(KC2, Nod2);
		
		sort(Nod1.begin(), Nod1.end());
		sort(Nod2.begin(), Nod2.end());
		set_difference(Nod1.begin(), Nod1.end(), Nod2.begin(), Nod2.end(), inserter(diff, diff.begin()));

		Nod1=Nod2;
		KC1=KC2;
		
		for(int i=0; i<diff.size(); i++)
		{
			KCor.push_back(vector<int>());
			KCor[h].push_back(diff[i]);
			KCor[h].push_back(K-1);
			h=h+1;
		}
		K=K+1;
	}	
}

/*Function GKCoreness: Computes the deepest GK-core a node belongs to. Input: graph g, output, list of pairs node--GKcoreness of the node*/

void Cores::GKCoreness(const vector<vector<int> > g, vector<vector<int> > &KCor)
{
	vector<vector<int> > KC1,KC2;
	vector<int> Nod1,Nod2, diff;
	Cores Core;
	KCor.clear();
	KC1.clear();
	
	KC1=g;
	Nodes(KC1,Nod1);
	int K=2;
	int h=0;
	while(KC1.size()>0)
	{
		diff.clear();
		KC2.clear();
		Nod2.clear();
		Core.GKcore(KC1,K,KC2);
		Nodes(KC2, Nod2);
		
		sort(Nod1.begin(), Nod1.end());
		sort(Nod2.begin(), Nod2.end());
		set_difference(Nod1.begin(), Nod1.end(), Nod2.begin(), Nod2.end(), inserter(diff, diff.begin()));

		Nod1=Nod2;
		KC1=KC2;
		
		for(int i=0; i<diff.size(); i++)
		{
			KCor.push_back(vector<int>());
			KCor[h].push_back(diff[i]);
			KCor[h].push_back(K-1);
			h=h+1;
		}
		K=K+1;
	}	
}

/*Function MCoreness: Computes the deepest M-core a node belongs to. Input: graph g, output, list of pairs node--Mcoreness of the node*/

void Cores::MCoreness(const vector<vector<int> > g, vector<vector<int> > &KCor)
{
	vector<vector<int> > KC1,KC2, B;
	vector<int> Nod1,Nod2, diff;
	Cores Core;
	KCor.clear();
	KC1.clear();
	
	KC1=g;
	Nodes(KC1,Nod1);
	int M=2;
	int h=0;
	while(KC1.size()>0)
	{
		diff.clear();
		KC2.clear();
		Nod2.clear();
		B.clear();
		First_Neighbors_List(KC1, B);
		Core.M_core(KC1,B,KC2, M);
		Nodes(KC2, Nod2);
		
		sort(Nod1.begin(), Nod1.end());
		sort(Nod2.begin(), Nod2.end());
		set_difference(Nod1.begin(), Nod1.end(), Nod2.begin(), Nod2.end(), inserter(diff, diff.begin()));

		Nod1=Nod2;
		KC1=KC2;
		
		for(int i=0; i<diff.size(); i++)
		{
			KCor.push_back(vector<int>());
			KCor[h].push_back(diff[i]);
			KCor[h].push_back(M-1);
			h=h+1;
		}
		M=M+1;
	}	
}

/**************************************************/
/***************	     END		***************/
/**************************************************/

}
#endif
