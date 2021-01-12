//
//  main.cpp
//  Git_Hub_Cores
//
//  Created by bernat on 12.01.21.
//  Copyright © 2021 Bernat. All rights reserved.
//
/******************************************************
 **   CORE EXTRACTING Functions                       **
 **   Bernat Corominas-Murtra                         **
 **   Last revision: 16.11.2017                       **
 **   Vienna Complexity Science Hub                     **
 **   Medical University of Vienna                     **
 **   Section for the Science of Complex Systems         **
 **     Includes:                                         **
 **        Graph Component structure                     **
 **        K-core                                         **
 **        Generalized K-core                             **
 **        M-core                                         **
 **        Weak Core                                     **
 **        Core Clusters (MCore of GKcore)                 **
 **        KCoreness                                     **
 **        GKCoreeness                                     **
 **        Mcoreness                                     **
 **                                                     **
 ******************************************************/

/******************************************************
 * Bibliography
 *
 *    K-core:
 *   K-core organization of complex networks
 *   S.N. Dorogovtsev, A.V. Goltsev, J.F.F. Mendes
 *    Phys. Rev. Lett. 96, 040601 (2006)
 *
 *   Generalized K-core:
 *   Detection of the Elite Structure in a Virtual Multiplex Social System by Means of a Generalised K-Core
 *   B. Corominas-Murtra, B. Fuchs, S. Thurner PloS one 9 (12), e112606 (2014)
 *
 *   M-core:
 *    Deciphering the global organization of clustering in real complex networks
 *    Pol Colomer-de-Simón, M. Ángeles Serrano, Mariano G. Beiró, J. Ignacio Alvarez-Hamelin & Marián Boguñá
 *    Scientific Reports 3, Article number: 2517 (2013)
 *
 *    Weak Core & Core_Cluster:
 *    The weak core and the structure of elites in social multiplex networks
 *    B. Corominas-Murtra and S. Thurner
 *    Springer Complexity: Interconnected Networks, Chapter 10, (2015)
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


//Input graph: Undirected graph. Data structure: list of pairs with integer ID.
//Example:
/*23 34
 45 78
 8     96
 65 8
 ...*/
//NO DUPLICATED EDGES ALLOWED: If e.g., 45 76 is in the list, 76 45 should NOT be in the list


#define K 2
#define M 1

int main (int argc, char * const argv[]) {
    
    //Salute
    cout<<endl;
    cout<<"********Cores*******"<<endl;
    cout<<endl;
    
    //Declarations
    vector<vector<int> > g, KC, GKC, MC, First_N, WeakC, CoreClus, KCore_seq, GKCore_seq, MCore_seq;
    Cores_BCM::Cores Core;
    
    
    //Input file names
    ifstream InputNet("Net_Example_1.txt");//Put the name of the input file here
    
    //Output file names
    //Output subgraphs: List of pairs of nodes
    ofstream KC_file("Kcore.txt");
    ofstream GKC_file("GKcore.txt");
    ofstream MC_file("Mcore.txt");
    ofstream WKC_file("Weak_core.txt");
    ofstream Cluster_Cores_file("Cluster_Core.txt");
    
    //Output list: Node ID---Coreness of the Node
    ofstream KCoreness_sequence_file("Coreness_sequence.txt");
    ofstream GKCoreness_sequence_file("GCoreness_sequence.txt");
    ofstream MCoreness_sequence_file("MCoreness_sequence.txt");
    
    //Getting data and storing it in g and First_N
    Cores_BCM::getGraph(InputNet, g);
    
    cout<<"Edges...: "<<g.size()<<endl;
    cout<<"First_neighbours..."<<endl;
    Cores_BCM::First_Neighbors_List (g, First_N); //only needed for the computation of the M-core
    cout<<"Nodes...: "<< First_N.size()<<endl;
    
    //Computing the cores
    cout<<"Kcore..."<<endl;
    Core.Kcore(g,K,KC);
    cout<<"GKcore..."<<endl;
    Core.GKcore(g,K,GKC);
    cout<<"Mcore..."<<endl;
    Core.M_core (g, First_N, MC, M);
    cout<<"WeaKcore..."<<endl;
    Core.Weak_core(g,  K,  M, WeakC);
    cout<<"CoreClusters..."<<endl;
    Core.Core_clusters(g,  K,  M, CoreClus);
    cout<<"Kcoreness..."<<endl;
    Core.KCoreness(g, KCore_seq);
    cout<<"GKcoreness..."<<endl;
    Core.GKCoreness(g,GKCore_seq);// Caution: time expensive for massive graphs
    cout<<"Mcoreness..."<<endl;
    Core.MCoreness(g,MCore_seq);
    
    //Printing results to files
    Cores_BCM::printGraph(KC_file,KC);
    Cores_BCM::printGraph(GKC_file,GKC);
    Cores_BCM::printGraph(MC_file,MC);
    Cores_BCM::printGraph(WKC_file,WeakC);
    Cores_BCM::printGraph(Cluster_Cores_file, CoreClus);
    Cores_BCM::printGraph(KCoreness_sequence_file, KCore_seq);
    Cores_BCM::printGraph(GKCoreness_sequence_file, GKCore_seq);
    Cores_BCM::printGraph(MCoreness_sequence_file, MCore_seq);
    
    return 0;
}

