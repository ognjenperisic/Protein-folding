// Point header file

#ifndef POINT_H
#define POINT_H

#include <fstream>
#include <direct.h>
#include <conio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>
#include <iostream>
#include <vector>

using namespace std;

class point
{

 public:

     point();
     point(int);     
     int set_variables(int, int);
     int fill_chain(int);
     int rotation();
     int initial_rotation();
     int end_rotation(int k = -1);
     int corner_move(int k);
     int pivot(int k);
     int check_structure(int k, int x, int y);
     int check();
     int crank(int k);
     int energy();
     int write_chain();
     int check_l(int r, int b);

     //point operator=(const point & a);
     

     static int max_counts;
     static int cooling_step;
     static int splot; 
     static float minimum;
     static vector<int> aa;

     const static int chain_size = 36;
     const static int burn = 1000;
     
 


//private:
     //int aa[chain_size];
     //vector<int> aa(chain_size);
     
     //int coord[chain_size][2];
     vector<vector <int>> coord;
      
};


#endif

