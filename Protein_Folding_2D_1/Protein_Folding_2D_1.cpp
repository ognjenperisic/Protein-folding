// Protein_Folding_2D_1.cpp 
// Written by Ognjen Perisic in 2004
// Modified by Ognjen Perisic in 2012/2015

#include "stdafx.h"
#include <fstream>
#include <direct.h>
#include <conio.h>
#include <stdlib.h>
#include <time.h>
#include <cmath>

using namespace std;

#define work

const int chain_size = 36;
const int burn = 1000;
const float minimum =-13;

class point
{
  public:
      int aa[chain_size];
      int coord[chain_size][2];
      int fill_chain(int);
      int rotation();
      int initial_rotation();
      int end_rotation(int k = -1);
      int corner_move(int k);
      int pivot(int k);
      int check_structure(int k, int x, int y);
      int check();
      int crank(int i);
      int energy();
      int write_chain();
      int check_l(int r, int b);

      static int max_counts;
      static int cooling_step;
      static int splot;
};


int point::max_counts = 0;
int point::cooling_step = 0;
int point::splot = 0;

int point::write_chain()
{
    FILE* fout;

    fout=fopen("conf.txt","w");
    
    for (int i = 0; i < chain_size; i++)
        fprintf(fout,"%3d %3d  %3d\n", coord[i][0], coord[i][1], aa[i]);
    
    fclose(fout);
    return 0;
}


int point::check_l(int r, int b)
{
    int k = 0;
    int ret = -1;
    int max = -50000, min = 10000;
    int rt;
    
    for (int i = 0; i < chain_size; i++){
        if (coord[i][r] >= max)
            max = coord[i][r];
        
        if (coord[i][r] <= min)
            min = coord[i][r];		
    }
    
    rt = max - min;
    
    if (rt <= 10) {
        printf("oo\n");
    }    
    
    return rt;
}


int point::check()
{
    int i, j;
    int ok = 1;	
    
    for (i = 0; i < chain_size; i++)
        for (j = 0; j < chain_size; j++) {
            if (i != j) {
                if ((coord[i][0]==coord[j][0]) & (coord[i][1]==coord[j][1])) {
                    ok = 0;
                    //break;
                }
            }
        }
        
    return ok;
}

int point::energy()
{
    float energy = 0;
    int i, j, dx, dy;
    
    for (i = 0; i<chain_size; i++)
        for (j = 0; j<chain_size; j++) {
            if ((i!=j) & (abs(i - j)>2)) {
                dx = abs(coord[i][0] - coord[j][0]);
                dy = abs(coord[i][1] - coord[j][1]);
                if ((dx + dy) == 1) {
                    energy = energy - aa[i]*aa[j];
                }
            }
        }
        
        energy = energy*0.5;
    return energy;
}


int point::fill_chain(int r)
{
    FILE* fin;
    int type;
    int i = 0;
    int x, y, a;
    
    fin = fopen("chain.txt", "r");
    while (fscanf(fin, "%d", &type)!=EOF) {
        aa[i] = type;
        coord[i][0] = i + 1 + 100;
        coord[i][1] = chain_size/2 + 100;
        i++;
    };
    printf("\nFinished defining linear chain.txt\n\n");


    fclose(fin);
    	
    if (r == 1) {
        i = 0;
        fin = fopen("circular_chain.txt", "r");
        
        while (fscanf(fin, "%d %d", &x, &y) != EOF) {
            coord[i][0] = x + 100;
            coord[i][1] = y + 100;
            i++;
        };
        fclose(fin);

        printf("\nFinished reading circular chain.txt\n\n");
    }

    if (r == 2) {
        i = 0;
        fin = fopen("rectangular_chain.txt", "r");
        
        while (fscanf(fin, "%d %d", &x, &y) != EOF) {
            coord[i][0] = x + 100;
            coord[i][1] = y + 100;
            i++;
        };
        fclose(fin);

        printf("\nFinished reading rectangular chain.txt\n\n");
    }

    
    return 0;
}


int point::rotation()
{
    int i = rand() / (RAND_MAX / chain_size + 1);
    return 0;
}


int point::end_rotation(int k)
{
    int sgn;
    int nxt;
    int a;
    int help_x, help_y;
    int out = 0;
    int m = rand();
    
    m = rand();
    m = rand();
    int i = m/(RAND_MAX /10 + 1);
    
    if (k<0)
        if (i<5) {
            k = 0;
            nxt = 1;
        }
        else {
            k = chain_size - 1;
            nxt = -1;
        }
    else {
        if (k==0)
            nxt = 1;
        else
            nxt = -1;
    }
    
    i = rand()/(RAND_MAX /10 + 1);
    
    if (i <= 5)
        sgn = -1;  // one of two possible movements
    else
        sgn = 1;
    
    int d1 = abs(abs(coord[k][0] - coord[k + nxt][0]) - abs(coord[k + 2*nxt][0] - coord[k + nxt][0]));
    int d2 = abs(abs(coord[k][1] - coord[k + nxt][1]) - abs(coord[k + 2*nxt][1] - coord[k + nxt][1]));
    int option  = d1*d2;  // straight or L shape;  0 or 1
    int option2 = abs(coord[k + 2*nxt][0] - coord[k + nxt][0]); 

    // up-down or left-right L shape	
    if (!option) // straight
    {
        int t1 = coord[k][0];
        int t2 = coord[k][1];
        help_x = coord[k + nxt][0] + abs(option - 1)*sgn*(coord[k + nxt][1] - t2);
        help_y = coord[k + nxt][1] + abs(option - 1)*sgn*(coord[k + nxt][0] - t1);
        
        if (check_structure(k, help_x, help_y)) {
            coord[k][0] = help_x;
            coord[k][1] = help_y;
            out = 1;
        }
        else {
            sgn= -sgn;
            help_x = coord[k + nxt][0] + abs(option - 1)*sgn*(coord[k + nxt][1]-t2);
            help_y = coord[k + nxt][1] + abs(option - 1)*sgn*(coord[k + nxt][0]-t1);
            
            if (check_structure(k, help_x, help_y)) {
                coord[k][0] = help_x;
                coord[k][1] = help_y;
                out = 1;
            }
        }
    }
    else // L shape
    {
        a = (sgn + 1)/2;
        if (!option2) {
            help_x = coord[k + nxt][option2] + abs((a - option2))*(coord[k + nxt][option2] - coord[k][option2]);
            help_y = coord[k + nxt][abs(option2 - 1)] + abs((a - abs(option2 - 1)))*(coord[k + nxt][abs(option2 - 1)]-coord[k + 2*nxt][abs(option2 - 1)]);
        }
        else {
            help_x = coord[k + nxt][abs(option2 - 1)] + abs((a - abs(option2 - 1)))*(coord[k + nxt][abs(option2 - 1)]-coord[k + 2*nxt][abs(option2 - 1)]);
            help_y = coord[k + nxt][option2] + abs((a - option2))*(coord[k + nxt][option2] - coord[k][option2]);
        }
        
        if (check_structure(k, help_x, help_y)) {
            coord[k][0] = help_x;
            coord[k][1] = help_y;
			out = 1;
        }
        else  // if the first movement is not possible
        {
            sgn = -sgn;
            a = (sgn + 1)/2;
            if (!option2) {
                help_x = coord[k + nxt][option2] + abs((a - option2))*(coord[k + nxt][option2] - coord[k][option2]);
                help_y = coord[k + nxt][abs(option2 - 1)] + abs((a - abs(option2 - 1)))*(coord[k + nxt][abs(option2 - 1)] - coord[k + 2*nxt][abs(option2 - 1)]);
            }
            else {
                help_x = coord[k + nxt][abs(option2 - 1)] + abs((a - abs(option2 - 1)))*(coord[k + nxt][abs(option2 - 1)] - coord[k + 2*nxt][abs(option2 - 1)]);
                help_y = coord[k + nxt][option2] + abs((a - option2))*(coord[k + nxt][option2] - coord[k][option2]);
            }

			if (check_structure(k, help_x, help_y)) {
                coord[k][0] = help_x;
                coord[k][1] = help_y;
                
                out = 1;
            }
        }
    }
    
    return out;
}


int point::check_structure(int k, int x, int y)
{
    int check = 1;
    
    for (int i = 0; i<chain_size; i++) {
        if (i != k) {
            if ((coord[i][0] == x) & (coord[i][1] == y))
                check = 0;
        }
    }

	return check;
}


int point::corner_move(int k)
{
    int d1 = abs(abs(coord[k][0] - coord[k - 1][0]) - abs(coord[k + 1][0] - coord[k][0]));
    int d2 = abs(abs(coord[k][1] - coord[k - 1][1]) - abs(coord[k + 1][1] - coord[k][1]));
    int option = d1*d2;  // straight or L shape  0 or 1
    
    if (option) {
        int shp = abs(coord[k - 1][0] - coord[k][0]);
        shp = 1 - 2*shp;
        
        int help_x = coord[k + shp][0];
        int help_y = coord[k - shp][1];
        if (check_structure(k, help_x, help_y)) {
            coord[k][0] = help_x;
            coord[k][1] = help_y;
            
            return 1;
        }
        else
            return 0;						
    }
    else
        return 0;
}


int point::pivot(int k)
{
    
    int nxt = 1, i, dst1, dst2, sgn;
    double rnd;
    rnd = ((double)rand() / ((double)RAND_MAX + 1) * 11);
    point hlp;
    
    for (int i = 0; i<chain_size; i++) {
        hlp.coord[i][0] = this->coord[i][0];
        hlp.coord[i][1] = this->coord[i][1];
        hlp.aa[i] = this->aa[i];
    }
    
    int d1 = abs(abs(coord[k][0] - coord[k - nxt][0]) - abs(coord[k + nxt][0] - coord[k][0]));
    int d2 = abs(abs(coord[k][1] - coord[k - nxt][1]) - abs(coord[k + nxt][1] - coord[k][1]));
    int option = d1*d2;  // straight or L shape;  0 or 1
    int option2 = abs(coord[k + 2*nxt][0] - coord[k + nxt][0]);

    // up-down or left-right L shape
    if (!option) {
        rnd = ((double)rand() / ((double)RAND_MAX + 1) );
        if (rnd <= 0.5)
            sgn = 1;
        else
            sgn = -1;
        
        for (i = 0; i<k; i++) {
            dst1 = coord[k][0] - coord[i][0];
            dst2 = coord[k][1] - coord[i][1];
            coord[i][1] = coord[k][1] + sgn*dst1;
            coord[i][0] = coord[k][0] - sgn*dst2;
        }
    }
    
    if (check() & (!option)) {
        return 1;
    }
    else {
        if (!option) {
            for (int i = 0; i<chain_size; i++) {
                this->coord[i][0] = hlp.coord[i][0];
                this->coord[i][1] = hlp.coord[i][1];
                this->aa[i] = hlp.aa[i];
            }
            sgn= -sgn;
            
            for (i = 0; i<k; i++) {
                dst1 = coord[k][0] - coord[i][0];
                dst2 = coord[k][1] - coord[i][1];
                coord[i][1] = coord[k][1] + sgn*dst1;
                coord[i][0] = coord[k][0] - sgn*dst2;
            }
            
            if (check()) {
                return 1;
            }
            else
                return 0;
        }
    }
    
    // L shape
    if (option) {
        
        rnd = ((double)rand()/((double)RAND_MAX + 1) );
        
        if (rnd <= 0.5)
            sgn = 1;
        else
            sgn = -1;
        
        if (sgn == 1) {
            
            rnd = ((double)rand()/((double)RAND_MAX + 1) );
            if (rnd <= 0.5) {
                for (i = 0; i<k; i++){
                    
                     dst1 = coord[k][0] - coord[i][0];						
                     dst2 = coord[k][1] - coord[i][1];
                     coord[i][0] = coord[k][0] - dst1;
                     coord[i][1] = coord[k][1] + dst2;                    
                }

				if (check())
                    return 1;
                else
                    return 0;
            }
            else
            {
                for (i = 0; i<k; i++) {                    
                    dst1 = coord[k][0] - coord[i][0];
					dst2 = coord[k][1] - coord[i][1];
                    coord[i][0] = coord[k][0] + dst2;
                    coord[i][1] = coord[k][1] + dst1;                    
                }

                if (check())
                    return 1;
                else
                    return 0;
            }
        }
    }
}


int point::crank(int k)
{
    int pos = 1;
    int mov;
    int cord = 0;
    int help_1, help_2;
    bool cr = false;
    
    while ((cord<2) & (!cr)) {
        if ((coord[k][cord] == coord[k - 1*pos][cord]) & ((coord[k][cord] + 1*pos) == coord[k + 1*pos][cord]) & (coord[k + 2*pos][cord]==coord[k+1*pos][cord])&
            ((coord[k - 2*pos][cord] + pos) == coord[k][cord]) & ((coord[k + pos][cord] + pos) == coord[k + 3*pos][cord])
            & (coord[k - pos][abs(1 - cord)] == coord[k + 2*pos][abs(1 - cord)]))
        {
            cr = true;
        }
        else  {
            pos = -1;
            
            if ((coord[k][cord] == coord[k - 1*pos][cord]) & ((coord[k][cord] + 1*pos) == coord[k + 1*pos][cord]) & (coord[k + 2*pos][cord] == coord[k + 1*pos][cord])&
                ((coord[k - 2*pos][cord] + pos) == coord[k][cord]) & ((coord[k + pos][cord] + pos) == coord[k + 3*pos][cord])
                & (coord[k - pos][abs(1 - cord)] == coord[k + 2*pos][abs(1 - cord)]))
            {
                cr = true;
            }
        }
        cord++;
    
    }
    
    if (cr) {
        if (cord > 1)
            cord = 0;
        else
            cord = 1;
        
        if (cord) {
            mov = coord[k][1] - coord[k - pos][1];
            help_1 = coord[k][1] - 2*mov;
            help_2 = coord[k + pos][1] - 2*mov;
            
            if (check_structure(k, coord[k][0], help_1) & (check_structure(k - pos, coord[k + pos][0],help_2))) {
                coord[k][1]       = help_1;
                coord[k + pos][1] = help_2;
            }
        }
        else {
            mov = coord[k][0] - coord[k - pos][0];
            help_1 = coord[k][0] - 2*mov;
            help_2 = coord[k + pos][0] - 2*mov;

            if (check_structure(k, help_1, coord[k][1]) & (check_structure(k + pos, help_2, coord[k + pos][1]))) {
                coord[k][0]       = help_1;
                coord[k + pos][0] = help_2;
            }
        }
    }
    
    return 0;
}



int main(int argc, char* argv[])
{
    int t = 0;
    double i;
    int rot;
    point a;
    point b;
    point c;
    double temp = 2;
    double rnd;
    double vl;
    
    FILE* fen;
    FILE* data;
    int min_e = 10;
    
    srand( (unsigned)time( NULL ) );
    

    point::max_counts = 2000;

    if (argc < 3){
        printf("2 arguments required,\n 1: type (0==linear, 1==circular, 2==rectangular),\n 2: number of steps\n\n");
        exit(1);
    }

    point::max_counts = atoi(argv[2]);



    point::cooling_step = point::max_counts / 220;
    point::splot = point::max_counts / 100;


    
    a.fill_chain(atoi(argv[1]));

    i = rand()/(RAND_MAX /chain_size + 1);
    if (i <= (chain_size/2))
        i = 0;
    else
        i = chain_size - 1;
    
    a.end_rotation(i);
    fen = fopen("energy.txt", "w");
    
    while (t < point::max_counts)  {

        if (fmod((float) t, (float) point::splot) == 0) {
            printf("Done : %3d %%\n", int(100*t/point::max_counts));
        }
        b = a;
        
        i = ((double)rand() / ((double)RAND_MAX + 1) * 11); 

        // 1-residue movement or a 2-residue movement
        if ((i <= 3)&(i > 1)) {
            i = rand()/(RAND_MAX /chain_size + 1);

            if ((i != 0) & (i != (chain_size - 1))) {
                rot = b.corner_move(i);
            }
            else {
                b.end_rotation(i);
            }
        }
        
        else if (i>3) {
            i = rand()/(RAND_MAX /chain_size + 1);
            b.crank(i);
        }
        else {
            c = b;
            i = rand()/(RAND_MAX /chain_size + 1);

            if ((i > 0) & (i < chain_size))
                if (c.pivot(i))			
                    b = c;
        }
        
        if (b.energy() < minimum) {
            if (b.energy() < min_e) {
                min_e = b.energy();
                b.write_chain();
                printf("Energy = %d\n", b.energy());
                data = fopen("data.txt", "w");
                fprintf(data, "Energy :%3d\n", b.energy());
                fprintf(data, "Step   :%10d\n", t);
                fprintf(data, "Temp   :%f\n", temp);
                fclose(data);
                fprintf(fen, "%10d  %3d\n", t, a.energy());
            }
        }
        
        if (b.energy() <= a.energy())
            a = b;
        else {
            rnd = (double)((double)rand() / ((double)RAND_MAX + 1) * 1);
            vl = exp((double(-(b.energy() - a.energy())))/temp);
            if (rnd < vl) {
                a = b;
            }
        }
        
        t++;

        if ((fmod((float) t, (float) point::cooling_step)==0) & (t>burn))
            temp = temp*0.99;

        if ((fmod((float)t, (float)floor(point::max_counts*0.0025)) == 0) & (t>burn)) {
            fprintf(fen,"%10d  %3d\n",t,a.energy());
        }
    }

	if (a.energy() < min_e) {
        a.write_chain();
        printf("Energy = %d\n", a.energy());
        data = fopen("data.txt", "w");
        fprintf(data, "Energy :%3d\n", b.energy());
        fprintf(data, "Step   :%10d\n", t);
        fprintf(data, "Temp   :%f\n", temp);
        fclose(data);
        fprintf(fen, "%10d  %3d\n", t, a.energy());
        fclose(fen);
    }

    //getche();
    return 0;
}