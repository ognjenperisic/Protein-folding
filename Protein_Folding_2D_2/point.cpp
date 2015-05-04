#include "stdafx.h"
#include "point.h"

using namespace std;

int point::max_counts    = 0;
int point::cooling_step  = 0;
int point::splot         = 0;
float point::minimum     = -13;

vector<int> point::aa;
    

int point::set_variables(int max_steps, int cool_step)
{
    point::max_counts   = max_steps;
    point::cooling_step = point::max_counts/cool_step;
    point::splot        = point::max_counts/100;

    return 0;
}

point::point()
{
}

point::point(int size)
{
    int i;

    point::aa.resize(chain_size);

    coord.resize(chain_size);
    for(i=0; i<chain_size; i++)
        coord[i].resize(2);
}

int point::fill_chain(int r)
{
    FILE* fin;
    int type;
    int i = 0;
    int x, y;
    
    point::aa.resize(point::chain_size);

    coord.resize(chain_size);
    for(i=0; i<chain_size; i++)
        coord[i].resize(2);


    fin = fopen("chain.txt", "r");

    i = 0;
        
    while (fscanf(fin, "%d", &type)!=EOF) {
        point::aa[i] = type;        
        coord[i][0] = i + 1 + 100;        
        coord[i][1] = chain_size/2 + 100;
        i++;
    };
    printf("\nFinished defining linear chain.txt\n\n");

    fclose(fin);
    	
    if (r==1) {
        i = 0;
        fin = fopen("circular_chain.txt", "r");
        
        while (fscanf(fin, "%d %d", &x, &y)!=EOF) {
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
    
    for (int i = 0; i < chain_size; i++) {
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
    int ok = 1;	

    vector<vector<int>>::iterator pr1;
    vector<vector<int>>::iterator pr2;    

                
    for(pr1 = this->coord.begin(); (pr1 != this->coord.end()) & (ok==1); pr1++)
        for(pr2 = this->coord.begin(); (pr2 != this->coord.end()) & (ok==1); pr2++) {
            if (pr1 != pr2) {
                if((pr1->at(0) == pr2->at(0)) & (pr1->at(1) == pr2->at(1))) {
                    ok = 0;                    
                }
            }
        }
        
    return ok;
}

int point::energy()
{
    float energy = 0;
    int i, j, dx, dy;
    int temp_aa_i, temp_aa_j;

    for (i = 0; i<chain_size; i++)
        for (j = 0; j<chain_size; j++) {
            if ((i!=j) & (abs(i - j)>2)) {
                dx = abs(coord[i][0] - coord[j][0]);
                dy = abs(coord[i][1] - coord[j][1]);


                if ((dx + dy) == 1) {
                    temp_aa_i = point::aa[i];
                    temp_aa_j = point::aa[j];

                    energy = energy - point::aa[i]*point::aa[j];
                }
            }
        }
        
        energy = energy*0.5;

    return energy;
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
    
    point hlp(36);
        

    for (int i = 0; i<chain_size; i++) {
        hlp.coord[i][0] = this->coord[i][0];
        hlp.coord[i][1] = this->coord[i][1];
        
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
            for (int j = 0; j<chain_size; j++) {
                this->coord[j][0] = hlp.coord[j][0];
                this->coord[j][1] = hlp.coord[j][1];
                
            }
            sgn= -sgn;
            
            for (int j = 0; j<k; j++) {
                dst1 = coord[k][0] - coord[j][0];
                dst2 = coord[k][1] - coord[j][1];
                coord[j][1] = coord[k][1] + sgn*dst1;
                coord[j][0] = coord[k][0] - sgn*dst2;
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
                for (i = 0; i<k; i++) {
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
            else {
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
    
    while ((cord<2) & (!cr))
    {
        if ((coord[k][cord] == coord[k - 1*pos][cord]) & ((coord[k][cord] + 1*pos) == coord[k + 1*pos][cord]) & (coord[k + 2*pos][cord]==coord[k+1*pos][cord])&
            ((coord[k - 2*pos][cord] + pos) == coord[k][cord]) & ((coord[k + pos][cord] + pos) == coord[k + 3*pos][cord])
            & (coord[k - pos][abs(1 - cord)] == coord[k + 2*pos][abs(1 - cord)]))
        {
            cr = true;
        }
        else {
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
    
    if (cr){
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

