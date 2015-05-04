// Protein_Folding_2D_2.cpp 
// Written by Ognjen Perisic in 2004
// Modified by Ognjen Perisic in 2012/2015

#include "stdafx.h"
#include "point.h"

using namespace std;




int main(int argc, char* argv[])
{
    int t = 0;
    double i;
    int rot;
    point a;
    point b(36);
    point c(36);
    double temp = 2;
    double rnd;
    double vl;
    
    FILE* fen;
    FILE* data;
    int min_e = 10;

    int pivot_test = 0;
    int residue;
    
    srand( (unsigned)time( NULL ) );
       
    if (argc < 3){
        printf("2 arguments required,\n 1: type (0==linear, 1==circular, 2==rectangular),\n 2: number of steps\n\n");
        exit(1);
    }
   


    a.fill_chain(atoi(argv[1]));

    a.set_variables(atoi(argv[2]), 200);

    i = rand()/(RAND_MAX /point::chain_size + 1);

    if (i <= (point::chain_size/2))
        i = 0;
    else
        i = point::chain_size - 1;

    a.end_rotation(i);


    fen = fopen("energy.txt", "w");

    while (t < point::max_counts) {
        if (fmod((float) t, (float) point::splot) == 0) {
            printf("%3d %% done!\n", int(100*t/point::max_counts));
        }

        
        b = a;
                    
        i = ((double)rand() / ((double)RAND_MAX + 1) * 11); 

        
        // 1-residue movement or a 2-residue movement
        if ((i <= 3)&(i > 1)) {
        
            i = rand()/(RAND_MAX /point::chain_size + 1);
            
            if ((i != 0) & (i != (point::chain_size - 1))) {        
                rot = b.corner_move((int) i);        
            }
            else {
                b.end_rotation((int)i);
            }
        }
        
        else if (i>3) {
            i = rand()/(RAND_MAX /point::chain_size + 1);

            while((i<=2)|(i>=(point::chain_size-3))) {
                i = rand()/(RAND_MAX /point::chain_size + 1);
            }
                        
            b.crank((int) i);
                        
        }
        else {            
            c = b;

            i = rand()/(RAND_MAX /point::chain_size + 1);
                        
            while (i>=(point::chain_size-2)|(i==0)) {
                i = rand()/(RAND_MAX /point::chain_size + 1);
            }

           
            pivot_test = c.pivot((int) i);
            
            if (pivot_test>0) {
                b = c;
            }
           
        }
        


        if (b.energy() < point::minimum) {
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

        if (b.energy() <= a.energy()) {
            a = b;  // switch configurations     
        }
        else {
            rnd = (double)((double)rand() / ((double)RAND_MAX + 1) * 1);
            vl = exp((double(-(b.energy() - a.energy())))/temp);

            if (rnd < vl) {
                a = b;
            }
        }
        
        t++;  // increase the time step

        if((t % 1000)==0)
            cout << "t = " << t << endl;

        if ((fmod((float) t, (float) point::cooling_step)==0) & (t > point::burn))
            temp = temp*0.99;

        if ((fmod((float)t, (float) floor(point::max_counts*0.0025)) == 0) & (t > point::burn)) {
            fprintf(fen,"%10d  %3d\n",t,a.energy());
        }
    } // while (t < point::max_counts)

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
        
    return 0;
}