cd c:\something
clear all
close all
clc
fig = 0;

number_of_residues = 36;
radius = 8;

x = zeros(1, number_of_residues);
y = x;



x(1) = floor(radius*sin(0*2*pi/number_of_residues)) + 100;
y(1) = floor(radius*cos(0*2*pi/number_of_residues)) + 100;

for i = 2:number_of_residues
    
    s = sin(i*2*pi/number_of_residues);
    c = cos(i*2*pi/number_of_residues);
    if s*c<0
        if s<0
            dx = 1;
            dy = 1;
        else
            dx = -1;
            dy = -1;
        end            
    else
        if s>0
            dx = 1;
            dy = -1;
        else
            dx = -1;
            dy = 1;
        end
    end
   
    pick = rand(1);
    if (pick <= 0.5)
        x(i) = x(i-1) + dx;
        y(i) = y(i-1);
    else
        x(i) = x(i-1);
        y(i) = y(i-1) + dy;
    end       
   
end

fig = fig + 1;
figure(fig)
plot(x,y,'o-')
axis square
grid on

fid = fopen('circular_chain.txt', 'w');
for i=1:number_of_residues
   fprintf(fid,'%6d  %6d\n', x(i), y(i));
end

fclose(fid);

if(3>1)
   x(1:9)   = 1:9;
   x(10:18) = 9;
   x(19:27) = 8:-1:0;
   x(28:36) = 0;
   
   y(1:9)   = 9;
   y(10:18) = 8:-1:0;
   y(19:27) = 0;
   y(28:36) = 1:9;
   
   a = 1:36;
end




fig = fig + 1;
figure(fig)
plot(x,y,'o-')
axis square
axis([-2 10 -2 10])
grid on


fid = fopen('rectangular_chain.txt', 'w');
for i=1:number_of_residues
   fprintf(fid,'%6d  %6d\n', x(i), y(i));
end

fclose(fid);
