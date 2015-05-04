cd c:\something
clear all
close all
clc
fig = 0;


%load conf13d3.txt
load conf.txt
%conf = conf13d3
h = 0;
t = 1;


for i = 1:length(conf)
   if (conf(i, 3) == 1)
      h(t, 1) = conf(i, 1);
      h(t, 2) = conf(i, 2);
      t = t + 1;
   end      
end

fig = fig + 1;
figure(fig)
plot(conf(:, 1), conf(:, 2), '.')
hold on
plot(h(:, 1), h(:, 2), 'ro')
plot(h(:, 1), h(:, 2), 'r.')
plot(conf(:, 1), conf(:, 2), '-')
plot(h(:, 1), h(:, 2), 'ro')
plot(h(:, 1), h(:, 2), 'r.')
legend('hydrophilic residues', 'hydrophobic residues')


dst = 0.5*(max(conf(:, 2)) - min(conf(:, 2)));
dstx = max(conf(:, 1)) - min(conf(:, 1));
dsty = max(conf(:, 2)) - min(conf(:, 2));
mxr = 0.8*max(dstx, dsty);
xm = (max(conf(:, 1)) + min(conf(:, 1)))/2;
ym = (max(conf(:, 2)) + min(conf(:, 2)))/2;
axis([round(xm - mxr) round(xm + mxr) round(ym - mxr) round(ym + mxr)])


load energy.txt
fig = fig + 1;
figure(fig)
plot(energy(:,1), energy(:,2))
axis([-inf inf (min(energy(:,2))-1) (max(energy(:,2))+1)])
xlabel('Steps (number of sampling steps)')
ylabel('Energy')