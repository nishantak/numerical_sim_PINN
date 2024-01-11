function c1(N, cfl, problem, scheme)

[xmin,xmax,Tf,x,u] = initialize(N, problem);
umin = min(u); umax = max(u);

dx = (xmax - xmin)/(N-1);

res = zeros(N,1);

t = 0;
it= 0;
while t < Tf
   dt = cfl * dx / max(abs(DFLUX(u)));
   lambda = dt/dx;
   res(:) = 0;
   for j=1:N-1
      flux = num_flux( scheme, lambda, u(j), u(j+1) );
      res(j)   = res(j)   + flux;
      res(j+1) = res(j+1) - flux;
   end
   u(2:N-1) = u(2:N-1) - lambda * res(2:N-1);
   t = t + dt;
   it= it+ 1;
   fprintf(1,'it = %d, t = %e, range= %e, %e\n', it, t, min(u), max(u))
   if mod(it,5) == 0
      [xe,ue] = exact_solution(problem, xmin, xmax, t);
      plot(xe,ue,x,u,'o','LineWidth',1.5);
      axis([xmin xmax umin umax])
      pause(1)
   end
end

function [xmin,xmax,Tf,x,u] = initialize(N, problem)

u = zeros(N,1);

if problem==1
   Tf   = 1.5;
   xmin = -2.0;
   xmax =  2.0;
   x = linspace(xmin, xmax, N);
   for j=1:N
      if x(j) <= -0.25
         u(j) = 1.0;
      elseif x(j) > -0.25 && x(j) <= 0.25
         u(j) = 3/4 - x(j);
      else
         u(j) = 0.5;
      end
   end
elseif problem==2
   Tf   = 1.0;
   xmin = -2.0;
   xmax =  2.0;
   ul = -0.5;
   ur = +1.0;
   x = linspace(xmin, xmax, N);
   for j=1:N
      if x(j) <= 0.0
         u(j) = ul;
      else
         u(j) = ur;
      end
   end
elseif problem==3
   Tf   = 1.25;
   xmin = 0.0;
   xmax = 3.0;
   ul = 1.0;
   ur = 1.0;
   x = linspace(xmin, xmax, N);
   for j=1:N
      if x(j) < 0 || x(j) > 1
         u(j) = 1;
      else
         u(j) = 1 + sin(pi*x(j));
      end
   end
else
   fprintf(1, 'Unknown problem\n');
   pause
end

function [xe,ue] = exact_solution(problem, xmin, xmax, t)
xe = linspace(xmin, xmax, 1000);
ue = zeros(size(xe));

if problem==1
   tc = 1;
   xc = 3/4;
   s  = 3/4;
   for j=1:length(xe)
      if t <= tc
         if xe(j) < t - 1/4
            ue(j) = 1.0;
         elseif xe(j) > 0.5*t + 1/4
            ue(j) = 0.5;
         else
            ue(j) = 3/4 - (xe(j) - 0.75*t)/(1-t);
         end
      else
         if xe(j) < xc + (t-tc)*s
            ue(j) = 1.0;
         else
            ue(j) = 0.5;
         end
      end
   end
elseif problem==2
   ul = -0.5;
   ur = +1.0;
   for j=1:length(xe)
      if xe(j) < ul*t
         ue(j) = ul;
      elseif xe(j) > ur*t
         ue(j) = ur;
      else
         ue(j) = xe(j)/t;
      end
   end
elseif problem==3
   ue = ones(size(xe));
else
   fprintf(1,'Unknown problem\n');
   pause
end

function flux = FLUX(u)
flux = 0.5*u^2;

function dflux = DFLUX(u)
dflux = u;

function a = ADFLUX(u,v)
if abs(u-v) > 1.e-14
   a = (FLUX(u) - FLUX(v))/(u-v);
else
   a = DFLUX(u);
end

function flux = num_flux(scheme, lambda, u, v)

if scheme==1
   flux = 0.5*(FLUX(u) + FLUX(v)) - (0.5/lambda)*(v - u);
elseif scheme==2
   a = ADFLUX(u,v);
   flux = 0.5*(FLUX(u) + FLUX(v)) - (0.5*lambda*a)*(FLUX(v) - FLUX(u));
elseif scheme==3
   a = ADFLUX(u,v);
   flux = 0.5*(FLUX(u) + FLUX(v)) - 0.5*abs(a)*(v - u);
elseif scheme==4
   if u < v
      u1 = max(0,u);
      u2 = min(0,v);
      flux = max( FLUX(u1), FLUX(u2) );
   else
      flux = max( FLUX(u), FLUX(v) );
   end
elseif scheme==5
elseif scheme==6
   delta = 0.5*max(0, abs(v-u));
   a = abs(ADFLUX(u,v));
   if a < delta
      a = delta;
   end
   flux = 0.5*(FLUX(u) + FLUX(v)) - 0.5*a*(v - u);
else
   fprintf(1, 'Unknown flux scheme %d\n', scheme);
   pause
end
