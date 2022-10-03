function dy = Marmottant_model(t,y,~,pdriv)

calling_obj = evalin('base','calling_obj');
c = calling_obj.c;
w = calling_obj.w;
P0 = calling_obj.P0;
R0 = calling_obj.R0;
rhoL = calling_obj.rhoL;
xi = calling_obj.xi;
muL = calling_obj.muL;
sigma0 = calling_obj.sigma0;
sigma_w = calling_obj.sigma_w;
gamma = calling_obj.gamma;
kappaS = calling_obj.kappaS;
Rbuckling = calling_obj.Rbuckling;
Rruptured = calling_obj.Rruptured;
broken = calling_obj.broken;

R = y(1);
Rdot = y(2);

F = griddedInterpolant(pdriv.t,pdriv.hyd);
P = F(t);

sigma = 0;

if broken == 0
    if R<Rbuckling
        sigma = 0;
    elseif R>=Rbuckling && R<=Rruptured
        sigma = xi*(R*R/Rbuckling/Rbuckling-1);
    elseif R>Rruptured
        sigma = sigma_w;
        broken = 1;
    end
else
     if R<Rbuckling
        sigma = 0;
     else
        sigma = xi*(R*R/Rbuckling/Rbuckling-1);
     end
     if sigma>sigma_w
         sigma = sigma_w;
     end
end    

term1 = ((P0 + 2*sigma0/R0)*((R/R0)^(-3*gamma))*(1-3*gamma*Rdot/c))/rhoL/R;
term2 = (-P0 -2*sigma/R -4*muL*Rdot/R -4*kappaS*Rdot/R/R -P)/rhoL/R;
term3 = -3*Rdot*Rdot/2/R;

Rdotdot = (term1+term2+term3); 

dy = [  y(2); Rdotdot ];
calling_obj.broken = broken;
assignin('base','calling_obj',calling_obj);
end