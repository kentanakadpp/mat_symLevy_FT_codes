function [x,y,z] = phi(t, alpha, beta) 
    if(t==0) 
        x = 1/(2+alpha+beta);
    else
        x = t/(1-exp(-2*t - alpha*(1-exp(-t)) - beta*(exp(t)-1)));
    end
    y = phid(t, alpha, beta);
    z = x - t;
end

function y = phid(t, alpha, beta)
    if(t==0) 
        y = (4+alpha^2+3*beta+beta^2+alpha*(5+2*beta))/(2*(2+alpha+beta)^2);
    else
        tmpexp = exp(-2*t - alpha*(1-exp(-t)) - beta*(exp(t)-1));
        one_m_tmpexp_2 = (1 - tmpexp)^2;
        if(one_m_tmpexp_2 > realmax)
            y = 0.0;
        else        
            y = (1 - tmpexp * (1 + t * (2 + alpha * exp(-t) + beta * exp(t)))) / one_m_tmpexp_2;
        end
    end
end
