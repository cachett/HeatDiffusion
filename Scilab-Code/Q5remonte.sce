// Script question 5, fonction de remonté calculant la solution X du système trans(L)X = Z

clc;


function X = remonte(Ldiag, Linf, Z)
    indiceMax = size(Ldiag, "c")
    //Calcul des coefficients, X est un vecteur ligne
    for i=indiceMax:-1:1
        if i==indiceMax then
            X(1,i)=Z(i)/Ldiag(i)
        else
            X(1,i)=(Z(i)-X(1,i+1)*Linf(i))/Ldiag(i)
        end
    end
endfunction
