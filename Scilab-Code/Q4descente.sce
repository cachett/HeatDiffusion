// Script question 4, fonction de descente calculant les coefficients de Z dans l'équation LZ=Y

clc;

function Z = descente(Ldiag, Linf, Y)
    //Nombre d'élément du vecteur ligne
    indiceMax = size(Ldiag, "c")
    //Calcul des coefficients, Z est un vecteur ligne 
    for i=1:indiceMax
        if i==1 then
            //Z(1) = Y(1)/L(1,1) 
            Z(1,i)=Y(i)/Ldiag(i)
        else
            Z(1,i)=(Y(i)-Z(1,i-1)*Linf(i-1))/Ldiag(i)
        end
    end
endfunction
