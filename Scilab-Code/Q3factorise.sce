// Script question 3, factorisation pour la méthode de Cholesky dans le cas tridiagonal

clc ; 


function [Ldiag, Linf] = factorise(Mdiag, Minf)
    //fonction qui renvoie deux vecteurs : 
    // linf contenant la sous diagonale de M de taille n-1 
    // et ldiag contenant la diagonale de M de taille n 
    
    //En Scilab il n'est pas nécessaire d'initialiser Ldiag et Linf. 
    indiceMaxDiag = size(Mdiag,"c")
    //Avec la méthode de cholesly, on calcule successivement les colonnes 
    for i=1:indiceMaxDiag
        if i==1 then //Initialisation du premier coef de la diagonale
            Ldiag(1,i)=sqrt(Mdiag(i))
        else
            Ldiag(1,i)=sqrt(Mdiag(i)-Linf(1,i-1).^2)
        end
        // On ne prends pas le dernier élement de la matrice M, Linf est de taille n-1 
        if i ~= indiceMaxDiag then
           Linf(1,i)=Minf(i)/Ldiag(1,i)
        end
    end
endfunction
