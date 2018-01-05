//Script pour la Q7, résoud numériquement le système trouvé pour i = 10, i=50, i=100 et u0=1.


exec("Q3factorise.sce")
exec("Q4descente.sce")
exec("Q5remonte.sce")

function [Adiag, Ainf, B] = genereMatricesAB(nbrPoints, l)//NbrPoints prend en compte les bornes!! Donc n = nbrPoints-2
    // on génére les matrices A et B trouvés à la Q6, connaissant l'expression de C
    deltax = 2*l/(nbrPoints-1)
    // Abcisse utilisé pour la calcul de la diagonale de A 
    abscisse = linspace(-l+deltax,l-deltax,nbrPoints-2)
    // Abscisse2 utilisé pour la calcul de la diagonale inférieur --> au total, n - 1 points c'est à dire nbrPoints - 3 
    abscisse2 = linspace(-l+2*deltax ,l-deltax,nbrPoints-3)
    //A(i, i) = C_(i+0.5) + C(i-0.5)
    Adiag = exp(-(abscisse +deltax/2)/l) + exp(-(abscisse -deltax/2)/l)
    //A(i, i-1) = -C(i-0.5)
    Ainf = -exp(-(abscisse2 -deltax/2)/l)
    //B est une matrice colonne à n points = nbrPoints - 2 
    B=zeros(1,nbrPoints-2)
    B(1,1) = exp(-(abscisse(1) -deltax/2)/l)

endfunction

function main()
    // On observe l'erreur entre la solution exacte et la solution trouvé numériquement pour différentes valeurs de pas
    for i=1:4
        //utilisation de Cholesky
        [Adiag, Ainf, B] = genereMatricesAB(10^i, 5)
        [Ldiag, Linf] = factorise(Adiag,Ainf)
        Z = descente(Ldiag, Linf, B)
        ordonnee = remonte(Ldiag, Linf, Z)
        affiche(ordonnee, i)
    end
endfunction


function affiche(ordonnee, i)
    // fonction qui affiche sur un graphique deux courbes : la solution explicite et la solution numérique. Donne aussi la norme infini de l'écart entre les deux
        abscisse= linspace(-5,5,10^i-2)
        ordonneeReel = (1/(exp(-1)-exp(1)))*(exp(abscisse/5)-exp(1))
        normeInf = max(abs(ordonnee-ordonneeReel))
        disp("Erreur Max" + string(normeInf))
        subplot(2,2,i)
        // résultat numérique en bleu
        plot(abscisse, ordonnee)
        //résultat exacte en rouge
        plot(abscisse, ordonneeReel,'r')
        title("Approximation de u pour un pas de: "+string(2*5/(10^i)))
        xlabel("x")
        ylabel("température")
endfunction
