// Script Q12, écris et représente graphiquement la fonction de cout J avec xd comme argument

exec("Q11flux.sce")

function norme = normeEuclidienneCarre(vecteur)
    //calcul de la norme euclidienne au carrée
    taille = size(vecteur, "c")
    norme = 0
    for i = 1:taille
        norme = norme +vecteur(i)^2
    end
endfunction



function Jxd=J(xd)
    // fonction de calcul de J
    Fcible = [-0.1,-0.18]
    // Fxd est une liste a deux élements [Finter, Ffin]
    Fxd = flux(xd)
    Jxd = normeEuclidienneCarre(Fxd-Fcible) / normeEuclidienneCarre(Fcible)
endfunction

function TraceJ()
    //fonction d'affichage du flux en fonction de xd 
        abscisse = linspace(-6,3,60)
    for i=1:60
        ordonnee(1,i) = J(abscisse(i))
    end
    plot(abscisse, ordonnee)
    title("Représentation fonction de coût J")
    xlabel("x")
    ylabel("J(x)")

endfunction

