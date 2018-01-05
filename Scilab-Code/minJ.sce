// Script Q14, permet de trouver le minimum de la fonction J par la méthode de Gauss-Newton

// les deux fonctions suivantes sont utilisés pour effectuer des tests

exec("Q11flux.sce")

function fx=f(x)
    fx = 2*x.^2 + x
endfunction

function FX=F(x)
    FX= [f(2*x),f(3*x)]
endfunction

function norme = normeEuclidienneCarre(vecteur)
    //calcul de la norme euclidienne au carré
    taille = size(vecteur, "c")
    norme = 0
    for i = 1 : taille
        norme = norme +vecteur(i)^2
    end
endfunction

function Jprimex=calculJprime(x, fonction, Fcible)
    // calcul de J', utile pour le critère d'arret
    Jprimex = (2*numderivative(fonction, x)*(fonction(x)-Fcible)')/normeEuclidienneCarre(Fcible)
endfunction

function xcourant=minNewtonJ(epsilon, mini, maxi, fonction, Fcible)
    //fonction de calcul du minimum de J
    // fonction est mis en paramètre afin de pouvoir tester MinNewtonJ sur d'autres exemples
    xcourant = (mini + maxi)/2
    Jprimex = calculJprime(xcourant, fonction, Fcible)
    while abs(Jprimex) > epsilon
        deltaxk = -(numderivative(fonction, xcourant)*((fonction(xcourant))-Fcible)')/(numderivative(fonction, xcourant)*numderivative(fonction,xcourant)')
        xcourant = xcourant + deltaxk
        Jprimex = calculJprime(xcourant, fonction, Fcible)
    end
endfunction


Fcible = [-0.1,-0.18]
res = minNewtonJ(10^-5, -6,3, flux, Fcible)
res
