//Script Q13, dichotomie afin de trouver le min de la fonction de coût J

exec("Q12.sci")

 //test de dichotomie
function testx=test(x)
    testx = (x-2).^2
endfunction
//retourne bien x = 2 


function minJ=dichotomie(fonction, epsilon, mini, maxi)
    //renvoie une approximation du minimum d'une fonction
    //fonction est un paramètre afin de pouvoir effectuer des tests
    while abs(maxi-mini)>epsilon
        abscisse = linspace(mini, maxi, 5)
        for i=1:5
            ordonnee(i) = fonction(abscisse(i))
        end
        // si J(x1)> J(x2)
        if ordonnee(1)>ordonnee(2) then
            // si J(x2) > J(x3)
            if ordonnee(2)>ordonnee(3) then
                mini = abscisse(2)
            // si J(x2) =< J(x3)
            else
                mini = abscisse(1)
                maxi = abscisse(3)
            end
         // si J(x1) =< J(x2)
        else
            maxi = abscisse(2)
        end
    end
    minJ = (maxi + mini)/2
endfunction

res=dichotomie(J, 10^-5, -6, 3)
disp(res)

