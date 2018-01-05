// Scripts Q11, renvoie le couple F(xd)=(f(tinter),f(tfin)) en fonction de la position xd du défaut
clc;

exec("Q3factorise.sce")
exec("Q4descente.sce")
exec("Q5remonte.sce")


//initialisation des variables
l = 10
a = 0.8
nbrPoints = 20
nbrTemps = 3000
T=60
deltax = 2*l/(nbrPoints-1)
deltat = T/nbrTemps
mu = deltat/deltax^2



function Cx = C(xd, x)
    //calcul la conductivité donnée par le modèle (1) suivant xd, au point x.
    Cx = 1 - a*exp((-(x-xd).^2 )/4)
endfunction

function [Mdiag, Minf, Ndiag, Ninf]=genereMatricesMN(xd)
    // génération de M = I + teta * mu * A et N = I - teta * mu * A
    // abcisse utilisée pour le calcul de la diagonale
    abscisse = linspace(-l+ deltax,l - deltax, nbrPoints - 2)
    // abscisse2 utilisée pour le calcul de la sous diagonale
    abscisse2 = linspace(-l + 2*deltax, l - deltax, nbrPoints -3)
    Adiag = C(xd, abscisse + deltax/2) + C(xd, abscisse - deltax/2)
    Ainf = -C(xd, abscisse2  - deltax/2)
    Mdiag = 1 + (1/2)*mu*Adiag
    Minf = (1/2)*mu*Ainf
    Ndiag = 1 - (1/2)*mu*Adiag
    Ninf = -(1/2)*mu*Ainf
    
endfunction


function MembreDroit = calculMembreDroit(U, Ndiag, Ninf, B)
    // calcul de NU^(k) + B(k)
    MembreDroit = Ndiag.* U  + [0,Ninf.*U(1:size(U,"c")-1)] + [Ninf .* U(2:size(U,"c")),0] + mu*B
endfunction


function [Uinter, Ufin]=iterTemps(Mdiag, Minf, Ndiag, Ninf)
    //permet de trouver u(x1,2/3T) et u(x1,T)
    U = zeros(1,nbrPoints-2)
    [MdiagFact, MinfFact] = factorise(Mdiag, Minf)
    B=zeros(1,nbrPoints-2)
    for k = 1 : nbrTemps - 1 
        //Premier terme de B : C1/2(u0(tk+1)téta + u0(tk)(1-téta)), le reste vaut 0
        //avec u0(t)= (t/T)^2 
        B(1,1) = C(xd,-l+deltax/2)*((((k+1)*deltat)/T)^2+((k*deltat)/T)^2)/2
        MembreDroit = calculMembreDroit(U, Ndiag, Ninf, B)
        // calcul de U^(k+1) = M^-1*N*U^(k) + M^-1 * B
        Z = descente(MdiagFact, MinfFact, MembreDroit)
        U = remonte(MdiagFact, MinfFact, Z)
        //On relève la coordonnée de U qui nous intéresse au moment opportun
        if k == 2*nbrTemps/3 then
            abscisse=linspace(-l, l, nbrPoints-2)
            plot(abscisse, U, "r")
             title("Représentation fonction de Uinter et Ufin avec xd = -6")
             xlabel("x")
             ylabel("U")

            Uinter = U(1,1)
        end
    end
    abscisse=linspace(-l, l, nbrPoints-2)
    plot(abscisse, U)
    Ufin = U(1,1)
endfunction


function Fxd = flux(xd)
    //On génère les matrices nécessaires dépendantes de xd
    [Mdiag, Minf, Ndiag, Ninf] = genereMatricesMN(xd)
    // On retrouve les valeurs u(x1, 2/3T) et u(x1,T) en itérant sur le vecteur U
    [Uinter, Ufin] = iterTemps(Mdiag, Minf, Ndiag, Ninf) 
//    disp("Uinter"+ string(Uinter))
    //On retourne les valeurs du flux sous forme d'un couple
    Finter = (C(xd,-l+deltax/2)*(Uinter-(4/9))/deltax) - ((2*deltax)/(3*T)) 
    Ffin = (C(xd,-l+(deltax/2))*(Ufin-1)/deltax) - (deltax/T)
    Fxd = [Finter, Ffin]
   disp(Fxd)
endfunction

flux(-6)

//function test()
    // fonction de test pour vérifier l'évolution du flux selon l'endroit de déformation
  //  for i=-20:20
    //    flux(i/2)
    //end
//endfunction


