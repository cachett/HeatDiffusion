// Scripts Q10, résoud numériquement l'équation (2) et montre la convergence vers la solution du problème stationnaire.

// On fixe u0(t)=1 qqsoit t, u^(0)(x)=0 qqsoit x dans ]-l,l[
// téta = 1/2 Cranck-Nicholson

clc;

exec("Q7pbstationnaire.sce")

//initialisation des variables
l = 10
nbrPoints = 500
deltax = 2*l/(nbrPoints-1)
deltat = 0.02
mu = deltat/deltax^2


function MembreDroit = calculMembreDroit(U, Ndiag, Ninf, B)
    // NU^(k) + B
        MembreDroit = Ndiag.* U  + [0,Ninf.*U(1:size(U,"c")-1)] + [Ninf .* U(2:size(U,"c")),0] + mu*B
endfunction


function affiche(i,abscisse, ordonnee)
    //fonction d'affichage pour les solutions numériques et stationnaires
        ordonneeReel = (1/(exp(-1)-exp(1)))*(exp(abscisse/l)-exp(1))
        subplot(2,2,i)
        //affichage de la solution numérique en bleu
        plot(abscisse, ordonnee)
        //affichage de la solution stationnaire en rouge
        plot(abscisse, ordonneeReel, 'r')
        title("Approximation de U au temps t = " + string((10^i)*deltat))
        xlabel("x")
        ylabel("température")
endfunction



function U=iter(i, Ndiag, Ninf, B, MdiagFact, MinfFact)
    // fonction qui retourne U^(k)
    U = zeros(1,nbrPoints-2)
    U(1,1) = 1
    for k=1:10^i
        // MembreDroite = NU^(k) + B 
        MembreDroit = calculMembreDroit(U, Ndiag, Ninf, B)
        Z = descente(MdiagFact, MinfFact, MembreDroit)
        // U^(k+1) = M^-1 N U^(k) + M^-1 B
        U = remonte(MdiagFact, MinfFact, Z)
    end
endfunction


function main()
    // On initialise les matrices M, N et B
    [Adiag, Ainf, B] = genereMatricesAB(nbrPoints, l)//NbrPoints prend en compte les bornes!! Donc n = nbrPoints-2
    Mdiag = 1 + (1/2)*mu*Adiag
    Minf = (1/2)*mu*Ainf
    Ndiag = 1 - (1/2)*mu*Adiag
    Ninf = -(1/2)*mu*Ainf
    [MdiagFact,MinfFact] = factorise(Mdiag, Minf)
    abscisse = linspace(-l,l,nbrPoints-2)
    // On calcul itérativement les U^(k) pour différent tk
    for i=1:4
        U=iter(i, Ndiag, Ninf, B, MdiagFact, MinfFact)
        affiche(i, abscisse, U)
    end
endfunction






