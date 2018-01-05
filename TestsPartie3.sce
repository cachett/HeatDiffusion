// Script de test pour la partie 3. 

exec("Q3factorise.sce")
exec("Q4descente.sce")
exec("Q5remonte.sce")

//A = [2,-1,0,0;
//    -1,2,-1,0;
//     0,-1,2,-1;
//     0,0,-1,2]
//Adiag = [2,2,2,2]
//Ainf = [-1,-1,-1]
//B = [1;2;3;4]
//
//disp("On resoud le system AX=B")
//disp(A)
//disp("Et B=")
//disp(B)
//
//disp("Avec la méthode classique")
//X = inv(A)*B
//disp(X)
//
//disp("Avec Cholesky")
//[Ldiag, Linf] = factorise(Adiag,Ainf)
//Z = descente(Ldiag, Linf, B)
//X2 = remonte(Ldiag, Linf, Z)

// affiche le temps de réalisation de la méthode de Cholesky 
function test_cholesky()
    Adiag = 2*ones(1,2000)
    Ainf = -ones(1,1999)
    B = ones(2000, 1)
    tic()
    [Ldiag, Linf] = factorise(Adiag,Ainf)
    for k=1:300
        Z = descente(Ldiag, Linf, B)
        X2 = remonte(Ldiag, Linf, Z)
    end
    temps=toc()
    disp("fini en "+string(temps) + " secondes")

endfunction



