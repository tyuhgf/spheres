##################################################################
##################################################################
##################################################################
####                                                          ####
####  GAP-Program                                             ####
####                                                          ####
####                 `BISTELLAR'                              ####
####                                                          ####
####  Version Nov/2003                                        ####
####  by Frank H. Lutz, TU Berlin, Germany                    ####
####                                                          ####
####  available at                                            ####
####                                                          ####
####     http://www.math.TU-Berlin.de/diskregeom/stellar/     ####
####                                                          ####
####  (First Version Nov/1997                                 ####
####  by Anders Björner, KTH Stockholm, Sweden,               ####
####     Frank H. Lutz, TU Berlin, Germany)                   ####
####                                                          ####
##################################################################
##################################################################
##################################################################


##################################################################
##################################################################
##                                                              ##
## The GAP program BISTELLAR is a heuristics to                 ##
##                                                              ##
## 1. recognize combinatorial d-spheres; and thus can be used   ##
##    to test whether a pure simplicial complex is a            ##
##    combinatorial manifold;                                   ##
##                                                              ##
## 2. obtain small or even minimal triangulations of simplicial ##
##    manifolds by reducing the number of vertices of a given   ##
##    triangulation;                                            ##
##                                                              ##
## 3. determine the homeomorphism type of a simplicial manifold ##
##    (use version BISTELLAR_EQUIVALENT from                    ##
##                                                              ##
##        http://www.math.TU-Berlin.de/diskregeom/stellar/      ##
##                                                              ##
##    for manifolds different from a sphere).                   ##
##                                                              ##
## The program is based on a simulated-annealing type           ##
## random strategy to perform bistellar moves which             ##
## locally modify the triangulation but do not change           ##
## the PL homeomorphism type of a manifold.                     ##
## (See the first reference below for further information.)     ##
##                                                              ##
##################################################################
##################################################################


##################################################################
##################################################################
##                                                              ##
## How to use the program?                                      ##
##                                                              ##
## The program BISTELLAR can easily be adapted to deal with     ##
## individual test objects.                                     ##
##                                                              ##
## (A) Before running the program, you can either write         ##
##     the facets of a pure simplicial complex that you         ##
##     want to examine to a separate file                       ##
##                                                              ##
##        `BISTELLAR.testobject'                                ## 
##                                                              ##
##     or you can select a predefined complex:                  ##
##                                                              ##
##################################################################
##################################################################



##################################################################
##################################################################
##                                                              ##
## (B) Adjust the following parameters:                         ##
##                                                              ##
##################################################################
##################################################################

rounds_max:=50000; # <--- maximal numbers of rounds 
                   #
level:=1;          # <--- 1 = examine complex
                   #      2 = examine vertex links 
                   #
verbose:=9;        # <--- 0 = no comments
                   #      1 = type of move, current g-vector
                   #          and number of rounds 
                   #          are printed to the screen
                   #      2 = type of move, current f- and g-vector
                   #          and number of rounds 
                   #          are printed to the screen 
                   #          and to the file `BISTELLAR.log'

##################################################################
##################################################################
##                                                              ##
## (C) Fix random sequence (optional):                          ##
##                                                              ##
##################################################################
##################################################################

### fix random sequence ###

R_N:=1;   # <--- 1,2,3,4,... different but fixed random sequences 
          #                  are used for distinct integers,
          #                  allowing you to reproduce any 
          #                  `random result'.
          #
###########

R_X:=[66318732,86395905,22233618,21989103,237245480,264566285,240037038,
      264902875,9274660,180361945,94688010,24032135,106293216,27264613,
      126456102,243761907,80312412,2522186,59575208,70682510,228947516,
      173992210,175178224,250788150,73030390,210575942,128491926,194508966,
      201311350,63569414,185485910,62786150,213986102,88913350,94904086,
      252860454,247700982,233113990,75685846,196780518,74570934,7958751,
      130274620,247708693,183364378,82600777,28385464,184547675,20423483,
      75041763,235736203,54265107,49075195,100648387,114539755]; 

##################################################################
##################################################################
##                                                              ##
## (D) Bistellar flip strategy (optional):                      ##
##                                                              ##
##        There is, a priory, no optimal flip strategy for all  ##
##        types of test objects. Thus, for every individual     ##
##        test object you might have to adapt the flip strategy ##
##        and, in particular, the values of the parameters      ##
##        `heating' and `relaxation' in the `Main Part'         ##
##        of the program. For dimensions 1-6, we have           ##
##        included strategies (along with possible              ##
##        alterations) that work well in many cases.            ##
##        For dimensions d >= 7 we provide a simple             ##
##        general strategy.                                     ##
##                                                              ##
## How to call the program?                                     ##
##                                                              ##
## (E) Install GAP on your computer.                            ##
##                                                              ##
## (F) Adapt the program BISTELLAR as specified in (A) - (D),   ##
##     save it as a file, start GAP, and call BISTELLAR         ##
##     by the command                                           ##
##                                                              ##
##        gap> Read("BISTELLAR");                               ##
##                                                              ##
## What is the output of the program?                           ##
##                                                              ##
## (G) The program stops as soon as the total number            ##
##     of rounds is executed or in the case that no further     ##
##     options for flips are available (no_options_flag=1).     ##
##                                                              ##
##     level=1: Every time BISTELLAR finds a triangulation      ##
##              smaller than all previous ones, this complex    ##
##              together with its f- and g-vector and the       ##
##              current number of rounds is printed to          ##
##              the file `BISTELLAR.out'.                       ##
##                                                              ##
##              If the final complex is the boundary            ##
##              of a simplex, then the test object              ##
##              is a PL-sphere. In this case (flag_sphere=1)    ##
##              the final g-vector is zero.                     ##
##                                                              ##
##     level=2: If all vertex links are bistellarly equivalent  ##
##              to the boundary of a simplex, then the          ##
##              examined complex is a combinatorial manifold    ##
##              (flag_manifold=1).                              ##
##                                                              ##
## (H) The main (global) variables that contain information     ##
##     about the present complex in every round are:            ##
##                                                              ##
##        faces[max] - maximal faces;                           ##
##        faces[i]   - (i-1)-dimensional faces;                 ##
##        f          - f-vector                                 ##
##        g          - g-vector.                                ##
##                                                              ##
## (I) gap> quit;                                               ##
##                                                              ##
##################################################################
##################################################################

##################################################################
##################################################################
##                                                              ##
## Exercises to get acquainted with the program:                ##
##                                                              ##
## (a) Use the program to obtain the unique minimal             ##
##     7-vertex-triangulation of the 2-dimensional torus.       ##
##                                                              ##
## (b) Show that the boundaries of the 4- and 5-dimensional     ##
##     cross-polytope are bistellarly equivalent to the         ##
##     boundary of a simplex.                                   ##
##                                                              ##
## (c) Repeat (b) for boundaries of cyclic polytopes.           ##
##                                                              ##
## (d) k-fold suspensions of the Poincaré homology 3-sphere     ##
##     give examples of non-PL (k+3)-spheres for k >= 2.        ##
##     In particular, the 2-fold suspension of the provided     ##
##     16-vertex triangulation of the Poincaré sphere yields    ##
##     a non-PL 5-sphere with 20 vertices.                      ##
##                                                              ##
##     Try to find a triangulation of this non-PL sphere with   ##
##     18 vertices.                                             ##
##                                                              ##
## (e) In fact, more economical one-point suspensions of the    ##
##     Poincaré homology 3-sphere yield non-PL d-spheres        ##
##     with d+13 vertices in dimensions d >= 5.                 ##
##                                                              ##
##     Open problem: Do there exist non-PL d-spheres with       ##
##                   less than d+13 vertices?                   ##
##                   (At least d+6 vertices are needed by       ##
##                   a result of Barnette and Gannon,           ##
##                   Disc. Math. 16, 291-298 (1976).)           ##
##                                                              ##
##################################################################
##################################################################


##################################################################
##################################################################
##                                                              ##
## For an exposition of the program BISTELLAR                   ##
## see the reference:                                           ##
##                                                              ##
##   A. Björner, F.H. Lutz.                                     ##
##   Simplicial manifolds, bistellar flips                      ##
##   and a 16-vertex triangulation of the                       ##
##   Poincaré homology 3-sphere.                                ##
##   Exp. Math. 9, 275-289 (2000).                              ##
##                                                              ##
## Further references:                                          ##
##                                                              ##
##   F.H. Lutz.                                                 ##
##   Triangulated Manifolds with Few Vertices and               ##
##   Vertex-Transitive Group Actions.                           ##
##   Dissertation. Shaker Verlag, Aachen, 1999.                 ##
##   (http://www.math.tu-berlin.de/~lutz/dissertation.ps.gz)    ##
##                                                              ##
##   The GAP Group.                                             ##
##   GAP -- Groups, Algorithms, and Programming,                ##
##   Version 4.3, 2002.                                         ##
##   (http://www.gap-system.org)                                ##
##                                                              ##
##################################################################
##################################################################






##################################################################
##################################################################
###              Global variables and functions                ###
##################################################################
##################################################################


# facets:=[];
dim:=0; max:=0; faces:=[]; maxvertex:=0;
f:=[]; h:=[]; g:=[];
chi:=0; sum:=0;
raw_options:=[];
options:=[];
ball_boundary_faces:=[];
randomelement:=[];



########################
###    Suspension    ###
########################

Suspension := function()

   local element, i, A, B, C;

   maxvertex:=0;
   for element in facets do
       for i in element do
           if i > maxvertex then
              maxvertex:=i;
           fi;
       od;
   od;

   A:=[];

   for element in facets do 
 
       B:=[maxvertex+1];
       UniteSet(B,element);
       AddSet(A,B);

       C:=[maxvertex+2];
       UniteSet(C,element);
       AddSet(A,C);

   od;

   facets:=[];
   UniteSet(facets,A);
       
end;



#########################################################################
###             Functions to compute the f,g,h-vectors                ###
#########################################################################


### Function  FacetsTofvector  computes the f-vector from the  facets ###
  
FacetsTofvector := function()

   local element, k;

   ## Computing the dimension dim=max-1 of the simplicial complex. ##

   max:=0;
   for element in facets do
       if Length(element) > max then
          max:=Length(element);
       fi;
   od; 

   dim:=max-1;

   if verbose = 1 then
      Print("dim = ",dim,"\n\n");
   elif verbose = 2 then
      Print("dim = ",dim,"\n\n");
      AppendTo(log_file,"dim = ",dim,"\n\n");
   elif verbose =9 then
      AppendTo(log_file,"[\n");
   fi;

   ## Computing all faces ##

   faces:=[];
   f:=[];

   chi:=0;
   sum:=0;

   for k in [1..max] do

       faces[k]:=[];
       for element in facets do
           if Length(element) >= k then
              UniteSet(faces[k],Combinations(element,k));
           fi;
       od;
       f[k]:=Length(faces[k]);

       if verbose = 1 then
          Print("f_",k-1," = ",f[k],"   ");
       elif verbose = 2 then
          Print("f_",k-1," = ",f[k],"   ");
          AppendTo(log_file,"f_",k-1," = ",f[k],"   ");
       fi;

       chi:=chi+((-1)^(k+1))*f[k];
       sum:=sum+f[k];

   od;

   if verbose = 1 then
      Print("\n\n","chi = ",chi,"    ","sum = ",sum,"\n");
   elif verbose = 2 then
      Print("\n\n","chi = ",chi,"    ","sum = ",sum,"\n");
      AppendTo(log_file,"\n\n","chi = ",chi,"    ","sum = ",sum,"\n");
   fi;

end;



### Function  fvector2hvector  converts the f-vector to the h-vector ###
###             ( f --> h )                                          ###
  
fvector2hvector := function()

   local k, i;

   h:=[];
  
   for k in [1..max] do
       h[k]:=(-1)^k*Binomial(max,k);
       for i in [1..max] do
           h[k]:=h[k]+(-1)^(k-i)*Binomial(max-i,max-k)*f[i];
       od;
   od;

end;



### Function  hvector2gvector  converts the h-vector to the g-vector ###
###             ( h --> g )                                          ###
  
hvector2gvector := function()

   local k;

   g:=[];

   if QuoInt(max,2) >= 1 then

      g[1]:=h[1]-1;

      if verbose = 1 then
         Print("\n","g_",1," = ",g[1],"   ");
      elif verbose = 2 then
         Print("\n","g_",1," = ",g[1],"   ");
         AppendTo(log_file,"\n","g_",1," = ",g[1],"   ");
      fi;
      
   fi;

   for k in [2..QuoInt(max,2)] do

       g[k]:=h[k]-h[k-1];

       if verbose = 1 then
          Print("g_",k," = ",g[k],"   ");
       elif verbose = 2 then
          Print("g_",k," = ",g[k],"   ");
          AppendTo(log_file,"g_",k," = ",g[k],"   ");
       fi;

   od;

   if verbose = 1 then
      Print("\n");
   elif verbose = 2 then
      Print("\n");
      AppendTo(log_file,"\n");
   fi;

end;




##################################################################
###       Initional options for moves (and reverse moves)      ###
###       ( Reverse_k_Move = (max-k-1)-Move )                  ###
##################################################################


InitionalMoveOptions := function(r)

    local testelement, count, linkface, maxface;
       
    raw_options[r+1]:=[];

    for testelement in faces[max-r] do
        linkface:=[];
        count:=0;
        for maxface in faces[max] do
            if IsSubset(maxface,testelement) = true  then
               count:=count+1;
               UniteSet(linkface,maxface);
            fi;
        od;
        if count = r + 1  then
           SubtractSet(linkface,testelement);
           AddSet(raw_options[r+1],[testelement,linkface]);
        fi;
    od;

end;



##################################################################
###       Include options for moves (and reverse moves)        ###
##################################################################


### Include_MoveOptions ###

Include_MoveOptions := function(r)

   local element;

   for element in raw_options[r+1] do
       if r = 0 then
          AddSet(options,element);
       elif not element[2] in faces[Length(element[2])] then
          AddSet(options,element);
       fi;
   od;

end;

### Include_ReverseMoveOptions ###

Include_ReverseMoveOptions := function(r)

   Include_MoveOptions(max-r-1);

end;



###################################################################
###                  Moves (and reverse moves)                  ###
###################################################################


### Ball_Boundary ###

Ball_Boundary := function()

   local element, j, count, linkface, maxface;

   for element in ball_boundary_faces do

       count:=0;
       linkface:=[];
       for maxface in faces[max] do
           if IsSubset(maxface,element) = true  then
              count:=count+1;
              UniteSet(linkface,maxface);
           fi;
       od;
       SubtractSet(linkface,element);

       j:=1;
       while j <= Length(raw_options[max-Length(element)+1]) do
             if element = raw_options[max-Length(element)+1][j][1] then
                RemoveSet(raw_options[max-Length(element)+1],
                          raw_options[max-Length(element)+1][j]);
                j:=Length(raw_options[max-Length(element)+1]);
             fi;
             j:=j+1;
       od;
       if count = max - Length(element) + 1 then
          AddSet(raw_options[max-Length(element)+1],[element,linkface]);
       fi;

   od;

end;



### 0_Move ###

0_Move := function()

   local element, i, A, linkface, maxface;

   maxvertex:=maxvertex+1;
   
   RemoveSet(faces[max],randomelement[1]);
   f[max]:=f[max]-1;
   RemoveSet(raw_options[1],randomelement);
 
   for i in [0..(max-1)] do
       for element in Combinations(randomelement[1],i) do

           A:=[maxvertex];
           UniteSet(A,element);  
           AddSet(faces[i+1],A);
           f[i+1]:=f[i+1]+1;

           linkface:=[];
           for maxface in Combinations(randomelement[1],max-1) do
               if IsSubset(maxface,element) = true  then
                  UniteSet(linkface,maxface);
               fi;
           od;
           SubtractSet(linkface,element);

           AddSet(raw_options[max-Length(A)+1],[A,linkface]);

       od;
   od;

   ball_boundary_faces:=[];
   for i in [1..(max-1)] do
       UniteSet(ball_boundary_faces,Combinations(randomelement[1],i));
   od;

   Ball_Boundary();

end;



### Move ###

Move := function(r)

   local element, i, j, A, count, maxface, linkface, 
         new_facets, ball_interior_faces;

   if r = 0 then

      0_Move();

   else

      for i in [0..r] do
          for element in Combinations(randomelement[2],i) do

              A:=[];
              UniteSet(A,randomelement[1]);
              UniteSet(A,element);
              RemoveSet(faces[max-r+i],A);
              f[max-r+i]:=f[max-r+i]-1;

              j:=1;
              while j <= Length(raw_options[max-Length(A)+1]) do
                    if A = raw_options[max-Length(A)+1][j][1] then
                       RemoveSet(raw_options[max-Length(A)+1],
                                 raw_options[max-Length(A)+1][j]);
                       j:=Length(raw_options[max-Length(A)+1]);
                    fi;
                    j:=j+1;
              od;
    
          od;
      od;    

      new_facets:=[];
      ball_interior_faces:=[];
   
      for i in [0..(max-r-1)] do
          for element in Combinations(randomelement[1],i) do

              A:=[];
              UniteSet(A,randomelement[2]);
              UniteSet(A,element);
              AddSet(faces[r+i+1],A);
              f[r+i+1]:=f[r+i+1]+1;

              if i = max - r - 1 then
                 AddSet(new_facets,A);
              else
                 AddSet(ball_interior_faces,A);
              fi;

          od;
      od;
   
      for i in [0..(max-r-1)] do
          for element in Combinations(randomelement[1],i) do

              A:=[];
              UniteSet(A,randomelement[2]);
              UniteSet(A,element);

              linkface:=[];
              for maxface in new_facets do
                  if IsSubset(maxface,A) = true  then
                     UniteSet(linkface,maxface);
                  fi;
              od;
              SubtractSet(linkface,A);  

              AddSet(raw_options[max-Length(A)+1],[A,linkface]);

          od;
      od;

      ball_boundary_faces:=[];
      for element in new_facets do
          for i in [1..(max-1)] do
              UniteSet(ball_boundary_faces,Combinations(element,i));
          od;
      od;
      SubtractSet(ball_boundary_faces,ball_interior_faces);

      Ball_Boundary();

   fi;

end;



####################################################################
###                    Available complexes                       ###
####################################################################


###  BdCyclicPolytope(d,n)  computes the boundary complex ###
###  of the d-dimensional cyclic polytope on  n  vertices ###

BdCyclicPolytope := function(d,n)

   local element, i, A, B; 

   if verbose = 1 then
      Print("\n\n","BdCyclicPolytope(",d,",",n,"):\n");
   elif verbose = 2 then
      Print("\n\n","BdCyclicPolytope(",d,",",n,"):\n");
      PrintTo(log_file,"\n\n","BdCyclicPolytope(",d,",",n,"):\n");
   fi;

   facets:=[];

   if IsInt(d/2) = true then
      for element in Combinations([1..(n-d/2)],d/2) do
          A:=[];
          for i in element do 
              AddSet(A,i+Position(element,i)-1);    
              AddSet(A,i+Position(element,i));
          od;
          AddSet(facets,A);
      od;
      for element in Combinations([1..(n-2-(d-2)/2)],(d-2)/2) do
          A:=[];
          for i in element do 
              AddSet(A,i+Position(element,i));    
              AddSet(A,i+Position(element,i)+1);
          od;
          AddSet(A,1);
          AddSet(A,n);
          AddSet(facets,A);
      od;
   else
      for element in Combinations([1..(n-1-(d-1)/2)],(d-1)/2) do
          A:=[];
          B:=[];
          for i in element do 
              AddSet(A,i+Position(element,i));    
              AddSet(A,i+Position(element,i)+1);
              AddSet(B,i+Position(element,i)-1);    
              AddSet(B,i+Position(element,i));
          od;
          AddSet(A,1);
          AddSet(B,n);
          AddSet(facets,A);
          AddSet(facets,B);
      od;
   fi;

end;



###  Function  BdCrossPolytope(d)  computes the boundary complex ###
###  of the d-dimensional crosspolytope                          ###

BdCrossPolytope := function(d)

   local element, a, i, B, C; 

   if verbose = 1 then 
      Print("\n\n","BdCrossPolytope(",d,"):\n");
   elif verbose = 2 then 
      Print("\n\n","BdCrossPolytope(",d,"):\n");
      PrintTo(log_file,"\n\n","BdCrossPolytope(",d,"):\n");
   fi;

   facets:=[];

   for a in [0..(2^d-1)] do
       B:=[];
       C:=[];
       for i in [0..(d-1)] do
           if QuoInt(a,2^(d-i-1)) = 1 then    
              a:=a-2^(d-i-1);
              AddSet(B,(d-i)+d+1);
              Add(C,1);
           else 
              AddSet(B,d-i+1);
              Add(C,0);
           fi;
       od;
       AddSet(facets,B);
   od;

end;



### Suspended_Poincare_Homology_3_Sphere ###

Suspended_Poincare_Homology_3_Sphere := function(s)

   local i;

   if verbose = 1 then
      Print("\n\n","Suspended_Poincare_Homology_3_Sphere(",s,"):\n");
   elif verbose = 2 then
      Print("\n\n","Suspended_Poincare_Homology_3_Sphere(",s,"):\n");
      PrintTo(log_file,"\n\n","Suspended_Poincare_Homology_3_Sphere(",s,"):\n");
   fi; 

   facets:=[[1,2,4,9],[1,2,4,15],[1,2,6,14],[1,2,6,15],[1,2,9,14],
            [1,3,4,12],[1,3,4,15],[1,3,7,10],[1,3,7,12],[1,3,10,15],
            [1,4,9,12],[1,5,6,13],[1,5,6,14],[1,5,8,11],[1,5,8,13],
            [1,5,11,14],[1,6,13,15],[1,7,8,10],[1,7,8,11],[1,7,11,12],
            [1,8,10,13],[1,9,11,12],[1,9,11,14],[1,10,13,15],[2,3,5,10],
            [2,3,5,11],[2,3,7,10],[2,3,7,13],[2,3,11,13],[2,4,9,13],
            [2,4,11,13],[2,4,11,15],[2,5,8,11],[2,5,8,12],[2,5,10,12],
            [2,6,10,12],[2,6,10,14],[2,6,12,15],[2,7,9,13],[2,7,9,14],
            [2,7,10,14],[2,8,11,15],[2,8,12,15],[3,4,5,14],[3,4,5,15],
            [3,4,12,14],[3,5,10,15],[3,5,11,14],[3,7,12,13],[3,11,13,14],
            [3,12,13,14],[4,5,6,7],[4,5,6,14],[4,5,7,15],[4,6,7,11],
            [4,6,10,11],[4,6,10,14],[4,7,11,15],[4,8,9,12],[4,8,9,13],
            [4,8,10,13],[4,8,10,14],[4,8,12,14],[4,10,11,13],[5,6,7,13],
            [5,7,9,13],[5,7,9,15],[5,8,9,12],[5,8,9,13],[5,9,10,12],
            [5,9,10,15],[6,7,11,12],[6,7,12,13],[6,10,11,12],[6,12,13,15],
            [7,8,10,14],[7,8,11,15],[7,8,14,15],[7,9,14,15],[8,12,14,15],
            [9,10,11,12],[9,10,11,16],[9,10,15,16],[9,11,14,16],[9,14,15,16],
            [10,11,13,16],[10,13,15,16],[11,13,14,16],[12,13,14,15],[13,14,15,16]];

   for i in [1..s] do
       Suspension();
   od;

end;



### 9_Vertex_Torus ###

9_Vertex_Torus := function()

   local element;

   if verbose = 1 then
      Print("\n\n","9_Vertex_Torus:\n");
   elif verbose = 2 then
      Print("\n\n","9_Vertex_Torus:\n");
      PrintTo(log_file,"\n\n","9_Vertex_Torus:\n");
   fi; 

   facets:=[[1,2,5],[1,2,6],[1,3,4],[1,3,9],[1,4,6],[1,5,9],
            [2,3,7],[2,3,8],[2,5,8],[2,6,7],[3,4,7],[3,8,9],
            [4,5,7],[4,5,8],[4,6,8],[5,7,9],[6,7,9],[6,8,9]];

end;



### TestObject ### 

TestObject := function()

   if verbose = 1 then
      Print("\n\n","TestObject:\n");
   elif verbose = 2 then 
      Print("\n\n","TestObject:\n");
      PrintTo(log_file,"\n\n","TestObject:\n");
   fi; 

   # Read(file);

end;





###################################################################
###################################################################
###                       Main Part                             ###
###################################################################
###################################################################



if type_of_object = 1 then
   TestObject();
elif type_of_object = 2 then
   BdCyclicPolytope(polytope_dimension,number_of_vertices);
elif type_of_object = 3 then
   BdCrossPolytope(polytope_dimension);
elif type_of_object = 4 then
   Suspended_Poincare_Homology_3_Sphere(number_of_suspensions);
else
   9_Vertex_Torus();
fi;

if level = 1 then
   flag_sphere:=1;
   n_objects:=1;
   if verbose = 1 then
      Print("\n\n");
   elif verbose = 2 then
      Print("\n\n");
      AppendTo(log_file,"\n\n","facets:=",facets,";\n\n\n");
   fi;
elif level = 2 then
   flag_manifold:=1;
   top_facets:=ShallowCopy(facets);
   nodes:=[];
   for element in top_facets do
       UniteSet(nodes,element);
   od;
   n_objects:=Length(nodes);
fi;

for counter in [1..n_objects] do

    if level = 2 then
       facets:=[];
       for element in top_facets do
           if nodes[counter] in element then
              copy_element:=ShallowCopy(element);
              RemoveSet(copy_element,nodes[counter]);
              AddSet(facets,copy_element);
           fi;
       od;
    fi;

    #############################################
    ### computing initional options for moves ###
    #############################################

    if level = 2 then
       if verbose = 1 then
          Print("\n\n","Examining vertex link ",nodes[counter],":\n\n");
       elif verbose = 2 then
          Print("\n\n","Examining vertex link ",nodes[counter],":\n\n");
          AppendTo(log_file,"\n\n","Examining vertex link ",nodes[counter],":\n\n");
       fi;
    fi;

    FacetsTofvector();
    fvector2hvector();
    hvector2gvector();

    if level = 1 then
       if verbose = 1 then
          Print("\n\n");
       elif verbose = 2 then
          Print("\n\n");
          AppendTo(log_file,"\n\n");
       fi;
    fi;

    maxvertex:=0;
    for element in faces[1] do
        if element[1] > maxvertex then
           maxvertex:=element[1];
        fi;
    od;

    for i in [0..(max-1)] do
        InitionalMoveOptions(i);
    od;

    ### initional parameters ###

    rounds:=1;
    heating:=0;
    relaxation:=0;
    no_options_flag:=0;
    flag_small:=0;
	
    g_min:=[];
    for i in g do
        Add(g_min,i);
    od;

    if level = 1 then
       PrintTo(out_file,
               # "#####\n#\n",
               # "# Triangulation produced by the GAP-program BISTELLAR:\n#\n",
               # "# f = ",f,"\n",
               # "# g = ",g,"\n#\n",
               # "# rounds = ",rounds,"\n",
               # "#\n#####\n\n",
               # "facets:=",faces[max],";\n");
               faces[max], "\n");
    fi;

    ######################################################
    ### Selection of strategies and execution of moves ###
    ######################################################

    while no_options_flag = 0 and rounds <= rounds_max do  

          ### strategy for selecting options ###   

          # <--- adapt bistellar flip strategies in dimensions 1-7
          #      or add a strategy for higher dimensions

          options:=[];

          if dim = 1 then

             Include_ReverseMoveOptions(0);
             if options = [] then
                no_options_flag:=1;
             fi;

          elif dim = 2 then
 
             Include_ReverseMoveOptions(0);
             if options = [] then
                Include_MoveOptions(1);
                if options = [] then
                   no_options_flag:=1;
                fi;
             fi;

          elif dim = 3 then

          #if rounds < 7 then
          #   Include_MoveOptions(0);
          #elif rounds < 90 then
          #   Include_MoveOptions(1);
          #elif IsInt(heating/10) = true and rounds <= 50000 then
          #   Include_MoveOptions(0);
          #   heating:=heating-1;
          #else
             if heating > 0 then
                Include_MoveOptions(1);
                if options = [] then
                   Include_ReverseMoveOptions(1);
                   heating:=0;
                fi;
                heating:=heating-1;
             else
                Include_ReverseMoveOptions(0);
                if options = [] then
                   Include_ReverseMoveOptions(1);
                   if options = [] then
                      Include_MoveOptions(1);
                      if relaxation = 10 then
                         heating:=15;
                         relaxation:=0;
                      fi;
                      relaxation:=relaxation+1;
                      if options = [] then
                         no_options_flag:=1;
                      fi;
                   fi;
                fi;
             fi;
          #fi;

          elif dim = 4 then

             if heating > 0 then
                #if IsInt(heating/20) = true and rounds <= 50000 then
                #   Include_MoveOptions(0);
                #else
                   Include_MoveOptions(1);
                   Include_MoveOptions(2);
                   if options = [] then
                      Include_ReverseMoveOptions(1);
                   fi;
                #fi;
                heating:=heating-1;
             else
                Include_ReverseMoveOptions(0);
                if options = [] then
                   Include_ReverseMoveOptions(1);
                   if options = [] then
                      Include_MoveOptions(1);
                      Include_MoveOptions(2);
                      if relaxation = 15 then 
                         heating:=20;
                         relaxation:=0;
                      fi;
                      relaxation:=relaxation+1;
                      if options = [] then
                         no_options_flag:=1;
                      fi;
                   fi;
                fi;
             fi;

          elif dim = 5 then

             #if rounds < 3 then
             #   Include_MoveOptions(0);
             #elif rounds < 20 then
             #   Include_MoveOptions(2);
             #   if options = [] then
             #      Include_MoveOptions(1);
             #   fi;
             #else
                if heating > 0 then
                   #if IsInt(heating/10) = true and rounds <= 50000 then
                   #   Include_MoveOptions(0);
                   #else
                      Include_MoveOptions(2);
                      if options = [] then
                         Include_MoveOptions(1);
                         if options = [] then
                            Include_ReverseMoveOptions(1);
                            Include_ReverseMoveOptions(2);
                            heating:=1;
                         fi;
                      fi;
                   #fi;
                   heating:=heating-1;
                else
                   Include_ReverseMoveOptions(0);
                   if options = [] then
                      Include_ReverseMoveOptions(1);
                      if options = [] then
                         Include_ReverseMoveOptions(2);
                         if options = [] then
                            Include_MoveOptions(1);
                            Include_MoveOptions(2);
                            if relaxation = 10 then 
                               heating:=20;
                               relaxation:=0;
                            fi;
                            relaxation:=relaxation+1;
                            if options = [] then
                               no_options_flag:=1;
                            fi;
                         fi;
                      fi;
                   fi;
                fi;
             #fi;

          elif dim = 6 then

             if heating > 0 then
                if IsInt(heating/2) then
                   Include_MoveOptions(2);
                   if options = [] then
                      Include_MoveOptions(1);
                      Include_MoveOptions(3);
                   fi;
                else
                   Include_MoveOptions(1);
                   Include_MoveOptions(2);
                   Include_MoveOptions(3); 
                fi;
                if options = [] then
                   Include_ReverseMoveOptions(1);
                   Include_ReverseMoveOptions(2);
                   heating:=0;
                   if options = [] then
                      no_options_flag:=1;
                   fi;
                fi;
                heating:=heating-1;
             else
                Include_ReverseMoveOptions(0);
                if options = [] then
                   Include_ReverseMoveOptions(1);
                   if options = [] then
                      Include_ReverseMoveOptions(2);
                      if options = [] then
                         Include_MoveOptions(1);
                         Include_MoveOptions(2);
                         Include_MoveOptions(3);
                         #if g[1] > 11 then
                         #   if relaxation = 35 then
                         #       heating:=100;
                         #       relaxation:=0;
                         #   fi;
                         #elif g[1] > 9 then
                         #   if relaxation = 30 then
                         #      heating:=75;
                         #      relaxation:=0;
                         #   fi;
                         #elif g[1] > 7 then
                         #   if relaxation = 25 then
                         #      heating:=50;
                         #      relaxation:=0;
                         #   fi;
                         #else
                            if relaxation = 20 then
                               heating:=25;
                               relaxation:=0;
                            fi;
                         #fi;
                         relaxation:=relaxation+1;
                         if options = [] then
                            no_options_flag:=1;
                         fi;
                      fi;
                   fi;
                fi;
             fi;

          else

             if heating > 0 then
                t:=Int((dim+1)/2)-1;
                while options = [] and t > 0 do
                      Include_MoveOptions(t);
                      t:=t-1;
                od;
                if options = [] then
                   for t in [1..Int(dim/2)] do
                       Include_ReverseMoveOptions(t);
                   od;
                   heating:=1;
                fi;
                heating:=heating-1;
             else
                t:=0;
                while options = [] and t <= Int((dim+1)/2)-1 do
                      Include_ReverseMoveOptions(t);
                      t:=t+1;
                od;
                if options = [] then
                   for t in [1..Int(dim/2)] do
                       Include_MoveOptions(t);
                   od;
                   if relaxation = 10 then 
                      heating:=20;
                      relaxation:=0;
                   fi;
                   relaxation:=relaxation+1;
                   if options = [] then
                      no_options_flag:=1;
                   fi;
                fi;
             fi;

          fi;

          ###############################
          ### choosing move at random ###
          ###############################

          if no_options_flag = 0 then
    
             randomelement:=RandomList(options);
        
             if Length(randomelement[1]) <= Int(max/2) then
                Move(max-Length(randomelement[1]));
                if verbose = 1 then
                   Print("\n","Reverse-",Length(randomelement[1])-1,"-Move   (",rounds," rounds)");
                elif verbose = 2 then
                   Print("\n","Reverse-",Length(randomelement[1])-1,"-Move   (",rounds," rounds)");
                   AppendTo(log_file,"\n","Reverse-",Length(randomelement[1])-1,"-Move   (",rounds," rounds)");
                   AppendTo(log_file,"\n","randomelement=",randomelement);
                elif verbose = 9 then
                   AppendTo(log_file, randomelement, ",\n");
                fi;
             else
                Move(max-Length(randomelement[1]));
                if verbose = 1 then
                   Print("\n",max-Length(randomelement[1]),"-Move   (",rounds," rounds)");
                elif verbose = 2 then
                   Print("\n",max-Length(randomelement[1]),"-Move   (",rounds," rounds)");
                   AppendTo(log_file,"\n",max-Length(randomelement[1]),"-Move   (",rounds," rounds)");
                   AppendTo(log_file,"\n","randomelement=",randomelement);
                elif verbose = 9 then
                   AppendTo(log_file, randomelement, ",\n");
                fi;
             fi;
       
             if verbose = 2 then
                Print("\n");
                AppendTo(log_file,"\n");
                for k in [1..max] do
                    Print("f_",k-1," = ",f[k],"   ");
                    AppendTo(log_file,"f_",k-1," = ",f[k],"   ");
                od;;
             fi;

             fvector2hvector();
             hvector2gvector();

             i:=0;
             while i < Int((dim+1)/2) do
                   i:=i+1;
                   if g[i] < g_min[i] then 
                      flag_small:=1;
                      i:=Int((dim+1)/2);
                   elif g[i] > g_min[i] then 
                      i:=Int((dim+1)/2);
                   fi;
             od;
             if flag_small = 1 then
                g_min:=[];
                for gg in g do
                    Add(g_min,gg);
                od;
                vertices:=[];
                for facet in faces[max] do
                    UniteSet(vertices,facet);
                od;
                relabeled_complex:=[];
                for facet in faces[max] do
                    relabeled_facet:=[];
                    for t in facet do
                        AddSet(relabeled_facet,Position(vertices,t));
                    od;
                    AddSet(relabeled_complex,relabeled_facet);
                od;
                if level = 1 then
                   PrintTo(out_file,
                           # "#####\n#\n",
                           # "# Triangulation produced by the GAP-program BISTELLAR:\n#\n",
                           # "# f = ",f,"\n",
                           # "# g = ",g,"\n#\n",
                           # "# rounds = ",rounds,"\n",
                           # "#\n#####\n\n",
                           # "facets:=",relabeled_complex,";\n");
                           relabeled_complex, "\n");
                fi;
                flag_small:=0;
             fi;

             rounds:=rounds+1;
    
          fi;

    od;

    if not Length(faces[max]) = dim+2 then
       if level = 1 then
          flag_sphere:=0;
       elif level = 2 then
          flag_manifold:=0; 
       fi;
    fi;

od;

if level = 1 and flag_sphere = 1 then
   if verbose = 1 then
      Print("\n\n","The examined complex is a sphere!!!\n\n\n");
   elif verbose = 2 then
      Print("\n\n","The examined complex is a sphere!!!\n\n\n");
      AppendTo(log_file,"\n\n","The test object is a sphere!!!\n\n\n");
   elif verbose = 9 then
      Print("The examined complex is a sphere!!!\n");
      AppendTo(log_file,"[]]\n");
   fi;
elif level = 2 and flag_manifold = 1 then
   if verbose = 1 then
      Print("\n\n","The examined complex is a combinatorial manifold!!!\n\n\n");
   elif verbose = 2 then
      Print("\n\n","The examined complex is a combinatorial manifold!!!\n\n\n");
      AppendTo(log_file,"\n\n","The test object is a combinatorial manifold!!!\n\n\n");
   elif verbose = 9 then
      Print("The test object is a combinatorial manifold!!!\n");
      AppendTo(log_file,"[]]\n");
   fi;
fi;

