Nsize = 5000;
"Init the matrix"
Timing[A = Table[2*KroneckerDelta[i,j] - KroneckerDelta[i-1,j] - KroneckerDelta[i+1,j],{i,1,Nsize},{j,1,Nsize}];]
"Init the vector"
Timing[b = Table[-1/((Nsize+2)*(Nsize+1)),{i,1,Nsize}];]
"Solve the linear system"
Timing[x = LinearSolve[A,b];]
