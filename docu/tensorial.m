(* :Context: tensorial` *)

(* :Title: Tensor formulation *)

(* :Author: Renán Cabrera  *)

(* :Summary: This is package for Tensor calculus


*)


(* :Package Version: 1.3       December 2000

This is one of a set of packages to be published in a book at the end of the year 2000.
This package was developed by the Author Lic. Renan Cabrera in collaboration with the 
"Academia Nacional de Ciencias de Bolivia" and can be used and distributed freely 
for non commercial purpose, always mentioning the author every time used.   

This software package and its accompanying documentation are provided without guarantee
of support or maintenace
*)

(*
   If you want to help me to improve this package or you have any suggestion or you want
   the last version please contact me at
    "renanbo@hotmail.com" or "rencabla@ceibo.entelnet.bo    "
*)

(*  Mathematica is a registered trademark of Wolfram Research, Inc.     *)

BeginPackage["tensorial`"]

Void = "  ";
Tensor::usage = "Tensor[A,sub,sup] represents a tensor with label A and indexes _sub and  ^sup"

SetSystem::usage = "SetSystem[Var,g,M] set up the geometrical space. Var is the label for the
coordinates, g is the label for the Metric Tensor and M is the matrix containing the 
particular values of the Metric Tensor"

ExpandIndex::usage = "ExpandIndex[   Tensor , ind_List ],  Expand a tensor over free indexes ind" 

SetTensorRule::usage = "SetTensorRule[Tensor,Values], Defines a set of rules for a expanded 
 Tensor with particular Values"

ExpandDoubleIndex::usgae = "ExpandDoubleIndex[ Tensor ,ind], Expand double indexes like a Sum over indexes ind"

Kronecker::usage = "Kronecker[i,j] ,is the ordinary Kroneckers Delta,
  Kronecker[ i,j... , k,l,m..... ] is the generalized Kroneckers Delta"

TensorSimplify::usage = "TensorSimplify[ Tensor ], simplify Tensor   "

ToCovariantIndex::usage = " ToCovariantIndex[w,ind,newind] Brings Tensor w to 
   Covariant Index ind with new label newind"

ToContravariantIndex::usage = " ToCovariantIndex[w,ind,newind] Brings Tensor w to 
   Contravariantvariant Index ind with new label newind"

LeviCivitaRule::usage = " LeviCivitaRule "
\[Epsilon] ::usage = " LeviCivita label"
\[Delta] ::usage = " Delta Kronecker label"

IndexDerivative::usage = " IndexDerivative[ T , Tensor[x,{Void},{i}] ] , Represents the
partial derivative of T respect the coordinate i "

CovariantDerivative::usage = " CovariantDerivative[ T , Tensor[x,{Void},{i}] ]  Represents
the Covariant Derivative respect the index i"

AbsoluteDerivative::usage = "AbsoluteDerivative[f,t,{x,g}], AbsoluteDerivative of f over x 
with g Metric Tensor and coordinates x"

IntrinsecDerivative::usage = " IntrinsecDerivative[ Tensor , parameter ,{x,g}] ,
IntrinsecDerivative with coordinates x and MetricTensor g  "

EvaluateIndexDerivative::usage = "EvaluateIndexDerivative[ T ], Evaluate IndexDerivatives in T"

Christoffel1::usage = " Cristoffel[g,x][{i,k},{a}] symbols of class 1 with metric Tensor g 
                        and coordinates x"
Christoffel2::usage = " Cristoffel2[g,x][{i,k},{a}] symbols of class 2 with metric Tensor g 
                        and coordinates x"

DummyVariableSimplify::usage = 
             "DummyVariableSimplify[ T  ], Simplify dummy variables in Tensor expresions T"

(*::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::::*)

Begin["`Private`"]

(*.............................................*)
(*...............Format Display................*)
(*.............................................*)

Format[ Tensor[A_] ]:= A;

Format[ Tensor[A_,sub_,sup_] ]:=
Subsuperscript[A, SequenceForm@@sub , SequenceForm@@sup ];

Unprotect[Power];
Format[ Power[A_,sup_] ]:=Superscript[ SequenceForm[" ",A  ," "]  ,sup];
Protect[Power];

(*.............................................*)
(*.................SetSystem...................*)
(*.............................................*)

SetSystem[var_,g_,mat_]:= Module[{i,j,w},
   NDim = Length[mat];
      
   gExpanded  = ExpandIndex[ Tensor[g,{i,j},{Void,Void}]  , {i,j} ]; 
   gInvExpanded = ExpandIndex[ Tensor[g,{Void,Void},{i,j}]  , {i,j} ]; 
    
   Evaluate[ gExpanded ]=    mat ;
   Evaluate[  gInvExpanded  ] = Inverse[mat] ;    


IDXX[ Tensor[var,{Void},{i_}] , Tensor[var,{Void},{j_}]  ] := 
      Kronecker[j,i];

IDXX[ IndexDerivative[F_,Tensor[var,{Void},{i_}] ]  ,       IndexDerivative[G_,Tensor[var,{Void},{j_}]]   ] :=  Kronecker[i,j];

TensorMetricRuleSimplify = { 
     w_ Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[w,i,j]/;MemberAllQ[w,i]  ,
     w_ Tensor[g,{i_,j_},{Void,Void}] :> ToCovariantIndex[w,j,i]/;MemberAllQ[w,j]  ,

     w_ Tensor[g,{Void,Void},{i_,j_}] :> ToContravariantIndex[w,i,j]/;MemberAllQ[w,i]  ,
     w_ Tensor[g,{Void,Void},{i_,j_}] :> ToContravariantIndex[w,j,i]/;MemberAllQ[w,j]  ,

        Tensor[g,{i_,Void},{Void,j_}] :> Kronecker[i,j],
        Tensor[g,{Void,i_},{j_,Void}] :> Kronecker[i,j]  } ;

CovariantDerivative[ Tensor[g,{i_,j_},{Void,Void}]  ,  __ ] = 0;
CovariantDerivative[ Tensor[g,{Void,Void},{i_,j_}] ,  __ ] = 0;

ToContravariantIndex[ w_  ,ind_ ,toind_] := 
   w Tensor[g,{Void,Void},{ind,toind}];
             ]			
(*.............................................*)
(*...............ExpandIndex...................*)
(*.............................................*)

ExpandIndex[ T_ , ind_] :=  
  Module[{ni,dum},
		
  dum = {#,Sequence@@{1,NDim}} &/@  ind ;
  Table[ T, Evaluate[Sequence@@dum] ]

]/;VectorQ[ind];

(*............................................*)	  
(*............ExpandDoubleIndex...............*)	
(*............................................*)

ExpandDoubleIndex[w_,ind_]:=ExpandDoubleIndexINT[ Expand[w] , ind ];

ExpandDoubleIndexINT[w_Plus,ind_]:=  Map[ExpandDoubleIndexINT[#,ind]& , w ];
ExpandDoubleIndexINT[w_,ind_List]:= Fold[ExpandDoubleIndexINT[#1,#2]& , w , ind]
		
ExpandDoubleIndexINT[w_,ind_] := w/;FreeAllQ[w,ind]

ExpandDoubleIndexINT[w_,ind_]:=Module[{q},
       Plus@@Table[ w /.ind->q, {q,1,NDim} ]  
				]; 

(*...........................................*)
(*...............TensorSimplify..............*)
(*...........................................*)

TensorSimplify[w_]:=Module[{Regla},
	Regla = {Sequence@@KroneckerRuleSimplify,
                 Sequence@@TensorMetricRuleSimplify
                 };
        w //.Regla
		]

(*............................................*)
(*...............Kronecker..LeviCivita........*)
(*............................................*)

Format[ Kronecker[x_,y_] ] =  
 Tensor[ \[Delta] ,x,y]; 

Kronecker[i_?NumericQ, i_?NumericQ] :=  1;
Kronecker[i_?NumericQ, j_?NumericQ] := 0;
Kronecker[i_,j_] :=  0/; NumericQ[i-j] ;
		
KroneckerRuleSimplify = 
{
T_ Kronecker[i_,j_]:> ( T/. i->j)/; TensorRank[Position[   T , i]]>=2
		,
T_ Kronecker[i_,j_]:> ( T/. j->i)/;TensorRank[Position[   T , j]]>=2
		};
(*..................*)

VectorNumericQ[w_List]:= And@@Map[ NumericQ , w ];

Tensor[ \[Epsilon] ,x_ ,{Void...}]  := 
      Signature[x]/; VectorNumericQ[x];

Tensor[ \[Epsilon] ,{Void...},x_]  := 
      Signature[x]/; VectorNumericQ[x];

Kronecker[x_List,y_List] := 
Tensor[ \[Epsilon],x, Table[Void,{NDim}] ]*
Tensor[\[Epsilon],Table[Void,{NDim}],y]/;VectorNumericQ[x]&&VectorNumericQ[y];

(*...........................................*)
(*...............SetTensorRule...............*)
(*...........................................*)

SetTensorRule[x_,y_]:= MapThread[Rule,{x,y}, TensorRank[y]  ]//Flatten	        

(*............................................*)
(*..........ToCovariantIndex..................*)
(*............................................*)

ToCovariantIndex[w_,ind_,toind_]:= w /;FreeAllQ[w,ind];

ToCovariantIndex[ Tensor[w_,sub_,sup_]  ,ind_,toind_ ]:=Tensor[w,sub,sup]/;FreeQ[sup,ind];

ToCovariantIndex[ Tensor[w_,sub_,sup_]  ,ind_,toind_ ]:=Module[{pos,SUB,SUP},
 pos =		Position[sup,ind];
 SUB = ReplacePart[sub,toind,pos];
 SUP = ReplacePart[sup,Void,pos];
		 Tensor[w,SUB,SUP]
		];

ToCovariantIndex[ T_Times  ,ind__ ]:=
  Map[ToCovariantIndex[ #  ,ind ]&, T ];

ToCovariantIndex[ w_Plus  ,ind__ ]:= 
    Map[ToCovariantIndex[ #  ,ind ]&,  Expand[w] ];

ToCovariantIndex[ w_  ,ind_List ,toind_List]:=
Fold[
ToCovariantIndex[ #1  ,Part[#2,1],Part[#2,2] ]&,w, Transpose[{ind ,toind}]];


(*...............*)

ToContravariantIndex[w_,ind_,toind_]:= w /;FreeAllQ[w,ind];

ToContravariantIndex[ Tensor[w_,sub_,sup_]  ,ind_,toind_ ]:=Module[{pos,SUB,SUP},
 pos =		Position[sub,ind];
 SUB = ReplacePart[sub,Void,pos];
 SUP = ReplacePart[sup,toind,pos];
		 Tensor[w,SUB,SUP]
		]/; MemberQ[sub,ind];

ToContravariantIndex[ T_Times  ,ind__ ]:=
  Map[ToContravariantIndex[ #  ,ind ]&, T ];

ToContravariantIndex[ w_Plus  ,ind__ ]:= 
    Map[ToContravariantIndex[ #  ,ind ]&,  Expand[w] ];

ToContravariantIndex[ w_  ,ind_List ,toind_List]:=
Fold[
ToContravariantIndex[ #1  ,Part[#2,1],Part[#2,2] ]&,w, Transpose[{ind ,toind}]];


(*............................................*)	
(*...............Index Derivative.............*)
(*............................................*)

Format[ IndexDerivative[ w_ , ind__] ] = HoldForm[ D[ w , ind ] ];
IndexDerivative[ 
	IndexDerivative[F_,u__], v__ ] := IndexDerivative[ F,u,v];

EvaluateIndexDerivative[w_]:=  w /. IndexDerivative -> IDXX /.IDXX->D;

IndexDerivative[ x_?NumericQ ,_] = 0; 

IDXX[u_?NumericQ,_]=0;

IDXX[u_,u_]=1;

IDXX[ F_[u_] , v_ ] := 
   IDXX[F[u],u] IDXX[u,v] /; (u=!=v && F=!=Tensor);
	
IDXX[ u_*v_ , z_] =  IDXX[u,z]v + u IDXX[v,z];

IDXX[ u_ + v_ , z_] =  IDXX[u,z] + IDXX[v,z];
		
IDXX[Power[u_,m_] ,v_ ] = m Power[u,m-1]IDXX[u,v]+ Power[u,m]Log[u]IDXX[m,v];

IDXX[ ArcCos[u_] , v_ ] = -1/Sqrt[1-u^2] IDXX[u,v];

IDXX[ ArcSin[u_] , v_ ] = 1/Sqrt[1-u^2] IDXX[u,v];

IDXX[ Sin[u_] , v_ ] = Cos[u] IDXX[u,v];

IDXX[ Cos[u_] , v_ ] = -Sin[u] IDXX[u,v];
		
IDXX[ u_ , x_ ]:=  0/; Position[u,Tensor]=={};

IDXX[u_ , x_, y__] := IDXX[ IDXX[u,x] , y ];
(*............................................*)
(*.............Covariant Derivative...........*)
(*............................................*)

Attributes[CovariantDerivative]={HoldFirst};
Attributes[CDCov]={HoldFirst};

Format[ CovariantDerivative[T_, x_ , ind_ , g_]  ]:= 
   Subscript[  T  ,SequenceForm[ ",",SequenceForm@@ind ] ];

Format[ 
CovariantDerivative[  CovariantDerivative[ T_, x_ , ind1_ , g_] , x_ , ind2_ , g_] 
		 ]:=   Subscript[  T  ,SequenceForm[ ",",SequenceForm@@{ind1,ind2} ] ];

Format[ 
CovariantDerivative[  
   CovariantDerivative[ 
        CovariantDerivative[ T_, x_ , ind0_ , g_], x_ , ind1_ , g_] , x_ , ind2_ , g_] 
		 ]:=   Subscript[  T  ,SequenceForm[ ",",SequenceForm@@{ind0,ind1,ind2} ] ];

(*.....................*)

CovariantDerivative[  w_ ,  x_ , j_List, g_ ] :=
     Fold[ CovariantDerivative[ #1 , x , #2 ,g ]& , w , j ];

CovariantDerivative[  w_Plus , x_ ,j_ ,g_ ]:= 
     CovariantDerivative[ # , x,j,g]& /@ w

CovariantDerivative[  w_*u_ ,  x_ , j_ , g_ ] := 
     CovariantDerivative[  w ,  x ,j , g   ]* u + 
    CovariantDerivative[  u ,  x ,j , g   ]* w ;

CovariantDerivative[  w_?NumericQ ,  Tensor[x_,{Void},{j_}], g_ ] = 0; 
CovariantDerivative[  Kronecker[i_,j_] ,  __ ] = 0;

CovariantDerivative[ Tensor[T_] , x_,i_,g_] := 
         IndexDerivative[ T , Tensor[x,{Void},{i}] ] ;

(*......................*)

CovariantDerivative[ Tensor[T_,sub_,sup_] , x_ , i_ , g_]:=
   Module[{DP,xx,ind},
	ind = Transpose[{ sub , sup , Range[Length[sub]] }];
	xx = Tensor[x,{Void},{i}];
   IndexDerivative[Tensor[T,sub,sup] , Tensor[x,{Void},{i}]  ] +
   Plus@@Map[ CDCov[ Tensor[T,sub,sup] ,#,xx,g]& ,ind]
		]/;AllNumberQ[{sub,sup,i}];

(*........................*)

CDCov[ Tensor[T_,subT_,supT_] , {sub_,Void,n_}, Tensor[x_,_,{i_}] ,g_]:=
	Module[{k},
      -ExpandDoubleIndex[ 
             
        Tensor[T,ReplacePart[subT,k,n],supT] * Christoffel2[g,x][{sub,i},{k}] 
        
             ,{k}]    ];
	
CDCov[ Tensor[T_,subT_,supT_] , {Void,sup_,n_},Tensor[x_,_,{i_}],g_]:=
        Module[{k},
        ExpandDoubleIndex[
             
      Tensor[T,subT, ReplacePart[supT,k,n] ]*Christoffel2[g,x][{i,k},{sup}]  

            ,{k}]     ];

(*............................................*)
(*.................Absolute Derivative........*)
(*............................................*)
Format[AbsoluteDerivative[w_,x_] ]= HoldForm[Dt[w,x]];
Format[AbsoluteDerivative[ AbsoluteDerivative[w_,x_] ,y_]] =  HoldForm[Dt[w,x,y]];
Format[AbsoluteDerivative[ AbsoluteDerivative[w_,x_] ,x_]] =  HoldForm[Dt[w,{x,2}]];

(*...................*)

AbsoluteDerivative[ u_?NumericQ , x_ ]:=  0
AbsoluteDerivative[ u_Symbol , x_ ]:= 
		 0/; MemberQ[ Attributes[u] ,Constant];

AbsoluteDerivative[u_ w_,x_]= 
  u AbsoluteDerivative[w,x]+w AbsoluteDerivative[u,x];

AbsoluteDerivative[u_ + w_,x_]= 
  AbsoluteDerivative[w,x] +  AbsoluteDerivative[u,x];

AbsoluteDerivative[ F_[u_] ,v_	]:= 
AbsoluteDerivative[F[u],u] AbsoluteDerivative[u,v] /; (u=!=v&&F=!=Tensor);

AbsoluteDerivative[Power[u_,m_] ,v_ ] = 
  m Power[u,m-1]AbsoluteDerivative[u,v] + 
    Power[u,m]Log[u]AbsoluteDerivative[m,v];

AbsoluteDerivative[ ArcCos[u_] , v_ ] = -1/Sqrt[1-u^2] AbsoluteDerivative[u,v];
AbsoluteDerivative[ ArcSin[u_] , v_ ] =  1/Sqrt[1-u^2] AbsoluteDerivative[u,v];
AbsoluteDerivative[ Cos[u_] , v_ ] = -Sin[u] AbsoluteDerivative[u,v]; 
AbsoluteDerivative[ Sin[u_] , v_ ] =  Cos[u] AbsoluteDerivative[u,v];

(*.................*)

AbsoluteDerivativeSimplifyRule :={

Sum[
IndexDerivative[  A_,  Tensor[x_,{Void},{k}] ] * 
        AbsoluteDerivative[ Tensor[x_,{Void},{k}] ,t_]
			,{k,1,NDim}]     ->	AbsoluteDerivative[A,t]};

(*............................................*)
(*............Intrinsec Derivative............*)
(*............................................*)

Attributes[IntrinsecDerivative]={HoldFirst};

Format[IntrinsecDerivative[ w_ ,t_,{x_,g_}]]= 
  Times[ SequenceForm[ \[Delta] ,w ] , Power[ SequenceForm[\[Delta],t],-1]];

IntrinsecDerivative[ Tensor[T_,sub_,sup_] , t_ ,{x_,g_}]:= 
Module[{xx,ind,i},
	ind = Transpose[{ sub , sup , Range[Length[sub]] }];
	xx = Tensor[x,{Void},{i}];            
 AbsoluteDerivative[Tensor[T,sub,sup] , t  ] +
 ExpandDoubleIndex[
 AbsoluteDerivative[ xx , t  ]*Plus@@Map[ CDCov[ Tensor[T,sub,sup] , # , xx , g]& ,ind]
                   ,{i}]
			]/;AllNumberQ[{sub,sup}];

IntrinsecDerivative[ Tensor[F_] ,t_,{x_,g_}] = AbsoluteDerivative[F,t];

(*............................................*)
(*................Dummy Variables.............*)
(*............................................*)

DoubleVoidRule = 
  Tensor[P_,{subP__},{supP1___,Void,supP2___}] -> Tensor[P,{subP},{supP1,fXX[Void],supP2}] ;

DoubleIndexSimplifyRule01=	
a_.*Tensor[A_,{subA__},{supA1___,i_,supA2___}]*Tensor[B_,{subB1___,i_,subB2___},{supB___}] +
b_.*Tensor[A_,{subAA__},{supAA1___,j_,supAA2___}]*Tensor[B_,{subBB1___,j_,subBB2___},{supBB___}] :> 
a*Tensor[A,{subA},{supA1,i,supA2}]  Tensor[B,{subB1,i,subB2},{supB}]+
b*Tensor[A,{subAA},{supAA1,fXX[i],supAA2} ]*Tensor[B,{subBB1,i,subBB2},{supBB}] ;

(*............*)
DoubleIndexSimplifyRule02 =
a_.*Tensor[ A_,{subA1___,i_,subA2___},{supA1___,i_,supA2___}] + 
b_.*Tensor[ A_,{subAA1___,k_,subAA2___},{supAA1___,k_,supAA2___}]	:>
a*Tensor[ A,{subA1,i,subA2},{supA1,i,supA2}]+ 
b*Tensor[ A,{subAA1,fXX[i],subAA2},{supAA1,i,supAA2}]/;Length[{subA2}]===Length[{subAA2}];

(*..............*)
DoubleIndexSimplifyRule03=	
a_. Tensor[A_,{subA__},{supA1___,i_,supA2___}]*Tensor[B_,{subB1___,i_,subB2___},{supB___}]+
b_. Tensor[A_,{subAA__},{supAA1___,j_,supAA2___}]*Tensor[B_,{subBB1___,j_,subBB2___},{supBB___}]* 
Tensor[P_,{subP__},{supP1___,i_,supP2___}]*Tensor[Q_,{subQ1___,i_,subQ2___},{supQ___}] :> 
a*Tensor[A,{subA},{supA1,i,supA2}]  Tensor[B,{subB1,i,subB2},{supB}] +
b*Tensor[A,{subAA},{supAA1,fXX[i],supAA2} ]*Tensor[B,{subBB1,i,subBB2},{supBB}]*Tensor[P,{subP},{supP1,j,supP2}]*  
Tensor[Q,{subQ1,j,subQ2},{supQ}] ;

(*............*)
DoubleIndexSimplifyRule04 =
a_.*Tensor[ A_,{subA1___,i_,subA2___},{supA1___,i_,supA2___}]+ 
b_.* Tensor[ A_,{subAA1___,k_,subAA2___},{supAA1___,k_,supAA2___}]*Tensor[P_,{subP__},{supP1___,i_,supP2___}]*  
Tensor[Q_,{subQ1___,i_,subQ2___},{supQ___}]  :> a*Tensor[ A,{subA1,i,subA2},{supA1,i,supA2}]+ 
b*Tensor[ A,{subAA1,fXX[i],subAA2},{supAA1,i,supAA2}]* Tensor[P,{subP},{supP1,k ,supP2}]*  
Tensor[Q,{subQ1,k,subQ2},{supQ}]/;Length[{subA2}]===Length[{subAA2}];	



DummyVariableSimplify[w_]:=Module[{},
Expand[w]//.DoubleVoidRule//.DoubleIndexSimplifyRule04//.
DoubleIndexSimplifyRule02//.DoubleIndexSimplifyRule03//.DoubleIndexSimplifyRule01/.fXX->Times];

(*............................................*)
(*.............Christofel....................*)
(*............................................*)

Format[Christoffel1[g_,x_][{i_,j_},{k_}]]=  
  SequenceForm[ "[",i,j,",",k,"]","[",g,"]" ];

Format[Christoffel2[g_,x_][{i_,j_},{k_}]]=  
  SequenceForm[ Tensor[ Global`\[CapitalGamma],{i,j},{k} ],"[",g,"]" ];

Christoffel1[g_,x_][{i_,j_},{k_}]:= Module[{p},
x[p_] = Tensor[x,{Void},{p}];
		
1/2( IndexDerivative[  Tensor[g,{i,k},{Void,Void}]  , x[j]  ] +
     IndexDerivative[  Tensor[g,{j,k},{Void,Void}] , x[i]  ]-
     IndexDerivative[  Tensor[g,{i,j},{Void,Void}] , x[k]  ] )		
		]/;IntegerQ[i]&&IntegerQ[j]&&IntegerQ[k];

(*..............*)

Christoffel2[g_,x_][{i_,j_},{k_}]:= Module[{aa},
ExpandDoubleIndex[
	Tensor[g,{Void,Void},{k,aa}] Christoffel1[g,x][{i,j},{aa}]      
       ,{aa}]
		]/;IntegerQ[i]&&IntegerQ[j]&&IntegerQ[k];

(*............................................*)
(*.............Member Functions...............*)
(*............................................*)

FreeTensorQ[w_]:= Position[w,Tensor]=={};
FreeAllQ[w_,i_]:= Position[w,i] == {};
FreeAllQ[w_,i_List]:=  And@@Map[FreeAllQ[w,#]&,i] ;
MemberAllQ[w_,in_]:= Position[w,in]!= {};
MemberAllQ[w_,in_List]:= And@@Map[ MemberAllQ[w,#]& , in ];
AllNumberQ[w_]:= NumericQ[  Apply[Plus,Flatten[w]/.Void->0]    ];	

(*..........................................*)
(*..............Protect.....................*)
(*..........................................*)

End[]
EndPackage[]