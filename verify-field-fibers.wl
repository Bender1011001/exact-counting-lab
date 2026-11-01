(* Independent full symbolic H(T) and Taylor expansions from exact inputs. *)
ClearAll["Global`*"];
SetDirectory[DirectoryName[$InputFileName]];
assert[test_, message_] := If[!TrueQ[test], Print["FAIL: ",message]; Exit[1]];
records=Import["field-wolfram-input.json","RawJSON"];
checks=0; results={};
Do[
  forms=Map[Function[form,Total[Map[Function[entry,
    x[entry["v"]] Sum[(entry["cs"][[j,1]]/entry["cs"][[j,2]]) T^(j-1),
                       {j,Length[entry["cs"]]}]],form]]],record["forms"],{2}];
  clauses=Table[1-Product[If[TrueQ[record["occ"][[e,k,2]]],1-forms[[e,k]],forms[[e,k]]],
                         {k,Length[forms[[e]]]}],{e,Length[forms]}];
  vars=Table[x[v],{v,0,record["n"]-1}];
  polynomial=Expand[Total[Table[(Times@@clauses)/.Thread[vars->bits],
                               {bits,Tuples[{0,1},record["n"]]}]]];
  Do[
    p=query["p"]; q=query["q"]; t=query["t"];
    cs=PadRight[CoefficientList[Normal[Series[polynomial/.T->t+z,{z,0,q}]],z],q+1];
    residue=Map[Mod[Numerator[#] PowerMod[Denominator[#],-1,p],p]&,cs];
    assert[residue===query["python"],{record["n"],t,p,residue,query["python"]}];
    checks++;
    ,{query,record["queries"]}];
  AppendTo[results,<|"n"->record["n"],"degree"->Exponent[polynomial,T],
                    "count_at_zero"->(polynomial/.T->0)|>];
  ,{record,records}];
report=<|"kernel"->$Version,"coefficient_vector_checks"->checks,
          "scalar_coefficients_checked"->3 checks,"passed"->True,"cases"->results|>;
Export["field-wolfram-results.json",report,"RawJSON"];
Print[ExportString[report,"RawJSON"]];
