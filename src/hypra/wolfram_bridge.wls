(* Hypra protocol 1. Own bridge only; HypExp/HPL are externally installed.
   Input expressions are JSON ASTs, never Wolfram source. No ToExpression. *)
$HistoryLength = 0;
$MaxExtraPrecision = 2000;
ClearAll[fail, decode, encode, integer, loadPackage, runRequest, node, argsOf, requireSeriesDepth];
fail[s_String] := Throw[s, "HypraError"];
integer[s_String] := Module[{negative, digits},
  If[!StringMatchQ[s, RegularExpression["-?(0|[1-9][0-9]*)"]], fail["Invalid integer"]];
  negative = StringStartsQ[s, "-"]; digits = If[negative, StringDrop[s, 1], s];
  If[StringLength[digits] > 100000, fail["Integer too large"]];
  If[negative, -1, 1] FromDigits[digits]
];
integer[_] := fail["Integer value must be a string"];
node[h_String, a_List] := <|"head" -> h, "args" -> a|>;
argsOf[n_Association] := Module[{a = Lookup[n, "args", Missing[]]},
  If[!ListQ[a], fail["AST args must be a list"]]; decode /@ a
];
$symbolsByName = <||>; $nodesBySymbol = <||>; $symbolAssumptions = {};
decode[n_Association] := Module[{h, a, name, s, facts, num, den},
  h = Lookup[n, "head", Missing[]];
  Switch[h,
    "Integer", integer[Lookup[n, "value", Missing[]]],
    "Rational", num = integer[Lookup[n, "numerator", Missing[]]];
      den = integer[Lookup[n, "denominator", Missing[]]];
      If[den == 0, fail["Zero denominator"]]; num/den,
    "Symbol", name = Lookup[n, "name", Missing[]];
      If[!StringQ[name] || !StringMatchQ[name, RegularExpression["[A-Za-z][A-Za-z0-9_]*"]], fail["Invalid symbol name"]];
      If[KeyExistsQ[$symbolsByName, name],
        s = $symbolsByName[name];
        If[Lookup[$nodesBySymbol[SymbolName[s]], "assumptions", <||>] =!= Lookup[n, "assumptions", <||>], fail["Conflicting assumptions for same symbol name"]];
        Return[s]];
      s = Unique["HypraInput`v"];
      AssociateTo[$symbolsByName, name -> s]; AssociateTo[$nodesBySymbol, SymbolName[s] -> n];
      facts = Lookup[n, "assumptions", <||>];
      If[!AssociationQ[facts], fail["Symbol assumptions must be an object"]];
      KeyValueMap[Function[{key, value},
        If[!MemberQ[{True, False}, value], fail["Symbol facts must be boolean"]];
        With[{condition = Switch[key,
          "real", Element[s, Reals], "integer", Element[s, Integers],
          "rational", Element[s, Rationals], "complex", Element[s, Complexes],
          "positive", s > 0, "negative", s < 0, "nonnegative", s >= 0,
          "nonpositive", s <= 0, "zero", s == 0, "nonzero", s != 0,
          "even", Element[s/2, Integers], "odd", Element[(s - 1)/2, Integers],
          _, True]},
          If[condition =!= True, AppendTo[$symbolAssumptions, If[value, condition, Not[condition]]]]
      ]], facts]; s,
    "Pi", Pi, "E", E, "ImaginaryUnit", I, "EulerGamma", EulerGamma,
    "Catalan", Catalan, "True", True, "False", False,
    "Add", Plus @@ argsOf[n], "Mul", Times @@ argsOf[n],
    "List", argsOf[n],
    "Pow", a = argsOf[n]; If[Length[a] != 2, fail["Pow arity"]]; Power @@ a,
    "Gamma", a = argsOf[n]; If[Length[a] != 1, fail["Gamma arity"]]; Gamma @@ a,
    "Log", a = argsOf[n]; If[Length[a] != 1, fail["Log arity"]]; Log @@ a,
    "Exp", a = argsOf[n]; If[Length[a] != 1, fail["Exp arity"]]; Exp @@ a,
    "PolyLog", a = argsOf[n]; If[Length[a] != 2, fail["PolyLog arity"]]; PolyLog @@ a,
    "PolyGamma", a = argsOf[n]; If[Length[a] != 2, fail["PolyGamma arity"]]; PolyGamma @@ a,
    "Zeta", a = argsOf[n]; If[!MemberQ[{1, 2}, Length[a]], fail["Zeta arity"]]; Zeta @@ a,
    "HypergeometricPFQ", a = argsOf[n];
      If[Length[a] != 3 || !ListQ[a[[1]]] || !ListQ[a[[2]]], fail["HypergeometricPFQ arity"]];
      HypergeometricPFQ @@ a,
    "HPL", a = argsOf[n]; If[Length[a] != 2 || !ListQ[a[[1]]] || !AllTrue[a[[1]], MemberQ[{-1, 0, 1}, #] &], fail["HPL requires expanded -1,0,1 alphabet"]]; HPL`HPL @@ a,
    "And", And @@ argsOf[n], "Or", Or @@ argsOf[n],
    "Not", a = argsOf[n]; If[Length[a] != 1, fail["Not arity"]]; Not @@ a,
    "Equal", Equal @@ argsOf[n], "Unequal", Unequal @@ argsOf[n],
    "Less", Less @@ argsOf[n], "LessEqual", LessEqual @@ argsOf[n],
    "Greater", Greater @@ argsOf[n], "GreaterEqual", GreaterEqual @@ argsOf[n],
    _, fail["Unsupported AST head"]
  ]
];
decode[_] := fail["AST node must be an object"];
encode[e_] := Module[{h = Head[e], compressed, expanded, a},
  Which[
    IntegerQ[e], <|"head" -> "Integer", "value" -> ToString[e, InputForm]|>,
    h === Rational, <|"head" -> "Rational", "numerator" -> ToString[Numerator[e], InputForm], "denominator" -> ToString[Denominator[e], InputForm]|>,
    e === Pi, <|"head" -> "Pi"|>, e === E, <|"head" -> "E"|>,
    e === I, <|"head" -> "ImaginaryUnit"|>, e === EulerGamma, <|"head" -> "EulerGamma"|>,
    e === Catalan, <|"head" -> "Catalan"|>,
    e === True, <|"head" -> "True"|>, e === False, <|"head" -> "False"|>,
    h === Symbol, If[!KeyExistsQ[$nodesBySymbol, SymbolName[e]], fail["Unknown output symbol"]]; $nodesBySymbol[SymbolName[e]],
    h === Complex, node["Add", {encode[Re[e]], node["Mul", {encode[Im[e]], <|"head" -> "ImaginaryUnit"|>}]}],
    h === Plus, node["Add", encode /@ (List @@ e)],
    h === Times, node["Mul", encode /@ (List @@ e)],
    h === Power, node["Pow", encode /@ (List @@ e)],
    h === List, node["List", encode /@ e],
    h === Gamma, node["Gamma", encode /@ (List @@ e)],
    h === Log, node["Log", encode /@ (List @@ e)],
    h === PolyLog, node["PolyLog", encode /@ (List @@ e)],
    h === PolyGamma, node["PolyGamma", encode /@ (List @@ e)],
    h === Zeta, node["Zeta", encode /@ (List @@ e)],
    h === HPL`HPL, a = List @@ e;
      If[Length[a] != 2 || !ListQ[a[[1]]] || !AllTrue[a[[1]], IntegerQ], fail["Unsupported HPL alphabet"]];
      compressed = a[[1]];
      expanded = Flatten[If[# == 0, {0}, Join[ConstantArray[0, Abs[#] - 1], {Sign[#]}]] & /@ compressed];
      node["HPL", {node["List", encode /@ expanded], encode[a[[2]]]}],
    True, fail["Unsupported or unevaluated output head: " <> ToString[h, InputForm]]
  ]
];
loadPackage[name_String, configured_] := Module[{path, candidates},
  path = If[StringQ[configured], ExpandFileName[configured], FindFile[name <> "`"]];
  If[StringQ[path] && DirectoryQ[path],
    candidates = {FileNameJoin[{path, name <> ".m"}], FileNameJoin[{path, "Kernel", "init.m"}], FileNameJoin[{path, "init.m"}]};
    path = SelectFirst[candidates, FileExistsQ, $Failed]
  ];
  If[!StringQ[path] || !FileExistsQ[path], fail["External package not found: " <> name]];
  If[Quiet[Check[Get[path]; True, False]] =!= True, fail["External package load failed: " <> name]];
  path
];
(* Never silently drop a shallower O-term. A parameter-dependent HypExp result
   must retain explicit SeriesData remainder evidence even if the truncated
   output is constant. Only originally eps-free atoms are exempt. *)
requireSeriesDepth[expr_, eps_Symbol, order_Integer, mustHaveRemainder_] := Module[{parts, data},
  parts = Cases[expr, sd_SeriesData :> sd, {0, Infinity}];
  If[parts === {} && TrueQ[mustHaveRemainder], fail["HypExp returned no explicit remainder evidence"]];
  Do[
    data = List @@ part;
    If[Length[data] != 6 || data[[1]] =!= eps || data[[2]] =!= 0 || !IntegerQ[data[[5]]] || !IntegerQ[data[[6]]] || data[[6]] != 1 || data[[5]] < order + 1,
      fail["Insufficient or non-Laurent series remainder"]], {part, parts}];
  True
];
runRequest[r_Association] := Module[{packages, hp, hl, meta, op, expression, eps, through, internal, assumptions, expanded, series, polynomial, hAtoms, replacements, converted},
  If[Lookup[r, "protocol", 0] =!= 1, fail["Protocol mismatch"]];
  packages = Lookup[r, "packages", <||>];
  hl = loadPackage["HPL", Lookup[packages, "hpl", Null]];
  hp = loadPackage["HypExp", Lookup[packages, "hypexp", Null]];
  If[Length[DownValues[HypExp`HypExp]] == 0, fail["HypExp entry point unavailable"]];
  If[Length[DownValues[HPL`HPL]] == 0, fail["HPL entry point unavailable"]];
  meta = <|"wolfram_version" -> $Version, "system_id" -> $SystemID,
    "hypexp_version" -> If[ValueQ[HypExp`$HypExpVersion], ToString[HypExp`$HypExpVersion, InputForm], "unreported (identified by source fingerprint)"],
    "hpl_version" -> If[ValueQ[HPL`$HPLVersion], ToString[HPL`$HPLVersion, InputForm], "unreported (identified by source fingerprint)"],
    "hypexp_path" -> hp, "hpl_path" -> hl,
    "hypexp_sha256" -> IntegerString[FileHash[hp, "SHA256"], 16, 64],
    "hpl_sha256" -> IntegerString[FileHash[hl, "SHA256"], 16, 64]|>;
  op = Lookup[r, "operation", ""];
  If[op === "info", Return[<|"protocol" -> 1, "status" -> "ok", "backend" -> meta|>]];
  If[op =!= "expand", fail["Unknown operation"]];
  through = Lookup[r, "through", Missing[]]; internal = Lookup[r, "internal_order", Missing[]];
  If[!IntegerQ[through] || !IntegerQ[internal] || internal < Max[0, through], fail["Invalid expansion order"]];
  expression = decode[Lookup[r, "expression", Missing[]]];
  eps = decode[Lookup[r, "epsilon", Missing[]]];
  If[Head[eps] =!= Symbol, fail["Regulator must be a symbol"]];
  assumptions = decode[Lookup[r, "assumptions", <|"head" -> "True"|>]];
  assumptions = And[assumptions, And @@ $symbolAssumptions];
  If[assumptions === False || TrueQ[Quiet[Refine[assumptions] === False]], fail["Contradictory assumptions"]];
  hAtoms = DeleteDuplicates[Cases[expression, _Hypergeometric2F1 | _HypergeometricPFQ, {0, Infinity}]];
  replacements = Table[
    converted = Quiet[Check[Assuming[assumptions, HypExp`HypExp[h, eps, internal]], $Failed]];
    If[converted === $Failed || !FreeQ[converted, _HypergeometricPFQ | _Hypergeometric2F1 | _HypExp`HypExp], fail["HypExp could not expand a supported atom"]];
    requireSeriesDepth[converted, eps, internal, !FreeQ[h, eps]];
    h -> Normal[converted], {h, hAtoms}];
  expanded = expression /. replacements;
  series = Quiet[Check[Assuming[assumptions, Series[expanded, {eps, 0, through}]], $Failed]];
  If[series === $Failed || !FreeQ[series, _Series | _Derivative | _ConditionalExpression | _DirectedInfinity | Indeterminate], fail["Incomplete or ambiguous Laurent expansion"]];
  requireSeriesDepth[series, eps, through, False];
  polynomial = Normal[series];
  If[!FreeQ[polynomial, _SeriesData | _HypergeometricPFQ | _Hypergeometric2F1], fail["Unresolved series output"]];
  <|"protocol" -> 1, "status" -> "ok", "backend" -> meta, "result" -> encode[polynomial],
    "complete" -> True, "achieved_through" -> through, "remainder_order" -> through + 1,
    "order_basis" -> "hypra-pole-bound-plus-hypexp-series",
    "certification" -> "Algorithmic truncation contract; not an independent proof of HypExp or endpoint integration",
    "internal_order" -> internal|>
];
Module[{argv, request, result, responsePath},
  argv = $ScriptCommandLine;
  If[Length[argv] < 2, Exit[2]];
  responsePath = Last[argv];
  request = Quiet[Check[Import[argv[[-2]], "RawJSON"], $Failed]];
  result = Catch[If[!AssociationQ[request], fail["Invalid JSON request"]]; runRequest[request], "HypraError", ( <|"protocol" -> 1, "status" -> "error", "error" -> #1|> & )];
  If[!AssociationQ[result], result = <|"protocol" -> 1, "status" -> "error", "error" -> "Unexpected bridge failure"|>];
  If[Quiet[Check[Export[responsePath, result, "RawJSON"]; True, False]] =!= True, Exit[3]];
];
Exit[0];
