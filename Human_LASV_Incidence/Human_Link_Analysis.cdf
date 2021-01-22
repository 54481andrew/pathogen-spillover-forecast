(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 11.2' *)

(***************************************************************************)
(*                                                                         *)
(*                                                                         *)
(*  Under the Wolfram FreeCDF terms of use, this file and its content are  *)
(*  bound by the Creative Commons BY-SA Attribution-ShareAlike license.    *)
(*                                                                         *)
(*        For additional information concerning CDF licensing, see:        *)
(*                                                                         *)
(*         www.wolfram.com/cdf/adopting-cdf/licensing-options.html         *)
(*                                                                         *)
(*                                                                         *)
(***************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1088,         20]
NotebookDataLength[     39125,       1119]
NotebookOptionsPosition[     32992,        990]
NotebookOutlinePosition[     33537,       1011]
CellTagsIndexPosition[     33494,       1008]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Human Lassa Case Rate", "Title",ExpressionUUID->"e8bd921e-cf66-4ebd-9513-78eb872ef4e0"],

Cell["\<\
This code shows the analyses performed on the system of ordinary differential \
equations that describe the spillover of Lassa virus into humans. \
\>", "Text",ExpressionUUID->"62b26a79-91fc-4292-bade-fd7745f6e7b1"],

Cell[CellGroupData[{

Cell["Define model", "Subchapter",ExpressionUUID->"bf5dc25c-c5a4-43fc-9614-1e3c3a1a6627"],

Cell["Define the system of ODEs that describe human infection. ", "Text",ExpressionUUID->"afbdb8fc-6b8e-47bd-a8cb-dbd223d92b9a"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"dsdt", " ", "=", " ", 
   RowBox[{"b", " ", "-", " ", 
    RowBox[{"d", " ", "s"}], " ", "-", " ", 
    RowBox[{"F", " ", "s"}], " ", "+", " ", 
    RowBox[{"\[Lambda]", " ", "r"}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"didt", " ", "=", " ", 
   RowBox[{
    RowBox[{"F", " ", "s"}], " ", "-", " ", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"d", " ", "+", " ", "\[Gamma]"}], ")"}], " ", "i"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"drdt", " ", "=", " ", 
   RowBox[{
    RowBox[{
     RowBox[{"(", 
      RowBox[{"1", "-", "\[Mu]"}], ")"}], " ", "\[Gamma]", " ", "i"}], " ", 
    "-", " ", 
    RowBox[{"d", " ", "r"}], " ", "-", " ", 
    RowBox[{"\[Lambda]", " ", "r"}]}]}], ";"}]}], "Input",ExpressionUUID->\
"c2034f9e-36c7-41d2-b8d0-3beb84a8be0b"]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[{
 "Derive relationship between ",
 StyleBox["seroprevalence",
  IgnoreSpellCheck->True],
 " and new case rate"
}], "Subchapter",ExpressionUUID->"4b2b9d26-da4c-4491-80b3-de7bbda8e7fc"],

Cell["Find the equilibrium.", "Text",ExpressionUUID->"0e34a607-f81a-4810-ad25-61c41f260938"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"eq", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"Solve", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"dsdt", "\[Equal]", "0"}], ",", " ", 
       RowBox[{"didt", "\[Equal]", "0"}], ",", " ", 
       RowBox[{"drdt", "\[Equal]", "0"}]}], "}"}], ",", " ", 
     RowBox[{"{", 
      RowBox[{"s", ",", "i", ",", "r"}], "}"}]}], "]"}], "]"}]}]], "Input",Exp\
ressionUUID->"93f7df5f-3eff-4d29-bb6c-999fff2366be"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"s", "\[Rule]", 
     FractionBox[
      RowBox[{"b", " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", "\[Gamma]"}], ")"}], " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
      RowBox[{
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "F"}], ")"}], " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "\[Gamma]"}], ")"}]}], "+", 
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}], 
       "+", 
       RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", "\[Mu]"}]}]]}], 
    ",", 
    RowBox[{"i", "\[Rule]", 
     FractionBox[
      RowBox[{"b", " ", "F", " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
      RowBox[{
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "F"}], ")"}], " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "\[Gamma]"}], ")"}]}], "+", 
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}], 
       "+", 
       RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", "\[Mu]"}]}]]}], 
    ",", 
    RowBox[{"r", "\[Rule]", 
     RowBox[{"-", 
      FractionBox[
       RowBox[{"b", " ", "F", " ", "\[Gamma]", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], 
       RowBox[{
        RowBox[{"d", " ", 
         RowBox[{"(", 
          RowBox[{"d", "+", "F"}], ")"}], " ", 
         RowBox[{"(", 
          RowBox[{"d", "+", "\[Gamma]"}], ")"}]}], "+", 
        RowBox[{"d", " ", 
         RowBox[{"(", 
          RowBox[{"d", "+", "F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}],
         "+", 
        RowBox[{
        "F", " ", "\[Gamma]", " ", "\[Lambda]", " ", "\[Mu]"}]}]]}]}]}], 
   "}"}], "}"}]], "Output",ExpressionUUID->"a7b364aa-15c1-442e-b579-\
621e4b233fd8"]
}, Open  ]],

Cell["\<\
Find seroprevalence at equilibrium (this is the data we have to work with). \
Note that this requires calculating the steady state population size in the \
presence of LASV virulence. \
\>", "Text",ExpressionUUID->"c2ea8db9-eed4-4832-8d72-a16e3823ce67"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Pop", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"s", " ", "+", " ", "i", " ", "+", " ", "r"}], " ", ")"}], "/.", 
    "eq"}], "]"}]}]], "Input",ExpressionUUID->"32742f3e-8886-422f-a2a1-\
20784ecfef8c"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox[
   RowBox[{"b", " ", 
    RowBox[{"(", 
     RowBox[{
      SuperscriptBox["d", "2"], "+", 
      RowBox[{"\[Gamma]", " ", "\[Lambda]"}], "+", 
      RowBox[{"d", " ", 
       RowBox[{"(", 
        RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "+", 
      RowBox[{"F", " ", 
       RowBox[{"(", 
        RowBox[{"\[Gamma]", "+", "\[Lambda]", "-", 
         RowBox[{"\[Gamma]", " ", "\[Mu]"}]}], ")"}]}]}], ")"}]}], 
   RowBox[{
    RowBox[{"d", " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "F"}], ")"}], " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "\[Gamma]"}], ")"}]}], "+", 
    RowBox[{"d", " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}], 
    "+", 
    RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", "\[Mu]"}]}]], 
  "}"}]], "Output",ExpressionUUID->"8afa6e5c-1107-40ea-8bcb-dbd1bc9d975a"]
}, Open  ]],

Cell["\<\
Variable [count] refers to the population size in a pixel of the population \
raster.  Solve for the birth rate [b] for a given value of [count].\
\>", "Text",ExpressionUUID->"0632205d-9488-4c23-8efd-ab33142bad23"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"Pop", "\[Equal]", "count"}], ",", " ", "b"}], "]"}]}]], "Input",E\
xpressionUUID->"f26625be-e98c-466f-874e-df94195b5e1e"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{"b", "\[Rule]", 
    FractionBox[
     RowBox[{"count", " ", 
      RowBox[{"(", 
       RowBox[{
        SuperscriptBox["d", "3"], "+", 
        RowBox[{
         SuperscriptBox["d", "2"], " ", "F"}], "+", 
        RowBox[{
         SuperscriptBox["d", "2"], " ", "\[Gamma]"}], "+", 
        RowBox[{"d", " ", "F", " ", "\[Gamma]"}], "+", 
        RowBox[{
         SuperscriptBox["d", "2"], " ", "\[Lambda]"}], "+", 
        RowBox[{"d", " ", "F", " ", "\[Lambda]"}], "+", 
        RowBox[{"d", " ", "\[Gamma]", " ", "\[Lambda]"}], "+", 
        RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", "\[Mu]"}]}], 
       ")"}]}], 
     RowBox[{
      SuperscriptBox["d", "2"], "+", 
      RowBox[{"d", " ", "F"}], "+", 
      RowBox[{"d", " ", "\[Gamma]"}], "+", 
      RowBox[{"F", " ", "\[Gamma]"}], "+", 
      RowBox[{"d", " ", "\[Lambda]"}], "+", 
      RowBox[{"F", " ", "\[Lambda]"}], "+", 
      RowBox[{"\[Gamma]", " ", "\[Lambda]"}], "-", 
      RowBox[{"F", " ", "\[Gamma]", " ", "\[Mu]"}]}]]}], "}"}], 
  "}"}]], "Output",ExpressionUUID->"d4a99042-dbfe-4563-a46b-ae24e1f002b4"]
}, Open  ]],

Cell["Calculate seroprevalence in terms of model parameters. ", "Text",ExpressionUUID->"6b092b83-ee5d-45b4-bedc-f12e405ef8cd"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Rstar", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"r", "/", 
      RowBox[{"(", "Pop", ")"}]}], "/.", "eq"}], "]"}], "//", 
   "Flatten"}]}]], "Input",ExpressionUUID->"2bad1404-bbfe-4377-b9c3-\
9be5d4530b36"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"-", 
   FractionBox[
    RowBox[{"F", " ", "\[Gamma]", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], 
    RowBox[{
     SuperscriptBox["d", "2"], "+", 
     RowBox[{"\[Gamma]", " ", "\[Lambda]"}], "+", 
     RowBox[{"d", " ", 
      RowBox[{"(", 
       RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "+", 
     RowBox[{"F", " ", 
      RowBox[{"(", 
       RowBox[{"\[Gamma]", "+", "\[Lambda]", "-", 
        RowBox[{"\[Gamma]", " ", "\[Mu]"}]}], ")"}]}]}]]}], "}"}]], "Output",E\
xpressionUUID->"ce9910af-960e-4069-8417-dcc388140caf"]
}, Open  ]],

Cell["\<\
Solve for the force of infection, F, in terms of seroprevalence.\
\>", "Text",ExpressionUUID->"7d9d9192-4ed3-456f-9d77-dde61685e10a"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"foi", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{"Solve", "[", 
     RowBox[{
      RowBox[{"hsero", " ", "\[Equal]", " ", "Rstar"}], ",", " ", "F"}], 
     "]"}], "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"9316df9a-\
cd22-498f-b8b6-6435e808a1e8"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"F", "\[Rule]", 
   RowBox[{"-", 
    FractionBox[
     RowBox[{"hsero", " ", 
      RowBox[{"(", 
       RowBox[{"d", "+", "\[Gamma]"}], ")"}], " ", 
      RowBox[{"(", 
       RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
     RowBox[{
      RowBox[{"d", " ", "hsero"}], "+", 
      RowBox[{"hsero", " ", "\[Lambda]"}], "+", 
      RowBox[{"\[Gamma]", " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "+", "hsero", "+", "\[Mu]", "-", 
         RowBox[{"hsero", " ", "\[Mu]"}]}], ")"}]}]}]]}]}], "}"}]], "Output",E\
xpressionUUID->"07582294-4d3f-4d5e-9212-829e86d47d28"]
}, Open  ]],

Cell["\<\
Rewrite the solution for [b], defined in [bsolve], with the force of \
infection terms F replaced with the expression above\
\>", "Text",ExpressionUUID->"7eac6851-b102-447b-a630-3bbd1ad5a950"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"bsolve", "/.", "foi"}], "]"}]}]], "Input",ExpressionUUID->\
"2dbd6c6f-97ed-4841-bd4f-2a4511237572"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{"b", "\[Rule]", 
    RowBox[{"-", 
     FractionBox[
      RowBox[{"count", " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", 
         RowBox[{"d", " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{"-", "1"}], "+", "hsero"}], ")"}], " ", "\[Mu]"}], "+", 
         RowBox[{"hsero", " ", "\[Lambda]", " ", "\[Mu]"}]}], ")"}]}], 
      RowBox[{
       RowBox[{"-", "1"}], "+", "\[Mu]"}]]}]}], "}"}], "}"}]], "Output",Expres\
sionUUID->"1ffee6b8-8aed-4a56-9280-bfacbb782b38"]
}, Open  ]],

Cell["\<\
The rate of new cases of Lassa is F*s. Calculate this, and rewrite everything \
in terms of seroprevalence.\
\>", "Text",ExpressionUUID->"a2cbf31a-22dc-4e2e-91e4-a8aa58df158c"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"NewCaseRate", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{"F", "/.", "foi"}], ")"}], "*", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"s", "/.", "eq"}], "/.", "foi"}], ")"}]}], "/.", "bsolve"}], 
    "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"e344f092-8501-4b81-\
886a-b595b9388254"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"-", 
   FractionBox[
    RowBox[{"count", " ", "hsero", " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "\[Gamma]"}], ")"}], " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
    RowBox[{"\[Gamma]", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}]]}], "}"}]], "Output",Expres\
sionUUID->"f3f64f48-e396-4b0b-b0b1-79ef760b751d"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Compare case rate estimates for virulence vs no virulence", "Subchapter",ExpressionUUID->"2572d4f9-224b-49c7-b5e6-86a655b367f8"],

Cell["\<\
Calculate the ratio of cases with, versus without, virulence. \
\>", "Text",ExpressionUUID->"be57a033-1922-4a08-bbdc-0de7a05bab83"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"\[Rho]", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"NewCaseRate", "/.", "eq"}], ")"}], "/", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"NewCaseRate", "/.", "eq"}], "/.", 
       RowBox[{"\[Mu]", "\[Rule]", "0"}]}], ")"}]}], " ", "]"}], "//", 
   "Flatten"}]}]], "Input",ExpressionUUID->"ddbc0717-8bfe-4f89-8870-\
e1b0df342e55"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox["1", 
   RowBox[{"1", "-", "\[Mu]"}]], "}"}]], "Output",ExpressionUUID->"95652392-\
fc06-422e-929b-01cb3848186c"]
}, Open  ]],

Cell[TextData[{
 "In words, including virulence increases the estimated rate of LASV cases  \
by ",
 Cell[BoxData[
  FormBox[
   RowBox[{"\[Rho]", " ", "=", " ", 
    SuperscriptBox[
     RowBox[{"(", 
      RowBox[{"1", "-", "\[Mu]"}], ")"}], 
     RowBox[{"-", "1"}]]}], TraditionalForm]],ExpressionUUID->
  "9e5ea74e-9487-4f64-8378-03e50bf6cb54"],
 ".  "
}], "Text",ExpressionUUID->"56599fca-1177-4675-bbac-69686eebea23"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"\[Rho]", "/.", 
  RowBox[{"\[Mu]", "\[Rule]", "0.02"}]}]], "Input",ExpressionUUID->"44a34a14-\
851a-42c1-9b40-2a97534f01d4"],

Cell[BoxData[
 RowBox[{"{", "1.0204081632653061`", "}"}]], "Output",ExpressionUUID->"b1c6d1b6-a6a2-465b-9000-1846f83f324b"]
}, Open  ]],

Cell["\<\
When \[Mu] = 0.02 (lower value from McCormick, Webb, 1987), increase is \
~1.02, meaning roughly 2% more cases occur. \
\>", "Text",ExpressionUUID->"4ffb9d56-9187-42aa-b8c2-38fe3f9b021b"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Compare estimates for reinfection vs no reinfection ", "Subchapter",ExpressionUUID->"20a2bd3f-2233-4454-a614-34d6bd8e0096"],

Cell["\<\
Calculate the ratio of cases with, versus without, reinfection. \
\>", "Text",ExpressionUUID->"0a8428b2-4919-445f-aa88-bda270e09c68"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"NewCaseRate", "/.", "eq"}], ")"}], "/", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"NewCaseRate", "/.", "eq"}], "/.", 
      RowBox[{"\[Lambda]", "\[Rule]", "0"}]}], ")"}]}], " ", "]"}], "//", 
  "Flatten"}]], "Input",ExpressionUUID->"4a900637-e7e7-4298-a70d-\
bfd610027a57"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox[
   RowBox[{"d", "+", "\[Lambda]"}], "d"], "}"}]], "Output",ExpressionUUID->\
"95288e70-f4a3-4af8-a3f4-65ba976b910c"]
}, Open  ]],

Cell[TextData[{
 "Reinfection multiplies the case rate by ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"d", "+", "\[Lambda]"}], "d"], TraditionalForm]],ExpressionUUID->
  "e17836b2-dc3d-460c-a48b-6ed3a509f55a"],
 ". Calculate the effect of reinfection for a mean lifespan of 50 years, and \
reinfection rate \[Lambda] = 0.064 / yr (from McCormick, Webb, 1987)."
}], "Text",ExpressionUUID->"b7f09cc0-7381-4115-a932-bc9ce3065a41"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Pars", "=", 
   RowBox[{"{", 
    RowBox[{"d", "\[Rule]", " ", "0.02"}], "}"}]}], ";"}]], "Input",Expression\
UUID->"9e0c04ef-017e-4376-b057-25f1dd3b5d6c"],

Cell["\<\
Ratio of estimates of new cases with (\[Lambda] = 0.064) and without (\
\[Lambda] = 0) seroreversion is\
\>", "Text",ExpressionUUID->"fc10009b-48de-4c99-a925-f8bcb4e9c63b"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"(", 
    RowBox[{"NewCaseRate", "/.", 
     RowBox[{"\[Lambda]", "\[Rule]", "0.064"}]}], ")"}], "/", 
   RowBox[{"(", 
    RowBox[{"NewCaseRate", "/.", 
     RowBox[{"\[Lambda]", "\[Rule]", "0"}]}], ")"}]}], "/.", 
  "Pars"}]], "Input",ExpressionUUID->"3f93b19e-c657-4271-acc9-643d36b44d33"],

Cell[BoxData[
 RowBox[{"{", "4.2`", "}"}]], "Output",ExpressionUUID->"ff1522a4-b20d-4629-aa94-3373ac8c0311"]
}, Open  ]],

Cell["\<\
Including reinfection multiplies the case estimates by a factor of 420%; \
including virulence increases estimates by 2%. \
\>", "Text",ExpressionUUID->"0e3c9cd6-dd3c-4a3b-b512-c9497faa5d20"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Human LASV cases with undetectable antibody immune class", "Title",ExpressionUUID->"1420339a-6f1f-4312-b2d3-5bb90cd4c53d"],

Cell["\<\
At the request of a reviewer, we also investigate the possibility that a \
fraction, a, of recovered individuals with immunity (R) transition into a \
class with undetectable immunity (C). This might be the case, for example, if \
immunity were caused primarily by a T-cell response. into humans. \
\>", "Text",ExpressionUUID->"4dbfbd23-a124-414b-91f2-e290b5a4bbe5"],

Cell[CellGroupData[{

Cell["Define model", "Subchapter",ExpressionUUID->"632f3d9d-12ae-4ca9-9a62-f069668a78e2"],

Cell["\<\
Define the system of ODEs that describe human infection.  Note the additional \
class, c, and the rate term a. \
\>", "Text",ExpressionUUID->"c3f12661-cdd3-4d43-9101-cf97ceb7b28d"],

Cell[BoxData[{
 RowBox[{
  RowBox[{"dsdt", " ", "=", " ", 
   RowBox[{"b", " ", "-", " ", 
    RowBox[{"d", " ", "s"}], " ", "-", " ", 
    RowBox[{"F", " ", "s"}], " ", "+", " ", 
    RowBox[{"\[Lambda]", " ", 
     RowBox[{"(", 
      RowBox[{"1", " ", "-", " ", "a"}], ")"}], " ", "r"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"didt", " ", "=", " ", 
   RowBox[{
    RowBox[{"F", " ", "s"}], " ", "-", " ", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"d", " ", "+", " ", "\[Gamma]"}], ")"}], " ", "i"}]}]}], 
  ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"drdt", " ", "=", " ", 
   RowBox[{
    RowBox[{
     RowBox[{"(", 
      RowBox[{"1", "-", "\[Mu]"}], ")"}], " ", "\[Gamma]", " ", "i"}], " ", 
    "-", " ", 
    RowBox[{"d", " ", "r"}], " ", "-", " ", 
    RowBox[{"\[Lambda]", " ", "r"}]}]}], ";"}], "\[IndentingNewLine]", 
 RowBox[{
  RowBox[{"dcdt", " ", "=", " ", 
   RowBox[{
    RowBox[{"a", " ", "lam", " ", "r"}], " ", "-", " ", 
    RowBox[{"d", " ", "c"}]}]}], ";"}]}], "Input",ExpressionUUID->"dc2fd5cd-\
5006-4506-9277-f4fb4eb4b0a2"]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[{
 "Derive relationship between ",
 StyleBox["seroprevalence",
  IgnoreSpellCheck->True],
 " and new case rate"
}], "Subchapter",ExpressionUUID->"4d9ab01c-3488-43d4-a66c-a533d09ca690"],

Cell["\<\
The remainder of this analysis is identical to what was performed above. Note \
that the equilibrium values of s are smaller when a > 0. Find the \
equilibrium. \
\>", "Text",ExpressionUUID->"37f4c21b-25e7-4fea-9d04-208fcba24398"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"eq", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"Solve", "[", 
    RowBox[{
     RowBox[{"{", 
      RowBox[{
       RowBox[{"dsdt", "\[Equal]", "0"}], ",", " ", 
       RowBox[{"didt", "\[Equal]", "0"}], ",", " ", 
       RowBox[{"drdt", "\[Equal]", "0"}], ",", " ", 
       RowBox[{"dcdt", " ", "\[Equal]", "0"}]}], "}"}], ",", " ", 
     RowBox[{"{", 
      RowBox[{"s", ",", "i", ",", "r", ",", "c"}], "}"}]}], "]"}], 
   "]"}]}]], "Input",ExpressionUUID->"e4dafdce-a95a-4d56-9cac-056895578223"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{
    RowBox[{"s", "\[Rule]", 
     FractionBox[
      RowBox[{"b", " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", "\[Gamma]"}], ")"}], " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
      RowBox[{
       SuperscriptBox["d", "3"], "+", 
       RowBox[{"d", " ", "F", " ", "\[Gamma]"}], "+", 
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{"F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}], "+", 
       RowBox[{
        SuperscriptBox["d", "2"], " ", 
        RowBox[{"(", 
         RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "+", 
       RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", 
        RowBox[{"(", 
         RowBox[{"a", "+", "\[Mu]", "-", 
          RowBox[{"a", " ", "\[Mu]"}]}], ")"}]}]}]]}], ",", 
    RowBox[{"i", "\[Rule]", 
     FractionBox[
      RowBox[{"b", " ", "F", " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
      RowBox[{
       SuperscriptBox["d", "3"], "+", 
       RowBox[{"d", " ", "F", " ", "\[Gamma]"}], "+", 
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{"F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}], "+", 
       RowBox[{
        SuperscriptBox["d", "2"], " ", 
        RowBox[{"(", 
         RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "+", 
       RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", 
        RowBox[{"(", 
         RowBox[{"a", "+", "\[Mu]", "-", 
          RowBox[{"a", " ", "\[Mu]"}]}], ")"}]}]}]]}], ",", 
    RowBox[{"r", "\[Rule]", 
     RowBox[{"-", 
      FractionBox[
       RowBox[{"b", " ", "F", " ", "\[Gamma]", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], 
       RowBox[{
        SuperscriptBox["d", "3"], "+", 
        RowBox[{"d", " ", "F", " ", "\[Gamma]"}], "+", 
        RowBox[{"d", " ", 
         RowBox[{"(", 
          RowBox[{"F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}], "+", 
        RowBox[{
         SuperscriptBox["d", "2"], " ", 
         RowBox[{"(", 
          RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "+", 
        RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", 
         RowBox[{"(", 
          RowBox[{"a", "+", "\[Mu]", "-", 
           RowBox[{"a", " ", "\[Mu]"}]}], ")"}]}]}]]}]}], ",", 
    RowBox[{"c", "\[Rule]", 
     RowBox[{"-", 
      FractionBox[
       RowBox[{"a", " ", "b", " ", "F", " ", "lam", " ", "\[Gamma]", " ", 
        RowBox[{"(", 
         RowBox[{
          RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], 
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{
          SuperscriptBox["d", "3"], "+", 
          RowBox[{"d", " ", "F", " ", "\[Gamma]"}], "+", 
          RowBox[{"d", " ", 
           RowBox[{"(", 
            RowBox[{"F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}], "+", 
          RowBox[{
           SuperscriptBox["d", "2"], " ", 
           RowBox[{"(", 
            RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "+", 
          RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", 
           RowBox[{"(", 
            RowBox[{"a", "+", "\[Mu]", "-", 
             RowBox[{"a", " ", "\[Mu]"}]}], ")"}]}]}], ")"}]}]]}]}]}], "}"}], 
  "}"}]], "Output",ExpressionUUID->"fb9b51c4-37b7-4b8f-9f33-b039b1cfd462"]
}, Open  ]],

Cell["\<\
Find seroprevalence at equilibrium (this is the data we have to work with). \
Note that this requires calculating the steady state population size in the \
presence of LASV virulence. \
\>", "Text",ExpressionUUID->"9482acc8-dd24-4bc6-a542-e4099cd9ddd7"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Pop", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
     "s", " ", "+", " ", "i", " ", "+", " ", "r", " ", "+", " ", "c"}], " ", 
     ")"}], "/.", "eq"}], "]"}]}]], "Input",ExpressionUUID->"03e0ed29-5af3-\
49f1-ba5f-60d71cd49e16"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox[
   RowBox[{"b", " ", 
    RowBox[{"(", 
     RowBox[{
      SuperscriptBox["d", "3"], "+", 
      RowBox[{"d", " ", "\[Gamma]", " ", "\[Lambda]"}], "+", 
      RowBox[{
       SuperscriptBox["d", "2"], " ", 
       RowBox[{"(", 
        RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "-", 
      RowBox[{"a", " ", "F", " ", "lam", " ", "\[Gamma]", " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], "+", 
      RowBox[{"d", " ", "F", " ", 
       RowBox[{"(", 
        RowBox[{"\[Gamma]", "+", "\[Lambda]", "-", 
         RowBox[{"\[Gamma]", " ", "\[Mu]"}]}], ")"}]}]}], ")"}]}], 
   RowBox[{"d", " ", 
    RowBox[{"(", 
     RowBox[{
      SuperscriptBox["d", "3"], "+", 
      RowBox[{"d", " ", "F", " ", "\[Gamma]"}], "+", 
      RowBox[{"d", " ", 
       RowBox[{"(", 
        RowBox[{"F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}], "+", 
      RowBox[{
       SuperscriptBox["d", "2"], " ", 
       RowBox[{"(", 
        RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "+", 
      RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", 
       RowBox[{"(", 
        RowBox[{"a", "+", "\[Mu]", "-", 
         RowBox[{"a", " ", "\[Mu]"}]}], ")"}]}]}], ")"}]}]], "}"}]], "Output",\
ExpressionUUID->"7554bdac-e1c0-4011-9268-aa5f0fe2a240"]
}, Open  ]],

Cell["\<\
Variable [count] refers to the population size in a pixel of the population \
raster.  Solve for the birth rate [b] for a given value of [count].\
\>", "Text",ExpressionUUID->"994b8c52-78f9-4174-8b8d-6822cf9ef3d1"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"Pop", "\[Equal]", "count"}], ",", " ", "b"}], "]"}]}]], "Input",E\
xpressionUUID->"3a249cc3-f0a3-4e96-8918-5268e71da9ca"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{"b", "\[Rule]", 
    FractionBox[
     RowBox[{"count", " ", "d", " ", 
      RowBox[{"(", 
       RowBox[{
        SuperscriptBox["d", "3"], "+", 
        RowBox[{
         SuperscriptBox["d", "2"], " ", "F"}], "+", 
        RowBox[{
         SuperscriptBox["d", "2"], " ", "\[Gamma]"}], "+", 
        RowBox[{"d", " ", "F", " ", "\[Gamma]"}], "+", 
        RowBox[{
         SuperscriptBox["d", "2"], " ", "\[Lambda]"}], "+", 
        RowBox[{"d", " ", "F", " ", "\[Lambda]"}], "+", 
        RowBox[{"d", " ", "\[Gamma]", " ", "\[Lambda]"}], "+", 
        RowBox[{"a", " ", "F", " ", "\[Gamma]", " ", "\[Lambda]"}], "+", 
        RowBox[{"F", " ", "\[Gamma]", " ", "\[Lambda]", " ", "\[Mu]"}], "-", 
        RowBox[{
        "a", " ", "F", " ", "\[Gamma]", " ", "\[Lambda]", " ", "\[Mu]"}]}], 
       ")"}]}], 
     RowBox[{
      SuperscriptBox["d", "3"], "+", 
      RowBox[{
       SuperscriptBox["d", "2"], " ", "F"}], "+", 
      RowBox[{
       SuperscriptBox["d", "2"], " ", "\[Gamma]"}], "+", 
      RowBox[{"d", " ", "F", " ", "\[Gamma]"}], "+", 
      RowBox[{"a", " ", "F", " ", "lam", " ", "\[Gamma]"}], "+", 
      RowBox[{
       SuperscriptBox["d", "2"], " ", "\[Lambda]"}], "+", 
      RowBox[{"d", " ", "F", " ", "\[Lambda]"}], "+", 
      RowBox[{"d", " ", "\[Gamma]", " ", "\[Lambda]"}], "-", 
      RowBox[{"d", " ", "F", " ", "\[Gamma]", " ", "\[Mu]"}], "-", 
      RowBox[{"a", " ", "F", " ", "lam", " ", "\[Gamma]", " ", "\[Mu]"}]}]]}],
    "}"}], "}"}]], "Output",ExpressionUUID->"d2f76664-93ff-44cb-9e2a-\
c98c820d3499"]
}, Open  ]],

Cell["Calculate seroprevalence in terms of model parameters. ", "Text",ExpressionUUID->"92d03e44-da90-42a5-9ada-ca151ee398e8"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Rstar", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"r", "/", 
      RowBox[{"(", "Pop", ")"}]}], "/.", "eq"}], "]"}], "//", 
   "Flatten"}]}]], "Input",ExpressionUUID->"b5fb890b-ffe5-4ef2-b98f-\
313adab66bb9"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"-", 
   FractionBox[
    RowBox[{"d", " ", "F", " ", "\[Gamma]", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], 
    RowBox[{
     SuperscriptBox["d", "3"], "+", 
     RowBox[{"d", " ", "\[Gamma]", " ", "\[Lambda]"}], "+", 
     RowBox[{
      SuperscriptBox["d", "2"], " ", 
      RowBox[{"(", 
       RowBox[{"F", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}], "-", 
     RowBox[{"a", " ", "F", " ", "lam", " ", "\[Gamma]", " ", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], "+", 
     RowBox[{"d", " ", "F", " ", 
      RowBox[{"(", 
       RowBox[{"\[Gamma]", "+", "\[Lambda]", "-", 
        RowBox[{"\[Gamma]", " ", "\[Mu]"}]}], ")"}]}]}]]}], "}"}]], "Output",E\
xpressionUUID->"dcb14bbf-9d32-4dfd-8a21-f84fb858ed3c"]
}, Open  ]],

Cell["\<\
Solve for the force of infection, F, in terms of seroprevalence. Note that \
the equilibrium value of foi is made larger by the parameter a. \
\>", "Text",ExpressionUUID->"be2995b9-ab4b-4ce6-913f-538535ff91bf"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"foi", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{"Solve", "[", 
     RowBox[{
      RowBox[{"hsero", " ", "\[Equal]", " ", "Rstar"}], ",", " ", "F"}], 
     "]"}], "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"04d0dada-\
7e7e-4425-8d25-d7b112ff88f3"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"F", "\[Rule]", 
   RowBox[{"-", 
    FractionBox[
     RowBox[{"d", " ", "hsero", " ", 
      RowBox[{"(", 
       RowBox[{"d", "+", "\[Gamma]"}], ")"}], " ", 
      RowBox[{"(", 
       RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
     RowBox[{
      RowBox[{
       SuperscriptBox["d", "2"], " ", "hsero"}], "+", 
      RowBox[{"d", " ", "hsero", " ", "\[Lambda]"}], "-", 
      RowBox[{"a", " ", "hsero", " ", "lam", " ", "\[Gamma]", " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], "+", 
      RowBox[{"d", " ", "\[Gamma]", " ", 
       RowBox[{"(", 
        RowBox[{
         RowBox[{"-", "1"}], "+", "hsero", "+", "\[Mu]", "-", 
         RowBox[{"hsero", " ", "\[Mu]"}]}], ")"}]}]}]]}]}], "}"}]], "Output",E\
xpressionUUID->"dd72a1c9-1b24-4a28-9805-e12ae0baf4be"]
}, Open  ]],

Cell["\<\
Rewrite the solution for [b], defined in [bsolve], with the force of \
infection terms F replaced with the expression above\
\>", "Text",ExpressionUUID->"290ebe75-75a4-4d3f-8d73-abbe35b16634"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"bsolve", "/.", "foi"}], "]"}]}]], "Input",ExpressionUUID->\
"fdb2205c-59df-4a0a-8315-55c4ff532101"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"{", 
   RowBox[{"b", "\[Rule]", 
    RowBox[{"-", 
     FractionBox[
      RowBox[{"count", " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", 
         RowBox[{"a", " ", "hsero", " ", 
          RowBox[{"(", 
           RowBox[{"lam", "-", "\[Lambda]"}], ")"}], " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}], "+", 
         RowBox[{"d", " ", 
          RowBox[{"(", 
           RowBox[{
            RowBox[{"-", "1"}], "+", "hsero"}], ")"}], " ", "\[Mu]"}], "+", 
         RowBox[{"hsero", " ", "\[Lambda]", " ", "\[Mu]"}]}], ")"}]}], 
      RowBox[{
       RowBox[{"-", "1"}], "+", "\[Mu]"}]]}]}], "}"}], "}"}]], "Output",Expres\
sionUUID->"005a9ede-0a6b-46b3-a5e8-8d0ab188305e"]
}, Open  ]],

Cell["\<\
The rate of new cases of Lassa is F*s. Calculate this, and rewrite everything \
in terms of seroprevalence. The decrease in s at equilibrium that is caused \
by a is exactly cancelled by the increase in force of infection that is \
caused by a. \
\>", "Text",ExpressionUUID->"256d7114-002d-4029-b566-c3e82ac8bf06"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"NewCaseRate", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{"F", "/.", "foi"}], ")"}], "*", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"s", "/.", "eq"}], "/.", "foi"}], ")"}]}], "/.", "bsolve"}], 
    "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"13aca1fb-b81e-434f-\
98fa-90bb0b3e99d9"],

Cell[BoxData[
 RowBox[{"{", 
  RowBox[{"-", 
   FractionBox[
    RowBox[{"count", " ", "hsero", " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "\[Gamma]"}], ")"}], " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
    RowBox[{"\[Gamma]", " ", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"-", "1"}], "+", "\[Mu]"}], ")"}]}]]}], "}"}]], "Output",Expres\
sionUUID->"f940c71d-2bb9-4d0a-a322-d2de4561ccd9"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1440, 780},
WindowMargins->{{219, Automatic}, {Automatic, 118}},
Visible->True,
ScrollingOptions->{"VerticalScrollRange"->Fit},
ShowCellBracket->Automatic,
CellContext->Notebook,
TrackCellChangeTimes->False,
SpellingDictionaries->{"CorrectWords"->{"seroprevalence"}},
FrontEndVersion->"11.2 for Linux x86 (64-bit) (September 10, 2017)",
StyleDefinitions->"Default.nb"
]
(* End of Notebook Content *)

(* Internal cache information *)
(*CellTagsOutline
CellTagsIndex->{}
*)
(*CellTagsIndex
CellTagsIndex->{}
*)
(*NotebookFileOutline
Notebook[{
Cell[CellGroupData[{
Cell[1510, 35, 93, 0, 98, "Title",ExpressionUUID->"e8bd921e-cf66-4ebd-9513-78eb872ef4e0"],
Cell[1606, 37, 225, 3, 35, "Text",ExpressionUUID->"62b26a79-91fc-4292-bade-fd7745f6e7b1"],
Cell[CellGroupData[{
Cell[1856, 44, 89, 0, 65, "Subchapter",ExpressionUUID->"bf5dc25c-c5a4-43fc-9614-1e3c3a1a6627"],
Cell[1948, 46, 128, 0, 35, "Text",ExpressionUUID->"afbdb8fc-6b8e-47bd-a8cb-dbd223d92b9a"],
Cell[2079, 48, 831, 24, 78, "Input",ExpressionUUID->"c2034f9e-36c7-41d2-b8d0-3beb84a8be0b"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2947, 77, 198, 5, 65, "Subchapter",ExpressionUUID->"4b2b9d26-da4c-4491-80b3-de7bbda8e7fc"],
Cell[3148, 84, 92, 0, 35, "Text",ExpressionUUID->"0e34a607-f81a-4810-ad25-61c41f260938"],
Cell[CellGroupData[{
Cell[3265, 88, 469, 12, 31, "Input",ExpressionUUID->"93f7df5f-3eff-4d29-bb6c-999fff2366be"],
Cell[3737, 102, 1988, 60, 57, "Output",ExpressionUUID->"a7b364aa-15c1-442e-b579-621e4b233fd8"]
}, Open  ]],
Cell[5740, 165, 263, 4, 35, "Text",ExpressionUUID->"c2ea8db9-eed4-4832-8d72-a16e3823ce67"],
Cell[CellGroupData[{
Cell[6028, 173, 274, 7, 31, "Input",ExpressionUUID->"32742f3e-8886-422f-a2a1-20784ecfef8c"],
Cell[6305, 182, 917, 26, 62, "Output",ExpressionUUID->"8afa6e5c-1107-40ea-8bcb-dbd1bc9d975a"]
}, Open  ]],
Cell[7237, 211, 224, 3, 35, "Text",ExpressionUUID->"0632205d-9488-4c23-8efd-ab33142bad23"],
Cell[CellGroupData[{
Cell[7486, 218, 220, 5, 31, "Input",ExpressionUUID->"f26625be-e98c-466f-874e-df94195b5e1e"],
Cell[7709, 225, 1146, 29, 62, "Output",ExpressionUUID->"d4a99042-dbfe-4563-a46b-ae24e1f002b4"]
}, Open  ]],
Cell[8870, 257, 126, 0, 35, "Text",ExpressionUUID->"6b092b83-ee5d-45b4-bedc-f12e405ef8cd"],
Cell[CellGroupData[{
Cell[9021, 261, 276, 8, 31, "Input",ExpressionUUID->"2bad1404-bbfe-4377-b9c3-9be5d4530b36"],
Cell[9300, 271, 634, 18, 57, "Output",ExpressionUUID->"ce9910af-960e-4069-8417-dcc388140caf"]
}, Open  ]],
Cell[9949, 292, 143, 2, 35, "Text",ExpressionUUID->"7d9d9192-4ed3-456f-9d77-dde61685e10a"],
Cell[CellGroupData[{
Cell[10117, 298, 311, 8, 31, "Input",ExpressionUUID->"9316df9a-cd22-498f-b8b6-6435e808a1e8"],
Cell[10431, 308, 631, 18, 57, "Output",ExpressionUUID->"07582294-4d3f-4d5e-9212-829e86d47d28"]
}, Open  ]],
Cell[11077, 329, 202, 3, 35, "Text",ExpressionUUID->"7eac6851-b102-447b-a630-3bbd1ad5a950"],
Cell[CellGroupData[{
Cell[11304, 336, 192, 4, 31, "Input",ExpressionUUID->"2dbd6c6f-97ed-4841-bd4f-2a4511237572"],
Cell[11499, 342, 551, 16, 56, "Output",ExpressionUUID->"1ffee6b8-8aed-4a56-9280-bfacbb782b38"]
}, Open  ]],
Cell[12065, 361, 186, 3, 35, "Text",ExpressionUUID->"a2cbf31a-22dc-4e2e-91e4-a8aa58df158c"],
Cell[CellGroupData[{
Cell[12276, 368, 405, 12, 31, "Input",ExpressionUUID->"e344f092-8501-4b81-886a-b595b9388254"],
Cell[12684, 382, 433, 13, 57, "Output",ExpressionUUID->"f3f64f48-e396-4b0b-b0b1-79ef760b751d"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[13166, 401, 134, 0, 65, "Subchapter",ExpressionUUID->"2572d4f9-224b-49c7-b5e6-86a655b367f8"],
Cell[13303, 403, 141, 2, 35, "Text",ExpressionUUID->"be57a033-1922-4a08-bbdc-0de7a05bab83"],
Cell[CellGroupData[{
Cell[13469, 409, 422, 12, 31, "Input",ExpressionUUID->"ddbc0717-8bfe-4f89-8870-e1b0df342e55"],
Cell[13894, 423, 156, 4, 55, "Output",ExpressionUUID->"95652392-fc06-422e-929b-01cb3848186c"]
}, Open  ]],
Cell[14065, 430, 424, 12, 35, "Text",ExpressionUUID->"56599fca-1177-4675-bbac-69686eebea23"],
Cell[CellGroupData[{
Cell[14514, 446, 148, 3, 31, "Input",ExpressionUUID->"44a34a14-851a-42c1-9b40-2a97534f01d4"],
Cell[14665, 451, 123, 1, 35, "Output",ExpressionUUID->"b1c6d1b6-a6a2-465b-9000-1846f83f324b"]
}, Open  ]],
Cell[14803, 455, 197, 3, 35, "Text",ExpressionUUID->"4ffb9d56-9187-42aa-b8c2-38fe3f9b021b"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15037, 463, 129, 0, 65, "Subchapter",ExpressionUUID->"20a2bd3f-2233-4454-a614-34d6bd8e0096"],
Cell[15169, 465, 143, 2, 35, "Text",ExpressionUUID->"0a8428b2-4919-445f-aa88-bda270e09c68"],
Cell[CellGroupData[{
Cell[15337, 471, 379, 11, 31, "Input",ExpressionUUID->"4a900637-e7e7-4298-a70d-bfd610027a57"],
Cell[15719, 484, 160, 4, 54, "Output",ExpressionUUID->"95288e70-f4a3-4af8-a3f4-65ba976b910c"]
}, Open  ]],
Cell[15894, 491, 440, 9, 42, "Text",ExpressionUUID->"b7f09cc0-7381-4115-a932-bc9ce3065a41"],
Cell[16337, 502, 190, 5, 31, "Input",ExpressionUUID->"9e0c04ef-017e-4376-b057-25f1dd3b5d6c"],
Cell[16530, 509, 182, 3, 35, "Text",ExpressionUUID->"fc10009b-48de-4c99-a925-f8bcb4e9c63b"],
Cell[CellGroupData[{
Cell[16737, 516, 338, 9, 31, "Input",ExpressionUUID->"3f93b19e-c657-4271-acc9-643d36b44d33"],
Cell[17078, 527, 108, 1, 35, "Output",ExpressionUUID->"ff1522a4-b20d-4629-aa94-3373ac8c0311"]
}, Open  ]],
Cell[17201, 531, 201, 3, 35, "Text",ExpressionUUID->"0e3c9cd6-dd3c-4a3b-b512-c9497faa5d20"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[17451, 540, 128, 0, 98, "Title",ExpressionUUID->"1420339a-6f1f-4312-b2d3-5bb90cd4c53d"],
Cell[17582, 542, 376, 5, 58, "Text",ExpressionUUID->"4dbfbd23-a124-414b-91f2-e290b5a4bbe5"],
Cell[CellGroupData[{
Cell[17983, 551, 89, 0, 65, "Subchapter",ExpressionUUID->"632f3d9d-12ae-4ca9-9a62-f069668a78e2"],
Cell[18075, 553, 190, 3, 35, "Text",ExpressionUUID->"c3f12661-cdd3-4d43-9101-cf97ceb7b28d"],
Cell[18268, 558, 1084, 32, 101, "Input",ExpressionUUID->"dc2fd5cd-5006-4506-9277-f4fb4eb4b0a2"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19389, 595, 198, 5, 65, "Subchapter",ExpressionUUID->"4d9ab01c-3488-43d4-a66c-a533d09ca690"],
Cell[19590, 602, 240, 4, 35, "Text",ExpressionUUID->"37f4c21b-25e7-4fea-9d04-208fcba24398"],
Cell[CellGroupData[{
Cell[19855, 610, 539, 13, 31, "Input",ExpressionUUID->"e4dafdce-a95a-4d56-9cac-056895578223"],
Cell[20397, 625, 3354, 88, 109, "Output",ExpressionUUID->"fb9b51c4-37b7-4b8f-9f33-b039b1cfd462"]
}, Open  ]],
Cell[23766, 716, 263, 4, 35, "Text"],
Cell[CellGroupData[{
Cell[24054, 724, 301, 8, 31, "Input"],
Cell[24358, 734, 1354, 36, 66, "Output"]
}, Open  ]],
Cell[25727, 773, 224, 3, 35, "Text"],
Cell[CellGroupData[{
Cell[25976, 780, 220, 5, 31, "Input"],
Cell[26199, 787, 1603, 38, 62, "Output"]
}, Open  ]],
Cell[27817, 828, 126, 0, 35, "Text"],
Cell[CellGroupData[{
Cell[27968, 832, 276, 8, 31, "Input"],
Cell[28247, 842, 850, 23, 58, "Output"]
}, Open  ]],
Cell[29112, 868, 220, 3, 35, "Text"],
Cell[CellGroupData[{
Cell[29357, 875, 311, 8, 31, "Input"],
Cell[29671, 885, 856, 23, 58, "Output"]
}, Open  ]],
Cell[30542, 911, 202, 3, 35, "Text"],
Cell[CellGroupData[{
Cell[30769, 918, 192, 4, 31, "Input"],
Cell[30964, 924, 783, 22, 56, "Output"]
}, Open  ]],
Cell[31762, 949, 324, 5, 58, "Text"],
Cell[CellGroupData[{
Cell[32111, 958, 405, 12, 31, "Input"],
Cell[32519, 972, 433, 13, 57, "Output"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature AwD8G50i6RKaLBKZ6bDGGRFg *)
