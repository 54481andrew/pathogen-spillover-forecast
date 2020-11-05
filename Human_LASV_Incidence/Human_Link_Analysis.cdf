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
NotebookDataLength[     21626,        620]
NotebookOptionsPosition[     17430,        538]
NotebookOutlinePosition[     17939,        559]
CellTagsIndexPosition[     17896,        556]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Human Lassa Case Rate", "Title",ExpressionUUID->"b1fc522e-426d-48e3-8830-5065e1a46473"],

Cell["\<\
This code shows the analyses performed on the system of ordinary differential \
equations that describe the spillover of Lassa virus into humans. \
\>", "Text",ExpressionUUID->"4463146e-3fa4-495d-9069-1e27af69ddd3"],

Cell[CellGroupData[{

Cell["Define model", "Subchapter",ExpressionUUID->"2850ec8d-850b-4c7b-bbf4-a54ab4a2573d"],

Cell["Define the system of ODEs that describe human infection. ", "Text",ExpressionUUID->"42116ebe-3b81-4f22-aef9-2e8449d13856"],

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
"501b1d5c-bb4f-4508-b7d9-fdc8c6b6c12f"]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[{
 "Derive relationship between ",
 StyleBox["seroprevalence",
  IgnoreSpellCheck->True],
 " and new case rate"
}], "Subchapter",ExpressionUUID->"f64bab01-404a-47bc-8b67-10e9a5498c16"],

Cell["Find the equilibrium.", "Text",ExpressionUUID->"cb82cdd6-f66f-4f89-9951-0201366ec56a"],

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
ressionUUID->"70ae355e-aa2e-4145-90b4-e0ac2dc958c9"],

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
   "}"}], "}"}]], "Output",ExpressionUUID->"82b897dd-0c65-4a3a-87f5-\
0d4903f2b12a"]
}, Open  ]],

Cell["\<\
Find seroprevalence at equilibrium (this is the data we have to work with). \
Note that this requires calculating the steady state population size in the \
presence of LASV virulence. \
\>", "Text",ExpressionUUID->"cb7f7ef7-b7e4-4065-9403-7301132cffe9"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Pop", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"s", " ", "+", " ", "i", " ", "+", " ", "r"}], " ", ")"}], "/.", 
    "eq"}], "]"}]}]], "Input",ExpressionUUID->"d7a3ff5c-feef-4075-9b03-\
d46fc7466e20"],

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
  "}"}]], "Output",ExpressionUUID->"3251a29b-8488-4c3f-9775-65dc9fc2076e"]
}, Open  ]],

Cell["\<\
Variable [count] refers to the population size in a pixel of the population \
raster.  Solve for the birth rate [b] for a given value of [count].\
\>", "Text",ExpressionUUID->"cdc6993f-7004-47d1-ace1-7f2c6e9b097e"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"Pop", "\[Equal]", "count"}], ",", " ", "b"}], "]"}]}]], "Input",E\
xpressionUUID->"96c90aa5-14a9-4351-8242-e22420208945"],

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
  "}"}]], "Output",ExpressionUUID->"dd9a85f9-6e2e-4093-a22d-096190f1f033"]
}, Open  ]],

Cell["Calculate seroprevalence in terms of model parameters. ", "Text",ExpressionUUID->"ac79dd3a-3f4a-47e5-b8e7-d728e08ef146"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Rstar", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"r", "/", 
      RowBox[{"(", "Pop", ")"}]}], "/.", "eq"}], "]"}], "//", 
   "Flatten"}]}]], "Input",ExpressionUUID->"9c2c00dc-fa32-4286-80a1-\
3263a2103cf0"],

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
xpressionUUID->"73988002-d076-4c63-9523-d60c1ed96b77"]
}, Open  ]],

Cell["\<\
Solve for the force of infection, F, in terms of seroprevalence.\
\>", "Text",ExpressionUUID->"5e14727e-3044-4ed0-9b89-854712b174be"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"foi", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{"Solve", "[", 
     RowBox[{
      RowBox[{"hsero", " ", "\[Equal]", " ", "Rstar"}], ",", " ", "F"}], 
     "]"}], "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"1e0ed308-\
b0b6-4ee1-b54e-37cd098d68c1"],

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
xpressionUUID->"503b3566-5bdd-4fe4-a781-1047d6b5a4dc"]
}, Open  ]],

Cell["\<\
Rewrite the solution for [b], defined in [bsolve], with the force of \
infection terms F replaced with the expression above\
\>", "Text",ExpressionUUID->"499923f6-113b-40bc-88aa-0c9cc16f1f9a"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"bsolve", "/.", "foi"}], "]"}]}]], "Input",ExpressionUUID->\
"603288de-07f2-4679-b68c-186bc6c64763"],

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
sionUUID->"17a361df-e515-4901-bcf4-3809cc6e50c9"]
}, Open  ]],

Cell["\<\
The rate of new cases of Lassa is F*s. Calculate this, and rewrite everything \
in terms of seroprevalence.\
\>", "Text",ExpressionUUID->"2048ec8e-decb-49ab-8284-c3679ac251b0"],

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
    "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"5bbffc46-64b3-4d5f-\
b37b-d551b9a41390"],

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
sionUUID->"c31edf2b-1deb-4549-bd7a-96a7e31db74c"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Compare case rate estimates for virulence vs no virulence", "Subchapter",ExpressionUUID->"fbe98d2e-94ef-4721-a999-7930bb12f5d4"],

Cell["\<\
Calculate the ratio of cases with, versus without, virulence. \
\>", "Text",ExpressionUUID->"aafc5193-131c-4d88-924a-283a0c66c519"],

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
   "Flatten"}]}]], "Input",ExpressionUUID->"57d90329-d407-4cb4-8c83-\
b119fc4b0cfb"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox["1", 
   RowBox[{"1", "-", "\[Mu]"}]], "}"}]], "Output",ExpressionUUID->"e860b986-\
a46f-49e1-b6e7-7acc377f2028"]
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
  "f78cb987-241b-40b1-9653-f47665d81117"],
 ".  "
}], "Text",ExpressionUUID->"7343e1a3-a050-4aac-ada6-b28e012778e8"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"\[Rho]", "/.", 
  RowBox[{"\[Mu]", "\[Rule]", "0.02"}]}]], "Input",ExpressionUUID->"eb9ecf19-\
58bb-424d-9a76-a55dd6caa242"],

Cell[BoxData[
 RowBox[{"{", "1.0204081632653061`", "}"}]], "Output",ExpressionUUID->"a86f1624-19d4-4625-b30b-4e0b6ad0f016"]
}, Open  ]],

Cell["\<\
When \[Mu] = 0.02 (lower value from McCormick, Webb, 1987), increase is \
~1.02, meaning roughly 2% more cases occur. \
\>", "Text",ExpressionUUID->"548300cb-e584-49cf-9b95-945e2e674712"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Compare estimates for reinfection vs no reinfection ", "Subchapter",ExpressionUUID->"33c08662-df47-4bee-afe9-d6c05a40f8e7"],

Cell["\<\
Calculate the ratio of cases with, versus without, reinfection. \
\>", "Text",ExpressionUUID->"b345b890-d016-4a19-9716-0b48ba169a29"],

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
  "Flatten"}]], "Input",ExpressionUUID->"38826dc2-a6b0-4235-8609-\
dd8326fb6d11"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox[
   RowBox[{"d", "+", "\[Lambda]"}], "d"], "}"}]], "Output",ExpressionUUID->\
"15e82ce5-0c7f-4984-bba9-00449d021c7f"]
}, Open  ]],

Cell[TextData[{
 "Reinfection multiplies the case rate by ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"d", "+", "\[Lambda]"}], "d"], TraditionalForm]],ExpressionUUID->
  "91017dfc-fab7-4def-8873-24890f2fcd5e"],
 ". Calculate the effect of reinfection for a mean lifespan of 50 years, and \
reinfection rate \[Lambda] = 0.064 / yr (from McCormick, Webb, 1987)."
}], "Text",ExpressionUUID->"74989c06-28fa-4945-af15-4079d2da773b"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Pars", "=", 
   RowBox[{"{", 
    RowBox[{"d", "\[Rule]", " ", "0.02"}], "}"}]}], ";"}]], "Input",Expression\
UUID->"9f6e1e35-20ea-450c-b70f-0a73f6d8be2f"],

Cell["\<\
Ratio of estimates of new cases with (\[Lambda] = 0.064) and without (\
\[Lambda] = 0) seroreversion is\
\>", "Text",ExpressionUUID->"4cd517e5-b200-4095-8e2f-a7253696ab83"],

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
  "Pars"}]], "Input",ExpressionUUID->"28eac22c-8577-412a-808d-ae406f46bb35"],

Cell[BoxData[
 RowBox[{"{", "4.2`", "}"}]], "Output",ExpressionUUID->"a289e475-f8a1-435b-8311-4829de7eec70"]
}, Open  ]],

Cell["\<\
Including reinfection multiplies the case estimates by a factor of 420%; \
including virulence increases estimates by 2%. \
\>", "Text",ExpressionUUID->"ac93e7b2-fb5a-4397-b90a-b99444692a5a"]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1920, 1052},
Visible->True,
ScrollingOptions->{"VerticalScrollRange"->Fit},
ShowCellBracket->Automatic,
Deployed->True,
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
Cell[1510, 35, 93, 0, 98, "Title",ExpressionUUID->"b1fc522e-426d-48e3-8830-5065e1a46473"],
Cell[1606, 37, 225, 3, 35, "Text",ExpressionUUID->"4463146e-3fa4-495d-9069-1e27af69ddd3"],
Cell[CellGroupData[{
Cell[1856, 44, 89, 0, 65, "Subchapter",ExpressionUUID->"2850ec8d-850b-4c7b-bbf4-a54ab4a2573d"],
Cell[1948, 46, 128, 0, 35, "Text",ExpressionUUID->"42116ebe-3b81-4f22-aef9-2e8449d13856"],
Cell[2079, 48, 831, 24, 78, "Input",ExpressionUUID->"501b1d5c-bb4f-4508-b7d9-fdc8c6b6c12f"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2947, 77, 198, 5, 65, "Subchapter",ExpressionUUID->"f64bab01-404a-47bc-8b67-10e9a5498c16"],
Cell[3148, 84, 92, 0, 35, "Text",ExpressionUUID->"cb82cdd6-f66f-4f89-9951-0201366ec56a"],
Cell[CellGroupData[{
Cell[3265, 88, 469, 12, 31, "Input",ExpressionUUID->"70ae355e-aa2e-4145-90b4-e0ac2dc958c9"],
Cell[3737, 102, 1988, 60, 57, "Output",ExpressionUUID->"82b897dd-0c65-4a3a-87f5-0d4903f2b12a"]
}, Open  ]],
Cell[5740, 165, 263, 4, 35, "Text",ExpressionUUID->"cb7f7ef7-b7e4-4065-9403-7301132cffe9"],
Cell[CellGroupData[{
Cell[6028, 173, 274, 7, 31, "Input",ExpressionUUID->"d7a3ff5c-feef-4075-9b03-d46fc7466e20"],
Cell[6305, 182, 917, 26, 62, "Output",ExpressionUUID->"3251a29b-8488-4c3f-9775-65dc9fc2076e"]
}, Open  ]],
Cell[7237, 211, 224, 3, 35, "Text",ExpressionUUID->"cdc6993f-7004-47d1-ace1-7f2c6e9b097e"],
Cell[CellGroupData[{
Cell[7486, 218, 220, 5, 31, "Input",ExpressionUUID->"96c90aa5-14a9-4351-8242-e22420208945"],
Cell[7709, 225, 1146, 29, 62, "Output",ExpressionUUID->"dd9a85f9-6e2e-4093-a22d-096190f1f033"]
}, Open  ]],
Cell[8870, 257, 126, 0, 35, "Text",ExpressionUUID->"ac79dd3a-3f4a-47e5-b8e7-d728e08ef146"],
Cell[CellGroupData[{
Cell[9021, 261, 276, 8, 31, "Input",ExpressionUUID->"9c2c00dc-fa32-4286-80a1-3263a2103cf0"],
Cell[9300, 271, 634, 18, 57, "Output",ExpressionUUID->"73988002-d076-4c63-9523-d60c1ed96b77"]
}, Open  ]],
Cell[9949, 292, 143, 2, 35, "Text",ExpressionUUID->"5e14727e-3044-4ed0-9b89-854712b174be"],
Cell[CellGroupData[{
Cell[10117, 298, 311, 8, 31, "Input",ExpressionUUID->"1e0ed308-b0b6-4ee1-b54e-37cd098d68c1"],
Cell[10431, 308, 631, 18, 57, "Output",ExpressionUUID->"503b3566-5bdd-4fe4-a781-1047d6b5a4dc"]
}, Open  ]],
Cell[11077, 329, 202, 3, 35, "Text",ExpressionUUID->"499923f6-113b-40bc-88aa-0c9cc16f1f9a"],
Cell[CellGroupData[{
Cell[11304, 336, 192, 4, 31, "Input",ExpressionUUID->"603288de-07f2-4679-b68c-186bc6c64763"],
Cell[11499, 342, 551, 16, 56, "Output",ExpressionUUID->"17a361df-e515-4901-bcf4-3809cc6e50c9"]
}, Open  ]],
Cell[12065, 361, 186, 3, 35, "Text",ExpressionUUID->"2048ec8e-decb-49ab-8284-c3679ac251b0"],
Cell[CellGroupData[{
Cell[12276, 368, 405, 12, 31, "Input",ExpressionUUID->"5bbffc46-64b3-4d5f-b37b-d551b9a41390"],
Cell[12684, 382, 433, 13, 57, "Output",ExpressionUUID->"c31edf2b-1deb-4549-bd7a-96a7e31db74c"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[13166, 401, 134, 0, 65, "Subchapter",ExpressionUUID->"fbe98d2e-94ef-4721-a999-7930bb12f5d4"],
Cell[13303, 403, 141, 2, 35, "Text",ExpressionUUID->"aafc5193-131c-4d88-924a-283a0c66c519"],
Cell[CellGroupData[{
Cell[13469, 409, 422, 12, 31, "Input",ExpressionUUID->"57d90329-d407-4cb4-8c83-b119fc4b0cfb"],
Cell[13894, 423, 156, 4, 55, "Output",ExpressionUUID->"e860b986-a46f-49e1-b6e7-7acc377f2028"]
}, Open  ]],
Cell[14065, 430, 424, 12, 35, "Text",ExpressionUUID->"7343e1a3-a050-4aac-ada6-b28e012778e8"],
Cell[CellGroupData[{
Cell[14514, 446, 148, 3, 31, "Input",ExpressionUUID->"eb9ecf19-58bb-424d-9a76-a55dd6caa242"],
Cell[14665, 451, 123, 1, 35, "Output",ExpressionUUID->"a86f1624-19d4-4625-b30b-4e0b6ad0f016"]
}, Open  ]],
Cell[14803, 455, 197, 3, 35, "Text",ExpressionUUID->"548300cb-e584-49cf-9b95-945e2e674712"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15037, 463, 129, 0, 65, "Subchapter",ExpressionUUID->"33c08662-df47-4bee-afe9-d6c05a40f8e7"],
Cell[15169, 465, 143, 2, 35, "Text",ExpressionUUID->"b345b890-d016-4a19-9716-0b48ba169a29"],
Cell[CellGroupData[{
Cell[15337, 471, 379, 11, 31, "Input",ExpressionUUID->"38826dc2-a6b0-4235-8609-dd8326fb6d11"],
Cell[15719, 484, 160, 4, 54, "Output",ExpressionUUID->"15e82ce5-0c7f-4984-bba9-00449d021c7f"]
}, Open  ]],
Cell[15894, 491, 440, 9, 42, "Text",ExpressionUUID->"74989c06-28fa-4945-af15-4079d2da773b"],
Cell[16337, 502, 190, 5, 31, "Input",ExpressionUUID->"9f6e1e35-20ea-450c-b70f-0a73f6d8be2f"],
Cell[16530, 509, 182, 3, 35, "Text",ExpressionUUID->"4cd517e5-b200-4095-8e2f-a7253696ab83"],
Cell[CellGroupData[{
Cell[16737, 516, 338, 9, 31, "Input",ExpressionUUID->"28eac22c-8577-412a-808d-ae406f46bb35"],
Cell[17078, 527, 108, 1, 35, "Output",ExpressionUUID->"a289e475-f8a1-435b-8311-4829de7eec70"]
}, Open  ]],
Cell[17201, 531, 201, 3, 35, "Text",ExpressionUUID->"ac93e7b2-fb5a-4397-b90a-b99444692a5a"]
}, Open  ]]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature SxpzvyyXsmR97DwHiVeUDkKa *)
