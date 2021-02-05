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
NotebookDataLength[     40173,       1121]
NotebookOptionsPosition[     33050,        992]
NotebookOutlinePosition[     33595,       1013]
CellTagsIndexPosition[     33552,       1010]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Human Lassa spillover rate", "Title",ExpressionUUID->"5298ff9f-40e3-4ff3-9ff5-2440ed9a9cf1"],

Cell["\<\
This code shows the analyses performed on the system of ordinary differential \
equations that describe the spillover of Lassa virus into humans. \
\>", "Text",ExpressionUUID->"24b28408-ec66-4eae-b220-c852c93b1e71"],

Cell[CellGroupData[{

Cell["Define model", "Subchapter",ExpressionUUID->"8755134a-6eac-4bc9-9856-401e29b9e064"],

Cell["Define the system of ODEs that describe human infection. ", "Text",ExpressionUUID->"056026c3-8a33-49c7-b3b7-a68cbc0ac675"],

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
"dc9fe8b2-9529-46ab-aae2-06f339f999e9"]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[{
 "Derive relationship between ",
 StyleBox["seroprevalence",
  IgnoreSpellCheck->True],
 " and spillover rate"
}], "Subchapter",ExpressionUUID->"58bb189a-cc09-4d6a-a254-5b13ccf3bed0"],

Cell["Find the equilibrium.", "Text",ExpressionUUID->"d9462527-7939-40cf-ba01-7789f1fd5955"],

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
ressionUUID->"6b827d1f-61a7-48b4-acee-01acb3d3954e"],

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
   "}"}], "}"}]], "Output",ExpressionUUID->"94aeb62a-5a2a-4464-9b2f-\
1ca2eb1eb8e8"]
}, Open  ]],

Cell["\<\
Find seroprevalence at equilibrium (this is the data we have to work with). \
Note that this requires calculating the steady state population size in the \
presence of LASV virulence. \
\>", "Text",ExpressionUUID->"94f5f41b-d930-4e06-8c74-83953ea28a46"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Pop", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"s", " ", "+", " ", "i", " ", "+", " ", "r"}], " ", ")"}], "/.", 
    "eq"}], "]"}]}]], "Input",ExpressionUUID->"c12493b6-9170-4608-8801-\
8cb0460c139a"],

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
  "}"}]], "Output",ExpressionUUID->"2603acf5-ab57-44c5-bd37-c99dcfd7e528"]
}, Open  ]],

Cell["\<\
Variable [count] refers to the population size in a pixel of the population \
raster.  Solve for the birth rate [b] for a given value of [count].\
\>", "Text",ExpressionUUID->"d7d137d3-fb68-45c0-af29-0ef1092966c9"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"Pop", "\[Equal]", "count"}], ",", " ", "b"}], "]"}]}]], "Input",E\
xpressionUUID->"4125e9a8-b744-4712-827e-a2d34cef7ee8"],

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
  "}"}]], "Output",ExpressionUUID->"c2849a23-df3d-40e8-ba94-e2f1905a2bb3"]
}, Open  ]],

Cell["Calculate seroprevalence in terms of model parameters. ", "Text",ExpressionUUID->"2ebd94f1-c81d-4e6a-9f35-d2dcfaab10dd"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Rstar", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"r", "/", 
      RowBox[{"(", "Pop", ")"}]}], "/.", "eq"}], "]"}], "//", 
   "Flatten"}]}]], "Input",ExpressionUUID->"ab7fa8e1-43c0-4313-a12e-\
7fef26d8df28"],

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
xpressionUUID->"e734f0d3-edff-4e5c-b679-b1894b2acf5d"]
}, Open  ]],

Cell["\<\
Solve for the force of infection, F, in terms of seroprevalence.\
\>", "Text",ExpressionUUID->"c6115c4f-bb63-4a3c-b860-486cb0267a60"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"foi", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{"Solve", "[", 
     RowBox[{
      RowBox[{"hsero", " ", "\[Equal]", " ", "Rstar"}], ",", " ", "F"}], 
     "]"}], "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"55ad8c49-\
881b-441c-bf61-6f7bd56c9d5c"],

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
xpressionUUID->"b9052ca2-35f6-4b36-98f9-496a13888e33"]
}, Open  ]],

Cell["\<\
Rewrite the solution for [b], defined in [bsolve], with the force of \
infection terms F replaced with the expression above\
\>", "Text",ExpressionUUID->"43561019-7baf-4587-bfa1-2b47292636cd"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"bsolve", "/.", "foi"}], "]"}]}]], "Input",ExpressionUUID->\
"2faf48f4-6027-45bf-ae27-6f92e623d7e7"],

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
sionUUID->"e9388348-3989-43dd-9b08-f32ced23ddf5"]
}, Open  ]],

Cell["\<\
The rate of new infections of Lassa is F*s. Calculate this, and rewrite \
everything in terms of seroprevalence.\
\>", "Text",ExpressionUUID->"797fd1d9-78e5-4e83-8445-62a0c5784161"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"NewInfRate", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{"F", "/.", "foi"}], ")"}], "*", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"s", "/.", "eq"}], "/.", "foi"}], ")"}]}], "/.", "bsolve"}], 
    "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"67dd637d-88ee-48d2-\
894f-827b68820c69"],

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
sionUUID->"11a7d8ef-9f41-4ee8-8c8c-c345411e019f"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["\<\
Compare spillover rate estimates for virulence vs no virulence\
\>", "Subchapter",ExpressionUUID->"3f97b281-6ea6-4e21-8806-b978700e2f3b"],

Cell["\<\
Calculate the ratio of infections with, versus without, virulence. \
\>", "Text",ExpressionUUID->"0ed46c9a-8ac6-459e-8c9a-9fc91b86c2f2"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"\[Rho]", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"(", 
      RowBox[{"NewInfRate", "/.", "eq"}], ")"}], "/", 
     RowBox[{"(", 
      RowBox[{
       RowBox[{"NewInfRate", "/.", "eq"}], "/.", 
       RowBox[{"\[Mu]", "\[Rule]", "0"}]}], ")"}]}], " ", "]"}], "//", 
   "Flatten"}]}]], "Input",ExpressionUUID->"96fbff50-f1ee-4815-af19-\
eb7ad1af2fac"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox["1", 
   RowBox[{"1", "-", "\[Mu]"}]], "}"}]], "Output",ExpressionUUID->"311a258c-\
a2a2-4c69-b5cb-b0747fb9a29b"]
}, Open  ]],

Cell[TextData[{
 "In words, including virulence increases the estimated rate of LASV \
infections  by ",
 Cell[BoxData[
  FormBox[
   RowBox[{"\[Rho]", " ", "=", " ", 
    SuperscriptBox[
     RowBox[{"(", 
      RowBox[{"1", "-", "\[Mu]"}], ")"}], 
     RowBox[{"-", "1"}]]}], TraditionalForm]],ExpressionUUID->
  "9e5ea74e-9487-4f64-8378-03e50bf6cb54"],
 ".  "
}], "Text",ExpressionUUID->"28b46dfe-2e26-4720-a085-be09b6746c44"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"\[Rho]", "/.", 
  RowBox[{"\[Mu]", "\[Rule]", "0.02"}]}]], "Input",ExpressionUUID->"f246dc43-\
d23a-4efc-be3f-0a6a279defd0"],

Cell[BoxData[
 RowBox[{"{", "1.0204081632653061`", "}"}]], "Output",ExpressionUUID->"4ccee5c8-1f8c-41df-9e00-184eed409b65"]
}, Open  ]],

Cell["\<\
When \[Mu] = 0.02 (lower value from McCormick, Webb, 1987), increase is \
~1.02, meaning roughly 2% more infections occur. \
\>", "Text",ExpressionUUID->"bbeccc32-74d8-4f32-8ede-fec4631408cd"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Compare estimates for reinfection vs no reinfection ", "Subchapter",ExpressionUUID->"552fee79-005c-4d05-8c9f-927c576e1219"],

Cell["\<\
Calculate the ratio of infections with, versus without, reinfection. \
\>", "Text",ExpressionUUID->"a2f8e18f-10e4-4548-bd66-c27d082b85dd"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"NewInfRate", "/.", "eq"}], ")"}], "/", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"NewInfRate", "/.", "eq"}], "/.", 
      RowBox[{"\[Lambda]", "\[Rule]", "0"}]}], ")"}]}], " ", "]"}], "//", 
  "Flatten"}]], "Input",ExpressionUUID->"60ae074c-51bf-4baa-bcfa-\
f918f09de36b"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox[
   RowBox[{"d", "+", "\[Lambda]"}], "d"], "}"}]], "Output",ExpressionUUID->\
"ad5c81a8-c807-4132-bc76-be31d0823541"]
}, Open  ]],

Cell[TextData[{
 "Reinfection multiplies the infection rate by ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"d", "+", "\[Lambda]"}], "d"], TraditionalForm]],ExpressionUUID->
  "e17836b2-dc3d-460c-a48b-6ed3a509f55a"],
 ". Calculate the effect of reinfection for a mean lifespan of 50 years, and \
reinfection rate \[Lambda] = 0.064 / yr (from McCormick, Webb, 1987)."
}], "Text",ExpressionUUID->"c6401ba6-f4b0-452e-8856-ebb257e30266"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Pars", "=", 
   RowBox[{"{", 
    RowBox[{"d", "\[Rule]", " ", "0.02"}], "}"}]}], ";"}]], "Input",Expression\
UUID->"0ff25392-1e02-4f1c-8d5f-870cac6b9dcb"],

Cell["\<\
Ratio of estimates of new infections with (\[Lambda] = 0.064) and without (\
\[Lambda] = 0) seroreversion is\
\>", "Text",ExpressionUUID->"7ff1c590-5417-4dba-9902-ca5c4d0a2483"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{
  RowBox[{
   RowBox[{"(", 
    RowBox[{"NewInfRate", "/.", 
     RowBox[{"\[Lambda]", "\[Rule]", "0.064"}]}], ")"}], "/", 
   RowBox[{"(", 
    RowBox[{"NewInfRate", "/.", 
     RowBox[{"\[Lambda]", "\[Rule]", "0"}]}], ")"}]}], "/.", 
  "Pars"}]], "Input",ExpressionUUID->"80f24613-750f-420b-9e12-d680f5e59793"],

Cell[BoxData[
 RowBox[{"{", "4.2`", "}"}]], "Output",ExpressionUUID->"860c6b87-c7fe-46f7-b0ec-8b585810a8d5"]
}, Open  ]],

Cell["\<\
Including reinfection multiplies the infection estimates by a factor of 420%; \
including virulence increases estimates by 2%. \
\>", "Text",ExpressionUUID->"909ae6c2-db01-4a10-8a85-db987ba0339f"]
}, Open  ]]
}, Open  ]],

Cell[CellGroupData[{

Cell["Human LASV spillover with undetectable antibody immune class", "Title",ExpressionUUID->"7d0f7f18-6e0a-4e2a-9066-03cff71fcb29"],

Cell["\<\
At the request of a reviewer, we also investigate the possibility that a \
fraction, a, of recovered individuals with immunity (R) transition into a \
class with undetectable immunity (C). This might be the case, for example, if \
immunity were caused primarily by a T-cell response in humans. \
\>", "Text",ExpressionUUID->"6f31b8c4-bded-42fa-a450-6d1a6300b25d"],

Cell[CellGroupData[{

Cell["Define model", "Subchapter",ExpressionUUID->"88338340-9a48-4924-8459-9f219097d12d"],

Cell["\<\
Define the system of ODEs that describe human infection.  Note the additional \
class, c, and the rate term a. \
\>", "Text",ExpressionUUID->"a364b61c-5d98-40eb-b5c2-dc15b3ecf97f"],

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
    RowBox[{"d", " ", "c"}]}]}], ";"}]}], "Input",ExpressionUUID->"616c8ac9-\
7672-4c52-ade4-3b23849df1e4"]
}, Open  ]],

Cell[CellGroupData[{

Cell[TextData[{
 "Derive relationship between ",
 StyleBox["seroprevalence",
  IgnoreSpellCheck->True],
 " and infection rate"
}], "Subchapter",ExpressionUUID->"affc0497-f7e3-4cee-90d9-5078aa44af58"],

Cell["\<\
The remainder of this analysis is identical to what was performed above. Note \
that the equilibrium values of s are smaller when a > 0. Find the \
equilibrium. \
\>", "Text",ExpressionUUID->"ad6b50bc-9bb4-4fc5-a33e-2b34bb14fa1d"],

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
   "]"}]}]], "Input",ExpressionUUID->"be05f869-f9a4-4cf9-8773-0ff3728639d6"],

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
  "}"}]], "Output",ExpressionUUID->"bff78f72-1712-46a1-a146-09b293eba06c"]
}, Open  ]],

Cell["\<\
Find seroprevalence at equilibrium (this is the data we have to work with). \
Note that this requires calculating the steady state population size in the \
presence of LASV virulence. \
\>", "Text",ExpressionUUID->"6e75533e-d8c3-4692-92e9-e8b9eace5799"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Pop", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{
     "s", " ", "+", " ", "i", " ", "+", " ", "r", " ", "+", " ", "c"}], " ", 
     ")"}], "/.", "eq"}], "]"}]}]], "Input",ExpressionUUID->"90a41a83-9199-\
4efe-843f-417e09b98fdf"],

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
ExpressionUUID->"f7ff57b5-93ba-4931-ab89-1ae08920899f"]
}, Open  ]],

Cell["\<\
Variable [count] refers to the population size in a pixel of the population \
raster.  Solve for the birth rate [b] for a given value of [count].\
\>", "Text",ExpressionUUID->"8013d08e-3c70-4155-ab65-150dbe167d38"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"Solve", "[", 
   RowBox[{
    RowBox[{"Pop", "\[Equal]", "count"}], ",", " ", "b"}], "]"}]}]], "Input",E\
xpressionUUID->"74071871-abfb-4360-8856-91dbb825893b"],

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
    "}"}], "}"}]], "Output",ExpressionUUID->"de91ae4f-26d7-4be3-b18b-\
8024ead62fd6"]
}, Open  ]],

Cell["Calculate seroprevalence in terms of model parameters. ", "Text",ExpressionUUID->"8e8fea37-3976-4850-84c0-18b42548008b"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Rstar", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{"r", "/", 
      RowBox[{"(", "Pop", ")"}]}], "/.", "eq"}], "]"}], "//", 
   "Flatten"}]}]], "Input",ExpressionUUID->"d9f7fe85-21e2-410e-b51f-\
ccb8d1c73291"],

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
xpressionUUID->"a7f6cd0a-d839-496d-8316-e6e0c6649831"]
}, Open  ]],

Cell["\<\
Solve for the force of infection, F, in terms of seroprevalence. Note that \
the equilibrium value of foi is made larger by the parameter a. \
\>", "Text",ExpressionUUID->"a1f51013-e7af-4658-b112-472e5f41e89e"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"foi", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{"Solve", "[", 
     RowBox[{
      RowBox[{"hsero", " ", "\[Equal]", " ", "Rstar"}], ",", " ", "F"}], 
     "]"}], "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"00e919c2-\
e7c5-48ff-aff3-6ede77e6b2b0"],

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
xpressionUUID->"5d064e5b-11f0-4474-b39a-23355ebd63e6"]
}, Open  ]],

Cell["\<\
Rewrite the solution for [b], defined in [bsolve], with the force of \
infection terms F replaced with the expression above\
\>", "Text",ExpressionUUID->"79c38174-bf34-4b23-8e84-b9796aeaa4aa"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"bsolve", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{"bsolve", "/.", "foi"}], "]"}]}]], "Input",ExpressionUUID->\
"04633174-1cc2-4203-a690-b6013776d0c1"],

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
sionUUID->"d6625080-db41-4e35-b925-4d8d16785147"]
}, Open  ]],

Cell["\<\
The rate of new infections of Lassa is F*s. Calculate this, and rewrite \
everything in terms of seroprevalence. The decrease in s at equilibrium that \
is caused by a is exactly cancelled by the increase in force of infection \
that is caused by a. \
\>", "Text",ExpressionUUID->"cd9ce0e8-dd96-4f76-bad7-c82aef6f32da"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"NewInfRate", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{
     RowBox[{
      RowBox[{"(", 
       RowBox[{"F", "/.", "foi"}], ")"}], "*", 
      RowBox[{"(", 
       RowBox[{
        RowBox[{"s", "/.", "eq"}], "/.", "foi"}], ")"}]}], "/.", "bsolve"}], 
    "]"}], "//", "Flatten"}]}]], "Input",ExpressionUUID->"80e01cee-1c17-4ed8-\
8eb9-09fb52f4591f"],

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
sionUUID->"9d7e0dec-bb0a-47cf-8c73-41a642bdb752"]
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
Cell[1510, 35, 98, 0, 98, "Title",ExpressionUUID->"5298ff9f-40e3-4ff3-9ff5-2440ed9a9cf1"],
Cell[1611, 37, 225, 3, 35, "Text",ExpressionUUID->"24b28408-ec66-4eae-b220-c852c93b1e71"],
Cell[CellGroupData[{
Cell[1861, 44, 89, 0, 65, "Subchapter",ExpressionUUID->"8755134a-6eac-4bc9-9856-401e29b9e064"],
Cell[1953, 46, 128, 0, 35, "Text",ExpressionUUID->"056026c3-8a33-49c7-b3b7-a68cbc0ac675"],
Cell[2084, 48, 831, 24, 78, "Input",ExpressionUUID->"dc9fe8b2-9529-46ab-aae2-06f339f999e9"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2952, 77, 199, 5, 65, "Subchapter",ExpressionUUID->"58bb189a-cc09-4d6a-a254-5b13ccf3bed0"],
Cell[3154, 84, 92, 0, 35, "Text",ExpressionUUID->"d9462527-7939-40cf-ba01-7789f1fd5955"],
Cell[CellGroupData[{
Cell[3271, 88, 469, 12, 31, "Input",ExpressionUUID->"6b827d1f-61a7-48b4-acee-01acb3d3954e"],
Cell[3743, 102, 1988, 60, 57, "Output",ExpressionUUID->"94aeb62a-5a2a-4464-9b2f-1ca2eb1eb8e8"]
}, Open  ]],
Cell[5746, 165, 263, 4, 35, "Text",ExpressionUUID->"94f5f41b-d930-4e06-8c74-83953ea28a46"],
Cell[CellGroupData[{
Cell[6034, 173, 274, 7, 31, "Input",ExpressionUUID->"c12493b6-9170-4608-8801-8cb0460c139a"],
Cell[6311, 182, 917, 26, 62, "Output",ExpressionUUID->"2603acf5-ab57-44c5-bd37-c99dcfd7e528"]
}, Open  ]],
Cell[7243, 211, 224, 3, 35, "Text",ExpressionUUID->"d7d137d3-fb68-45c0-af29-0ef1092966c9"],
Cell[CellGroupData[{
Cell[7492, 218, 220, 5, 31, "Input",ExpressionUUID->"4125e9a8-b744-4712-827e-a2d34cef7ee8"],
Cell[7715, 225, 1146, 29, 62, "Output",ExpressionUUID->"c2849a23-df3d-40e8-ba94-e2f1905a2bb3"]
}, Open  ]],
Cell[8876, 257, 126, 0, 35, "Text",ExpressionUUID->"2ebd94f1-c81d-4e6a-9f35-d2dcfaab10dd"],
Cell[CellGroupData[{
Cell[9027, 261, 276, 8, 31, "Input",ExpressionUUID->"ab7fa8e1-43c0-4313-a12e-7fef26d8df28"],
Cell[9306, 271, 634, 18, 57, "Output",ExpressionUUID->"e734f0d3-edff-4e5c-b679-b1894b2acf5d"]
}, Open  ]],
Cell[9955, 292, 143, 2, 35, "Text",ExpressionUUID->"c6115c4f-bb63-4a3c-b860-486cb0267a60"],
Cell[CellGroupData[{
Cell[10123, 298, 311, 8, 31, "Input",ExpressionUUID->"55ad8c49-881b-441c-bf61-6f7bd56c9d5c"],
Cell[10437, 308, 631, 18, 57, "Output",ExpressionUUID->"b9052ca2-35f6-4b36-98f9-496a13888e33"]
}, Open  ]],
Cell[11083, 329, 202, 3, 35, "Text",ExpressionUUID->"43561019-7baf-4587-bfa1-2b47292636cd"],
Cell[CellGroupData[{
Cell[11310, 336, 192, 4, 31, "Input",ExpressionUUID->"2faf48f4-6027-45bf-ae27-6f92e623d7e7"],
Cell[11505, 342, 551, 16, 56, "Output",ExpressionUUID->"e9388348-3989-43dd-9b08-f32ced23ddf5"]
}, Open  ]],
Cell[12071, 361, 191, 3, 35, "Text",ExpressionUUID->"797fd1d9-78e5-4e83-8445-62a0c5784161"],
Cell[CellGroupData[{
Cell[12287, 368, 404, 12, 31, "Input",ExpressionUUID->"67dd637d-88ee-48d2-894f-827b68820c69"],
Cell[12694, 382, 433, 13, 57, "Output",ExpressionUUID->"11a7d8ef-9f41-4ee8-8c8c-c345411e019f"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[13176, 401, 147, 2, 65, "Subchapter",ExpressionUUID->"3f97b281-6ea6-4e21-8806-b978700e2f3b"],
Cell[13326, 405, 146, 2, 35, "Text",ExpressionUUID->"0ed46c9a-8ac6-459e-8c9a-9fc91b86c2f2"],
Cell[CellGroupData[{
Cell[13497, 411, 420, 12, 31, "Input",ExpressionUUID->"96fbff50-f1ee-4815-af19-eb7ad1af2fac"],
Cell[13920, 425, 156, 4, 55, "Output",ExpressionUUID->"311a258c-a2a2-4c69-b5cb-b0747fb9a29b"]
}, Open  ]],
Cell[14091, 432, 429, 12, 35, "Text",ExpressionUUID->"28b46dfe-2e26-4720-a085-be09b6746c44"],
Cell[CellGroupData[{
Cell[14545, 448, 148, 3, 31, "Input",ExpressionUUID->"f246dc43-d23a-4efc-be3f-0a6a279defd0"],
Cell[14696, 453, 123, 1, 35, "Output",ExpressionUUID->"4ccee5c8-1f8c-41df-9e00-184eed409b65"]
}, Open  ]],
Cell[14834, 457, 202, 3, 35, "Text",ExpressionUUID->"bbeccc32-74d8-4f32-8ede-fec4631408cd"]
}, Open  ]],
Cell[CellGroupData[{
Cell[15073, 465, 129, 0, 65, "Subchapter",ExpressionUUID->"552fee79-005c-4d05-8c9f-927c576e1219"],
Cell[15205, 467, 148, 2, 35, "Text",ExpressionUUID->"a2f8e18f-10e4-4548-bd66-c27d082b85dd"],
Cell[CellGroupData[{
Cell[15378, 473, 377, 11, 31, "Input",ExpressionUUID->"60ae074c-51bf-4baa-bcfa-f918f09de36b"],
Cell[15758, 486, 160, 4, 54, "Output",ExpressionUUID->"ad5c81a8-c807-4132-bc76-be31d0823541"]
}, Open  ]],
Cell[15933, 493, 445, 9, 42, "Text",ExpressionUUID->"c6401ba6-f4b0-452e-8856-ebb257e30266"],
Cell[16381, 504, 190, 5, 31, "Input",ExpressionUUID->"0ff25392-1e02-4f1c-8d5f-870cac6b9dcb"],
Cell[16574, 511, 187, 3, 35, "Text",ExpressionUUID->"7ff1c590-5417-4dba-9902-ca5c4d0a2483"],
Cell[CellGroupData[{
Cell[16786, 518, 336, 9, 31, "Input",ExpressionUUID->"80f24613-750f-420b-9e12-d680f5e59793"],
Cell[17125, 529, 108, 1, 35, "Output",ExpressionUUID->"860c6b87-c7fe-46f7-b0ec-8b585810a8d5"]
}, Open  ]],
Cell[17248, 533, 206, 3, 35, "Text",ExpressionUUID->"909ae6c2-db01-4a10-8a85-db987ba0339f"]
}, Open  ]]
}, Open  ]],
Cell[CellGroupData[{
Cell[17503, 542, 132, 0, 98, "Title",ExpressionUUID->"7d0f7f18-6e0a-4e2a-9066-03cff71fcb29"],
Cell[17638, 544, 373, 5, 58, "Text",ExpressionUUID->"6f31b8c4-bded-42fa-a450-6d1a6300b25d"],
Cell[CellGroupData[{
Cell[18036, 553, 89, 0, 65, "Subchapter",ExpressionUUID->"88338340-9a48-4924-8459-9f219097d12d"],
Cell[18128, 555, 190, 3, 35, "Text",ExpressionUUID->"a364b61c-5d98-40eb-b5c2-dc15b3ecf97f"],
Cell[18321, 560, 1084, 32, 101, "Input",ExpressionUUID->"616c8ac9-7672-4c52-ade4-3b23849df1e4"]
}, Open  ]],
Cell[CellGroupData[{
Cell[19442, 597, 199, 5, 65, "Subchapter",ExpressionUUID->"affc0497-f7e3-4cee-90d9-5078aa44af58"],
Cell[19644, 604, 240, 4, 35, "Text",ExpressionUUID->"ad6b50bc-9bb4-4fc5-a33e-2b34bb14fa1d"],
Cell[CellGroupData[{
Cell[19909, 612, 539, 13, 31, "Input",ExpressionUUID->"be05f869-f9a4-4cf9-8773-0ff3728639d6"],
Cell[20451, 627, 3354, 88, 109, "Output",ExpressionUUID->"bff78f72-1712-46a1-a146-09b293eba06c"]
}, Open  ]],
Cell[23820, 718, 263, 4, 35, "Text",ExpressionUUID->"6e75533e-d8c3-4692-92e9-e8b9eace5799"],
Cell[CellGroupData[{
Cell[24108, 726, 301, 8, 31, "Input",ExpressionUUID->"90a41a83-9199-4efe-843f-417e09b98fdf"],
Cell[24412, 736, 1354, 36, 66, "Output",ExpressionUUID->"f7ff57b5-93ba-4931-ab89-1ae08920899f"]
}, Open  ]],
Cell[25781, 775, 224, 3, 35, "Text",ExpressionUUID->"8013d08e-3c70-4155-ab65-150dbe167d38"],
Cell[CellGroupData[{
Cell[26030, 782, 220, 5, 31, "Input",ExpressionUUID->"74071871-abfb-4360-8856-91dbb825893b"],
Cell[26253, 789, 1603, 38, 62, "Output",ExpressionUUID->"de91ae4f-26d7-4be3-b18b-8024ead62fd6"]
}, Open  ]],
Cell[27871, 830, 126, 0, 35, "Text",ExpressionUUID->"8e8fea37-3976-4850-84c0-18b42548008b"],
Cell[CellGroupData[{
Cell[28022, 834, 276, 8, 31, "Input",ExpressionUUID->"d9f7fe85-21e2-410e-b51f-ccb8d1c73291"],
Cell[28301, 844, 850, 23, 58, "Output",ExpressionUUID->"a7f6cd0a-d839-496d-8316-e6e0c6649831"]
}, Open  ]],
Cell[29166, 870, 220, 3, 35, "Text",ExpressionUUID->"a1f51013-e7af-4658-b112-472e5f41e89e"],
Cell[CellGroupData[{
Cell[29411, 877, 311, 8, 31, "Input",ExpressionUUID->"00e919c2-e7c5-48ff-aff3-6ede77e6b2b0"],
Cell[29725, 887, 856, 23, 58, "Output",ExpressionUUID->"5d064e5b-11f0-4474-b39a-23355ebd63e6"]
}, Open  ]],
Cell[30596, 913, 202, 3, 35, "Text",ExpressionUUID->"79c38174-bf34-4b23-8e84-b9796aeaa4aa"],
Cell[CellGroupData[{
Cell[30823, 920, 192, 4, 31, "Input",ExpressionUUID->"04633174-1cc2-4203-a690-b6013776d0c1"],
Cell[31018, 926, 783, 22, 56, "Output",ExpressionUUID->"d6625080-db41-4e35-b925-4d8d16785147"]
}, Open  ]],
Cell[31816, 951, 329, 5, 58, "Text",ExpressionUUID->"cd9ce0e8-dd96-4f76-bad7-c82aef6f32da"],
Cell[CellGroupData[{
Cell[32170, 960, 404, 12, 31, "Input",ExpressionUUID->"80e01cee-1c17-4ed8-8eb9-09fb52f4591f"],
Cell[32577, 974, 433, 13, 90, "Output",ExpressionUUID->"9d7e0dec-bb0a-47cf-8c73-41a642bdb752"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}
]
*)

(* End of internal cache information *)

(* NotebookSignature Vu03y3LOM1jbZAwILG15mxtQ *)
