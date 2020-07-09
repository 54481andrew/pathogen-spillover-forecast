(* Content-type: application/vnd.wolfram.cdf.text *)

(*** Wolfram CDF File ***)
(* http://www.wolfram.com/cdf *)

(* CreatedBy='Mathematica 11.0' *)

(*************************************************************************)
(*                                                                       *)
(*  The Mathematica License under which this file was created prohibits  *)
(*  restricting third parties in receipt of this file from republishing  *)
(*  or redistributing it by any means, including but not limited to      *)
(*  rights management or terms of use, without the express consent of    *)
(*  Wolfram Research, Inc. For additional information concerning CDF     *)
(*  licensing and redistribution see:                                    *)
(*                                                                       *)
(*        www.wolfram.com/cdf/adopting-cdf/licensing-options.html        *)
(*                                                                       *)
(*************************************************************************)

(*CacheID: 234*)
(* Internal cache information:
NotebookFileLineBreakTest
NotebookFileLineBreakTest
NotebookDataPosition[      1064,         20]
NotebookDataLength[      9460,        357]
NotebookOptionsPosition[      8612,        304]
NotebookOutlinePosition[      9115,        325]
CellTagsIndexPosition[      9072,        322]
WindowFrame->Normal*)

(* Beginning of Notebook Content *)
Notebook[{

Cell[CellGroupData[{
Cell["Human Lassa Case Rate", "Title"],

Cell["\<\
This code shows the analyses performed on the system of ordinary differential \
equations in the Appendix. \
\>", "Text"],

Cell[CellGroupData[{

Cell["Define model", "Subchapter"],

Cell["Define the system of ODEs that describe human infection. ", "Text"],

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
    RowBox[{"\[Gamma]", " ", "i"}], " ", "-", " ", 
    RowBox[{"d", " ", "r"}], " ", "-", " ", 
    RowBox[{"\[Lambda]", " ", "r"}]}]}], ";"}]}], "Input"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Derive relationship between seroprevalence and new case rate", \
"Subchapter"],

Cell["Find the equilibrium.", "Text"],

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
      RowBox[{"s", ",", "i", ",", "r"}], "}"}]}], "]"}], "]"}]}]], "Input"],

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
         RowBox[{"d", "+", "F", "+", "\[Gamma]"}], ")"}], " ", 
        "\[Lambda]"}]}]]}], ",", 
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
         RowBox[{"d", "+", "F", "+", "\[Gamma]"}], ")"}], " ", 
        "\[Lambda]"}]}]]}], ",", 
    RowBox[{"r", "\[Rule]", 
     FractionBox[
      RowBox[{"b", " ", "F", " ", "\[Gamma]"}], 
      RowBox[{
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "F"}], ")"}], " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "\[Gamma]"}], ")"}]}], "+", 
       RowBox[{"d", " ", 
        RowBox[{"(", 
         RowBox[{"d", "+", "F", "+", "\[Gamma]"}], ")"}], " ", 
        "\[Lambda]"}]}]]}]}], "}"}], "}"}]], "Output"]
}, Open  ]],

Cell["\<\
Find seroprevalence at equilibrium (this is the data we have to work with)\
\>", "Text"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"Rstar", " ", "=", " ", 
  RowBox[{
   RowBox[{"r", "/", 
    RowBox[{"(", 
     RowBox[{"b", "/", "d"}], ")"}]}], "/.", "eq"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox[
   RowBox[{"d", " ", "F", " ", "\[Gamma]"}], 
   RowBox[{
    RowBox[{"d", " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "F"}], ")"}], " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "\[Gamma]"}], ")"}]}], "+", 
    RowBox[{"d", " ", 
     RowBox[{"(", 
      RowBox[{"d", "+", "F", "+", "\[Gamma]"}], ")"}], " ", "\[Lambda]"}]}]], 
  "}"}]], "Output"]
}, Open  ]],

Cell["\<\
Solve for the force of infection, F, in terms of seroprevalence.\
\>", "Text"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"foi", " ", "=", " ", 
  RowBox[{
   RowBox[{"FullSimplify", "[", 
    RowBox[{"Solve", "[", 
     RowBox[{
      RowBox[{"hsero", " ", "\[Equal]", " ", "Rstar"}], ",", " ", "F"}], 
     "]"}], "]"}], "//", "Flatten"}]}]], "Input"],

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
      RowBox[{"-", "\[Gamma]"}], "+", 
      RowBox[{"hsero", " ", 
       RowBox[{"(", 
        RowBox[{"d", "+", "\[Gamma]", "+", "\[Lambda]"}], ")"}]}]}]]}]}], 
  "}"}]], "Output"]
}, Open  ]],

Cell["\<\
The rate of new cases of Lassa is F*s. Calculate this, and rewrite everything \
in terms of seroprevalence.\
\>", "Text"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"NewCaseRate", " ", "=", " ", 
  RowBox[{"FullSimplify", "[", 
   RowBox[{
    RowBox[{"(", 
     RowBox[{"F", "/.", "foi"}], ")"}], "*", 
    RowBox[{"(", 
     RowBox[{
      RowBox[{"s", "/.", "eq"}], "/.", "foi"}], ")"}]}], "]"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox[
   RowBox[{"b", " ", "hsero", " ", 
    RowBox[{"(", 
     RowBox[{"d", "+", "\[Gamma]"}], ")"}], " ", 
    RowBox[{"(", 
     RowBox[{"d", "+", "\[Lambda]"}], ")"}]}], 
   RowBox[{"d", " ", "\[Gamma]"}]], "}"}]], "Output"]
}, Open  ]],

Cell["\<\
Find the overall effect of reinfection on case-rate estimate.\
\>", "Text"],

Cell[CellGroupData[{

Cell[BoxData[
 RowBox[{"NewCaseRate", "/", 
  RowBox[{"(", 
   RowBox[{"NewCaseRate", "/.", 
    RowBox[{"\[Lambda]", "\[Rule]", "0"}]}], ")"}]}]], "Input"],

Cell[BoxData[
 RowBox[{"{", 
  FractionBox[
   RowBox[{"d", "+", "\[Lambda]"}], "d"], "}"}]], "Output"]
}, Open  ]],

Cell[TextData[{
 "Reinfection multiplies the case rate by ",
 Cell[BoxData[
  FormBox[
   FractionBox[
    RowBox[{"d", "+", "\[Lambda]"}], "d"], TraditionalForm]]]
}], "Text"]
}, Open  ]],

Cell[CellGroupData[{

Cell["Estimate the effect of reinfection on total case numbers", "Subchapter"],

Cell["\<\
Assume the population of west Africa is 11 million, 16% of individuals have \
Lassa antibodies (average across all data), lifespan is 50 years, and Lassa \
recovery rate is 1 month. \
\>", "Text"],

Cell[BoxData[
 RowBox[{
  RowBox[{"Pars", "=", 
   RowBox[{"{", 
    RowBox[{
     RowBox[{"hsero", "\[Rule]", "0.16"}], ",", " ", 
     RowBox[{"d", "\[Rule]", " ", "0.02"}], ",", " ", 
     RowBox[{"\[Gamma]", "\[Rule]", " ", "12"}]}], "}"}]}], ";"}]], "Input"],

Cell["\<\
Find the ratio of estimates of new cases with (\[Lambda] = 0.064) and without \
(\[Lambda] = 0) seroreversion. \
\>", "Text"],

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
  "Pars"}]], "Input"],

Cell[BoxData[
 RowBox[{"{", "4.2`", "}"}]], "Output"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
},
WindowSize->{1013, 887},
WindowMargins->{{Automatic, 58}, {Automatic, 18}},
Visible->True,
ScrollingOptions->{"VerticalScrollRange"->Fit},
ShowCellBracket->Automatic,
CellContext->Notebook,
TrackCellChangeTimes->False,
FrontEndVersion->"11.0 for Mac OS X x86 (32-bit, 64-bit Kernel) (September \
21, 2016)",
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
Cell[1486, 35, 38, 0, 92, "Title"],
Cell[1527, 37, 131, 3, 30, "Text"],
Cell[CellGroupData[{
Cell[1683, 44, 34, 0, 63, "Subchapter"],
Cell[1720, 46, 73, 0, 30, "Text"],
Cell[1796, 48, 702, 20, 75, "Input"]
}, Open  ]],
Cell[CellGroupData[{
Cell[2535, 73, 84, 1, 63, "Subchapter"],
Cell[2622, 76, 37, 0, 30, "Text"],
Cell[CellGroupData[{
Cell[2684, 80, 412, 11, 32, "Input"],
Cell[3099, 93, 1547, 48, 47, "Output"]
}, Open  ]],
Cell[4661, 144, 98, 2, 30, "Text"],
Cell[CellGroupData[{
Cell[4784, 150, 163, 5, 32, "Input"],
Cell[4950, 157, 400, 13, 47, "Output"]
}, Open  ]],
Cell[5365, 173, 88, 2, 30, "Text"],
Cell[CellGroupData[{
Cell[5478, 179, 254, 7, 32, "Input"],
Cell[5735, 188, 457, 15, 47, "Output"]
}, Open  ]],
Cell[6207, 206, 131, 3, 30, "Text"],
Cell[CellGroupData[{
Cell[6363, 213, 269, 8, 32, "Input"],
Cell[6635, 223, 267, 8, 47, "Output"]
}, Open  ]],
Cell[6917, 234, 85, 2, 30, "Text"],
Cell[CellGroupData[{
Cell[7027, 240, 156, 4, 32, "Input"],
Cell[7186, 246, 103, 3, 46, "Output"]
}, Open  ]],
Cell[7304, 252, 176, 6, 39, "Text"]
}, Open  ]],
Cell[CellGroupData[{
Cell[7517, 263, 78, 0, 63, "Subchapter"],
Cell[7598, 265, 206, 4, 49, "Text"],
Cell[7807, 271, 263, 7, 32, "Input"],
Cell[8073, 280, 135, 3, 30, "Text"],
Cell[CellGroupData[{
Cell[8233, 287, 283, 9, 32, "Input"],
Cell[8519, 298, 53, 1, 32, "Output"]
}, Open  ]]
}, Open  ]]
}, Open  ]]
}
]
*)

(* NotebookSignature lxpyydmE4ZcRLAg4cIiwVI03 *)
