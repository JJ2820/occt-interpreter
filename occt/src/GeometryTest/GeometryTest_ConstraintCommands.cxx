// Created on: 1993-08-12
// Created by: Bruno DUMORTIER
// Copyright (c) 1993-1999 Matra Datavision
// Copyright (c) 1999-2014 OPEN CASCADE SAS
//
// This file is part of Open CASCADE Technology software library.
//
// This library is free software; you can redistribute it and/or modify it under
// the terms of the GNU Lesser General Public License version 2.1 as published
// by the Free Software Foundation, with special exception defined in the file
// OCCT_LGPL_EXCEPTION.txt. Consult the file LICENSE_LGPL_21.txt included in OCCT
// distribution for complete text of the license and disclaimer of any warranty.
//
// Alternatively, this file may be used under the terms of Open CASCADE
// commercial license or contractual agreement.

// JPI : Commande smooth transferee dans GeomliteTest
// PMN : Ajout de la commande smooth
// JCT : Correction d'un trap dans la commande gcarc

#include <Standard_Stream.hxx>

#include <GeometryTest.hxx>
#include <DrawTrSurf.hxx>
#include <Draw.hxx>

#include <Draw_Interpretor.hxx>
#include <Geom2dGcc_Circ2d2TanRad.hxx>
#include <Geom2dGcc_Circ2d3Tan.hxx>
#include <Geom2dGcc_Circ2dTanCen.hxx>
#include <Geom2dGcc_Lin2d2Tan.hxx>
#include <Geom2dGcc_Lin2dTanObl.hxx>
#include <Geom2dGcc.hxx>
#include <Geom2dGcc_QualifiedCurve.hxx>
#include <Geom2d_CartesianPoint.hxx>
#include <Geom2d_Circle.hxx>
#include <Geom2d_Line.hxx>
#include <GeomAPI_Interpolate.hxx>
#include <Geom2dAPI_Interpolate.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <Geom_BSplineCurve.hxx>
#include <TColgp_HArray1OfPnt2d.hxx>
#include <Geom2d_BSplineCurve.hxx>
#include <Geom_TrimmedCurve.hxx>
#include <TColgp_HArray1OfPnt.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <TColgp_Array1OfPnt2d.hxx>
#include <TColgp_HArray1OfVec.hxx>
#include <TColgp_Array1OfVec.hxx>
#include <TColStd_HArray1OfBoolean.hxx>
#include <AppParCurves_MultiBSpCurve.hxx>
#include <GC_MakeSegment.hxx>
#include <GC_MakeArcOfCircle.hxx>

#include <stdio.h>
// Standard_IMPORT Draw_Color DrawTrSurf_CurveColor(const Draw_Color);


static Standard_Integer solutions(Draw_Interpretor& di,
                                  Geom2dGcc_Circ2d2TanRad& ct3, const char* name) 
{
  char solname[200];

  // Draw_Color col = DrawTrSurf_CurveColor(Draw_Color(Draw_vert));
  // DrawTrSurf_CurveColor(col);

  if (ct3.IsDone()) {
    for (Standard_Integer i = 1 ; i <= ct3.NbSolutions() ; i++) {
      Handle(Geom2d_Circle) C = new Geom2d_Circle(ct3.ThisSolution(i));
      Sprintf(solname,"%s_%d",name,i);
      DrawTrSurf::Set(solname, C);
      di << solname << " ";
    }
    return 0;
  }
  else {
    di << "Circ2d2TanRad Not done";
    return 1;
  }
}

static Standard_Integer solutions(Draw_Interpretor& di,
                                  Geom2dGcc_Circ2d3Tan& ct3, const char* name) 
{
  char solname[200];

  // Draw_Color col = DrawTrSurf_CurveColor(Draw_Color(Draw_vert));
  // DrawTrSurf_CurveColor(col);

  if (ct3.IsDone()) {
    for (Standard_Integer i = 1 ; i <= ct3.NbSolutions() ; i++) {
      Handle(Geom2d_Circle) C = new Geom2d_Circle(ct3.ThisSolution(i));
      Sprintf(solname,"%s_%d",name,i);
      DrawTrSurf::Set(solname, C);
      di << solname << " ";
    }
    return 0;
  }
  else {
    di << "Circ2d3Tan Not done";
    return 1;
  }
}

//=======================================================================
//function : solutions
//purpose  : 
//=======================================================================
static Standard_Integer solutions(Draw_Interpretor& theDI,
                                  Geom2dGcc_Circ2dTanCen& theCt2,
                                  const char* theName)
{
  char solname[200];

  // Draw_Color col = DrawTrSurf_CurveColor(Draw_Color(Draw_vert));
  // DrawTrSurf_CurveColor(col);

  if (theCt2.IsDone())
  {
    for (Standard_Integer i = 1; i <= theCt2.NbSolutions(); i++)
    {
      Handle(Geom2d_Circle) C = new Geom2d_Circle(theCt2.ThisSolution(i));
      Sprintf(solname, "%s_%d", theName, i);
      DrawTrSurf::Set(solname, C);
      theDI << solname << " ";
    }
    return 0;
  }
  else
  {
    theDI << "Circ2dTanCen Not done";
    return 1;
  }
}

//=======================================================================
//function : Cirtang
//purpose  : 
//=======================================================================
static Standard_Integer Cirtang(Draw_Interpretor& theDI,
                                Standard_Integer theNArgs,
                                const char** theArgVals)
{
  if (theNArgs < 3)
  {
    theDI << "Use: " << theArgVals[0] << "result [-t <Tolerance>] -c <curve> -p <point> -r <Radius>...\n";
    return 1;
  }

  Standard_Real aTol = Precision::Confusion();
  Handle(Geom2d_Curve) aC[3];
  gp_Pnt2d aP[3];
  Standard_Real aRadius = -1.0;

  Standard_Integer aNbCurves = 0, aNbPnts = 0;

  for (Standard_Integer anArgID = 2; anArgID < theNArgs; anArgID++)
  {
    if (theArgVals[anArgID][0] != '-')
    {
      theDI << "Cannot interpret the argument #" << anArgID << " (" << theArgVals[anArgID] << ")\n";
      return 1;
    }
    else if (!strcmp(theArgVals[anArgID], "-c"))
    {
      if (aNbCurves >= 3)
      {
        theDI << "A lot of curves are given (not greater than 3 ones are expected)\n";
        return 1;
      }

      aC[aNbCurves] = DrawTrSurf::GetCurve2d(theArgVals[++anArgID]);
      if (aC[aNbCurves].IsNull())
      {
        theDI << "Error: " << theArgVals[anArgID] << " is not a curve\n";
        return 1;
      }

      aNbCurves++;
    }
    else if (!strcmp(theArgVals[anArgID], "-p"))
    {
      if (aNbPnts >= 3)
      {
        theDI << "A lot of points are given (not greater than 3 ones are expected)\n";
        return 1;
      }

      if (!DrawTrSurf::GetPoint2d(theArgVals[++anArgID], aP[aNbPnts]))
      {
        theDI << "Error: " << theArgVals[anArgID] << " is not a point\n";
        return 1;
      }

      aNbPnts++;
    }
    else if (!strcmp(theArgVals[anArgID], "-r"))
    {
      aRadius = Draw::Atof(theArgVals[++anArgID]);
    }
    else if (!strcmp(theArgVals[anArgID], "-t"))
    {
      aTol = Draw::Atof(theArgVals[++anArgID]);
    }
    else
    {
      theDI << "Unknown option " << theArgVals[anArgID] << "\n";
      return 1;
    }
  }

  if (aNbCurves == 3)
  {
    // C-C-C
    Geom2dGcc_Circ2d3Tan aCt3(Geom2dGcc::Unqualified(aC[0]),
                              Geom2dGcc::Unqualified(aC[1]),
                              Geom2dGcc::Unqualified(aC[2]),
                              aTol, 0, 0, 0);
    theDI << "Solution of type C-C-C is: ";
    return solutions(theDI, aCt3, theArgVals[1]);
  }
  else if (aNbCurves == 2)
  {
    if (aNbPnts >= 1)
    {
      // C-C-P
      Geom2dGcc_Circ2d3Tan aCt3(Geom2dGcc::Unqualified(aC[0]),
                                Geom2dGcc::Unqualified(aC[1]),
                                new Geom2d_CartesianPoint(aP[0]),
                                aTol, 0, 0);
      theDI << "Solution of type C-C-P is: ";
      return solutions(theDI, aCt3, theArgVals[1]);
    }
    else if (aRadius > 0)
    {
      // C-C-R
      Geom2dGcc_Circ2d2TanRad aCt3(Geom2dGcc::Unqualified(aC[0]),
                                   Geom2dGcc::Unqualified(aC[1]),
                                   aRadius, aTol);
      theDI << "Solution of type C-C-R is: ";
      return solutions(theDI, aCt3, theArgVals[1]);
    }

    theDI << "Error: Unsupported set of input data!\n";
    return 1;
  }
  else if (aNbCurves == 1)
  {
    if (aNbPnts == 2)
    {
      //C-P-P
      Geom2dGcc_Circ2d3Tan aCt3(Geom2dGcc::Unqualified(aC[0]),
                                new Geom2d_CartesianPoint(aP[0]),
                                new Geom2d_CartesianPoint(aP[1]),
                                aTol,0);
      theDI << "Solution of type C-P-P is: ";
      return solutions(theDI, aCt3, theArgVals[1]);
    }
    else if (aNbPnts == 1)
    {
      if (aRadius > 0.0)
      {
        //C-P-R
        Geom2dGcc_Circ2d2TanRad aCt3(Geom2dGcc::Unqualified(aC[0]),
                                     new Geom2d_CartesianPoint(aP[0]),
                                     aRadius, aTol);
        theDI << "Solution of type C-P-R is: ";
        return solutions(theDI, aCt3, theArgVals[1]);
      }
      else
      {
        // C-P
        Geom2dGcc_Circ2dTanCen aCt2(Geom2dGcc::Unqualified(aC[0]),
                                    new Geom2d_CartesianPoint(aP[0]), aTol);
        theDI << "Solution of type C-P is: ";
        return solutions(theDI, aCt2, theArgVals[1]);
      }
    }

    theDI << "Error: Unsupported set of input data!\n";
    return 1;
  }
  else if (aNbPnts >= 2)
  {
    if (aNbPnts == 3)
    {
      //P-P-P
      Geom2dGcc_Circ2d3Tan aCt3(new Geom2d_CartesianPoint(aP[0]),
                                new Geom2d_CartesianPoint(aP[1]),
                                new Geom2d_CartesianPoint(aP[2]),
                                aTol);
      theDI << "Solution of type P-P-P is: ";
      return solutions(theDI, aCt3, theArgVals[1]);
    }
    else if (aRadius > 0)
    {
      //P-P-R
      Geom2dGcc_Circ2d2TanRad aCt3(new Geom2d_CartesianPoint(aP[0]),
                                   new Geom2d_CartesianPoint(aP[1]),
                                   aRadius, aTol);
      theDI << "Solution of type P-P-R is: ";
      return solutions(theDI, aCt3, theArgVals[1]);
    }

    theDI << "Error: Unsupported set of input data!\n";
    return 1;
  }

  theDI << "Error: Unsupported set of input data!\n";
  return 1;
}


//=======================================================================
//function : lintang
//purpose  : 
//=======================================================================

static Standard_Integer lintang (Draw_Interpretor& di,Standard_Integer n, const char** a)
{
  if (n < 4) return 1;

  Handle(Geom2d_Curve) C1 = DrawTrSurf::GetCurve2d(a[2]);
  Handle(Geom2d_Curve) C2 = DrawTrSurf::GetCurve2d(a[3]);

  char solname[200];

  if (C1.IsNull() || C2.IsNull())
    return 1;

  // Draw_Color col = DrawTrSurf_CurveColor(Draw_Color(Draw_vert));

  if (n >= 5) {
    Handle(Geom2d_Line) L = Handle(Geom2d_Line)::DownCast(C2);
    if (L.IsNull()) {
      di << "Second argument must be a line";
      return 1;
    }
    Standard_Real ang = Draw::Atof(a[4]) * (M_PI / 180.0);
    Geom2dGcc_Lin2dTanObl ct3(Geom2dGcc::Unqualified(C1),
      L->Lin2d(),
      Precision::Angular(),
      (C1->FirstParameter()+C1->LastParameter())/2.,
      ang);
    if (ct3.IsDone()) {
      for (Standard_Integer i = 1 ; i <= ct3.NbSolutions() ; i++) {
        Handle(Geom2d_Line) LS = new Geom2d_Line(ct3.ThisSolution(i));
        Sprintf(solname,"%s_%d",a[1],i);
        char* temp = solname; // pour portage WNT
        DrawTrSurf::Set(temp,LS);
        di << solname << " ";
      }
    }
    else
      di << "Lin2dTanObl Not done\n";
  }
  else {
    Geom2dGcc_Lin2d2Tan ct3(Geom2dGcc::Unqualified(C1),
      Geom2dGcc::Unqualified(C2),
      Precision::Angular(),
      (C1->FirstParameter()+C1->LastParameter())/2.,
      (C2->FirstParameter()+C2->LastParameter())/2.);
    if (ct3.IsDone()) {
      for (Standard_Integer i = 1 ; i <= ct3.NbSolutions() ; i++) {
        Handle(Geom2d_Line) LS = new Geom2d_Line(ct3.ThisSolution(i));
        Sprintf(solname,"%s_%d",a[1],i);
        char* temp = solname; // pour portage WNT
        DrawTrSurf::Set(temp,LS);
        di << solname << " ";
      }
    }
    else
      di << "Lin2d2Tan Not done\n";
  }

  // DrawTrSurf_CurveColor(col);

  return 0;
}

//==================================================================================
static Standard_Integer interpol (Draw_Interpretor& di,Standard_Integer n, const char** a)
//==================================================================================
{
  return 0;
}

static Standard_Integer tanginterpol (Draw_Interpretor& di,
                                      Standard_Integer n, 
                                      const char** a)
{


  if (n < 4)
    return 1;

  Standard_Integer 
    ii,
    jj,
    //    num_knots,
    //    degree,
    num_tangents,
    num_read,
    num_start,
    num_parameters ;


  Standard_Real 
    //    delta,
    tolerance;
  //    parameter ;

  Standard_Boolean periodic_flag = Standard_False ;
  gp_Pnt a_point ;
  gp_Vec a_vector ;
  tolerance = 1.0e-5 ;




  Handle(Geom_BSplineCurve) NewCurvePtr ;




  num_read = 2 ;
  if (strcmp(a[num_read],"p") == 0) {
    periodic_flag = Standard_True ;
    num_read += 1 ;
  }
  num_parameters = Draw::Atoi(a[num_read]) ;

  if (num_parameters < 2) {
    num_parameters = 2 ;
  }
  if ( n <  num_parameters * 3 + num_read) {
    return 1 ;
  }
  Handle(TColgp_HArray1OfPnt)   PointsArrayPtr=
    new TColgp_HArray1OfPnt(1,num_parameters) ;

  num_tangents = ((n - num_read) / 3)  - num_parameters ;
  num_tangents = Max (0,num_tangents) ; 
  num_tangents = Min (num_parameters, num_tangents) ;
  ii = 1 ;
  num_start = num_read ;
  num_read += 1 ;
  while (num_read <= num_parameters * 3 + num_start ) {
    for (jj = 1 ; jj <= 3 ; jj++) {
      a_point.SetCoord(jj,Draw::Atof(a[num_read])) ;
      num_read += 1 ;
    }
    PointsArrayPtr->SetValue(ii,a_point) ;
    ii += 1 ;
  }
  GeomAPI_Interpolate anInterpolator(PointsArrayPtr,
    periodic_flag,
    tolerance) ; 

  if (num_tangents > 0) {
    TColgp_Array1OfVec TangentsArray(1,num_parameters) ;
    Handle(TColStd_HArray1OfBoolean) 
      TangentFlagsPtr =
      new TColStd_HArray1OfBoolean(1,num_parameters) ;

    for (ii = 1 ; ii <= num_tangents ; ii++) {
      TangentFlagsPtr->SetValue(ii,Standard_True) ;
    }
    for (ii = num_tangents + 1 ; ii <= num_parameters ; ii++) {
      TangentFlagsPtr->SetValue(ii,Standard_False) ;
    }
    ii = 1 ;
    while (ii <= num_tangents) {
      for (jj = 1 ; jj <= 3 ; jj++) {
        a_vector.SetCoord(jj,Draw::Atof(a[num_read])) ;
        num_read += 1 ;
      }
      TangentsArray.SetValue(ii,a_vector) ;
      ii += 1 ;
    }


    anInterpolator.Load(TangentsArray,
      TangentFlagsPtr) ;
  }
  anInterpolator.Perform() ;
  if (anInterpolator.IsDone()) {
    NewCurvePtr =
      anInterpolator.Curve() ;

    DrawTrSurf::Set(a[1],
      NewCurvePtr) ;
    di << a[2] << " " ;

  }
  return 0 ;
}

//==================================================================================
static Standard_Integer gcarc (Draw_Interpretor& di,Standard_Integer n, const char** a)
//==================================================================================
{
  if (n >= 5) {
    gp_Pnt P1,P2,P3,P4;
    if (!strcmp(a[2], "seg")) {
      if (DrawTrSurf::GetPoint(a[3], P1)) {
        if (DrawTrSurf::GetPoint(a[4], P2)) {
          Handle(Geom_Curve) theline (GC_MakeSegment(P1,P2).Value());
          DrawTrSurf::Set(a[1], theline);
          return 1;
        }
      }
    }
    else if (!strcmp(a[2], "cir")) {
      if (DrawTrSurf::GetPoint(a[3], P1)) {
        if (DrawTrSurf::GetPoint(a[4], P2)) {
          if (DrawTrSurf::GetPoint(a[5], P3)) {
            //	    if (DrawTrSurf::GetPoint(a[6], P4)) {
            if (n>6) {
              DrawTrSurf::GetPoint(a[6], P4);
              gp_Vec V1 = gp_Vec(P2,P3);                                                    
              Handle(Geom_Curve)thearc (GC_MakeArcOfCircle(P1,V1,P4).Value());
              DrawTrSurf::Set(a[1], thearc);
              return 1;
            }
            else {
              Handle(Geom_Curve)thearc (GC_MakeArcOfCircle(P1,P2,P3).Value());
              DrawTrSurf::Set(a[1], thearc);
              return 1;
            }
          }
        }
      }
    }
  }
  di <<"give a name for arc and the type seg or cir then\n";
  di <<"give passing points p1 p2 for seg    p1 p2 p3 or p1 p2 p3 p4 for cir (p2 p3 is a tgtvec)!\n";
  return 0;
}

//=======================================================================
//function : ConstraintCommands
//purpose  : 
//=======================================================================


void  GeometryTest::ConstraintCommands(Draw_Interpretor& theCommands)
{

  static Standard_Boolean loaded = Standard_False;
  if (loaded) return;
  loaded = Standard_True;

  DrawTrSurf::BasicCommands(theCommands);

  const char* g;
  // constrained constructs
  g = "GEOMETRY Constraints";

  theCommands.Add("cirtang",
    "cirtang cname [-t <Tolerance>] -c <curve> -p <point> -r <Radius>...",
    __FILE__,
    Cirtang, g);

  theCommands.Add("lintan",
    "lintan lname curve1 curve2 [angle]",
    __FILE__,
    lintang,g);


  theCommands.Add("interpol",
    "interpol cname [fic]", 
    __FILE__,
    interpol, g);
  theCommands.Add("tanginterpol",
    "tanginterpol curve [p] num_points points [tangents] modifier  p = periodic",
    __FILE__,
    tanginterpol,g);

  theCommands.Add("gcarc",
    "gcarc name seg/cir p1 p2 p3 p4",
    __FILE__,
    gcarc,g);
}
