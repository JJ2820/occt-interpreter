#ifndef E0_IO_INTERROGATE_H
#define E0_IO_INTERROGATE_H

#include <Precision.hxx>
#include <TopoDS.hxx>
#include <TopoDS_Edge.hxx>
#include <TopoDS_Face.hxx>
#include <TopoDS_Compound.hxx>
#include <BRep_Builder.hxx>
#include <BRepTools.hxx>
#include <BRep_Tool.hxx>
#include <TopExp_Explorer.hxx>
#include <TopLoc_Location.hxx>
#include <Poly_Triangulation.hxx>
#include <TColgp_Array1OfPnt.hxx>
#include <BRepMesh_IncrementalMesh.hxx>
#include <Geom_Surface.hxx>
#include "data.hpp"
#include "commonIO.hpp" // /home/jj/ideeza/occt-interpreter/shape-io/commonIO.hpp
#include <BRepBndLib.hxx>
#include <Geom_Plane.hxx>
#include "surfaceIO.hpp" //home/jj/ideeza/occt-interpreter/shape-io/surfaceIO.hpp


namespace e0 {
namespace io {

// Memory-optimized tessellation writer
void writeFaceTessellation(const Handle(Poly_Triangulation)& aTr, 
                         const TopLoc_Location& aLocation,
                         const Handle(Geom_Surface)& aSurface, 
                         DATA& tessOut) {
  
  // Add bounds checking
  if (aTr.IsNull() || aTr->NbTriangles() == 0 || aTr->NbNodes() == 0) {
    return;
  }

  const Poly_Array1OfTriangle& triangles = aTr->Triangles();
  Standard_Integer nnn = aTr->NbTriangles();
  bool isPlane = aSurface->IsKind("Geom_Plane");

  // Pre-transform points with bounds checking
  std::vector<gp_Pnt> transformedPoints;
  transformedPoints.reserve(aTr->NbNodes());
  for(Standard_Integer i = 1; i <= aTr->NbNodes(); i++) {
    if (i > aTr->NbNodes()) break;  // Extra safety check
    transformedPoints.push_back(aTr->Node(i).Transformed(aLocation));
  }

  // Pre-compute normals with bounds checking
  std::vector<gp_Vec> normals;
  if (!isPlane) {
    normals.reserve(aTr->NbNodes());
    for(Standard_Integer i = 1; i <= aTr->NbNodes(); i++) {
      if (i > aTr->NbNodes()) break;  // Extra safety check
      try {
        gp_Pnt2d uv = aTr->UVNode(i);
        gp_Pnt dummy;
        gp_Vec d1u, d1v;
        aSurface->D1(uv.X(), uv.Y(), dummy, d1u, d1v);
        gp_Vec normal = d1u.Crossed(d1v);
        if (normal.Magnitude() > Precision::Confusion()) {
          normal.Multiply(1.0 / normal.Magnitude());
        } else {
          normal = gp_Vec(0, 0, 1);  // Default normal if calculation fails
        }
        normals.push_back(normal.Transformed(aLocation));
      } catch (Standard_Failure const&) {
        normals.push_back(gp_Vec(0, 0, 1));  // Default normal on failure
      }
    }
  }

  // Process triangles with bounds checking
  const Standard_Integer BATCH_SIZE = 1000;
  for(Standard_Integer batch = 0; batch < nnn; batch += BATCH_SIZE) {
    Standard_Integer batchEnd = std::min(batch + BATCH_SIZE, nnn);
    
    for(Standard_Integer nt = batch + 1; nt <= batchEnd; nt++) {
      try {
        Standard_Integer n1, n2, n3;
        triangles(nt).Get(n1, n2, n3);
        
        // Validate indices
        if (n1 <= 0 || n2 <= 0 || n3 <= 0 || 
            n1 > aTr->NbNodes() || n2 > aTr->NbNodes() || n3 > aTr->NbNodes()) {
          continue;  // Skip invalid triangles
        }
        
        n1--; n2--; n3--;  // Convert to 0-based index

        // Validate transformed points access
        if (n1 >= transformedPoints.size() || n2 >= transformedPoints.size() || 
            n3 >= transformedPoints.size()) {
          continue;
        }

        DATA def = Array();
        DATA tr = Array();
        
        tr.append(pntWrite(transformedPoints[n1]));
        tr.append(pntWrite(transformedPoints[n2]));
        tr.append(pntWrite(transformedPoints[n3]));
        
        def.append(tr);

        if (!isPlane && !normals.empty()) {
          // Validate normals access
          if (n1 >= normals.size() || n2 >= normals.size() || n3 >= normals.size()) {
            continue;
          }
          
          DATA norms = Array();
          norms.append(dirWrite(normals[n1]));
          norms.append(dirWrite(normals[n2]));
          norms.append(dirWrite(normals[n3]));
          def.append(norms);
        }
        
        tessOut.append(def);
      } catch (Standard_Failure const&) {
        continue;  // Skip problematic triangles
      }
    }
  }
}

DATA interrogate(const TopoDS_Shape& aShape, Standard_Real aDeflection = 15, 
                Standard_Boolean INTERROGATE_STRUCT_ONLY = false) {
  DATA out = Object();
  
  try {
    // Validate input shape
    if (aShape.IsNull()) {
      throw Standard_Failure("Null shape provided");
    }

    // Use adaptive deflection with bounds
    Standard_Real actualDeflection = aDeflection;
    if (actualDeflection <= 0 || actualDeflection > 1000) {  // Add reasonable bounds
      //actualDeflection = calculateOptimalDeflection(aShape);
      actualDeflection = 15.0;;
    }

    // Clean existing triangulation
    BRepTools::Clean(aShape);

    // Perform incremental meshing with error handling
    try {
      BRepMesh_IncrementalMesh mesher(aShape, actualDeflection, 
        Standard_True,   // relative
        0.5,            // angular deflection
        Standard_False  // parallel
      );
    } catch (Standard_Failure const& e) {
      std::cerr << "Meshing failed: " << e.GetMessageString() << std::endl;
      throw;
    }

    DATA facesOut = Array();
    TopExp_Explorer aExpFace;

    for(aExpFace.Init(aShape, TopAbs_FACE); aExpFace.More(); aExpFace.Next()) {
      try {
        TopoDS_Face aFace = TopoDS::Face(aExpFace.Current());
        if (aFace.IsNull()) continue;

        TopLoc_Location aLocation;
        Handle(Poly_Triangulation) aTr = BRep_Tool::Triangulation(aFace, aLocation);

        if(aTr.IsNull() || aTr->NbTriangles() == 0 || aTr->NbNodes() == 0) continue;

        DATA faceOut = Object();
        Handle(Geom_Surface) aSurface = BRep_Tool::Surface(aFace);
        
        if (!aSurface.IsNull()) {
          if (aSurface->IsKind("Geom_Plane")) {
            faceOut["surface"] = surfaceWrite(Handle(Geom_Plane)::DownCast(aSurface));
          } else {
            faceOut["surface"] = {"TYPE", "UNKNOWN"};
          }
        }

        if (!INTERROGATE_STRUCT_ONLY) {
          DATA tessOut = Array();
          writeFaceTessellation(aTr, aLocation, aSurface, tessOut);
          faceOut["tess"] = tessOut;
        }

        faceOut["inverted"] = aFace.Orientation() == TopAbs_REVERSED;
        TopoDS_Face* persistFace = new TopoDS_Face(aFace);
        faceOut["ref"] = e0::io::getStableRefernce(aFace);
        faceOut["ptr"] = ((std::uintptr_t)persistFace);
        
        facesOut.append(faceOut);
      } catch (Standard_Failure const& e) {
        std::cerr << "Face processing failed: " << e.GetMessageString() << std::endl;
        continue;  // Skip problematic faces
      }
    }

    out["faces"] = facesOut;
  } catch (Standard_Failure const& e) {
    std::cerr << "Interrogation failed: " << e.GetMessageString() << std::endl;
    throw;
  }
  
  return out;
}

void UpdateTessellation(TopoDS_Shape& shape, double deflection) {
  // Clean existing triangulation
  BRepTools::Clean(shape);
  
  // Perform incremental meshing with optimized parameters
  BRepMesh_IncrementalMesh mesher(shape, deflection, 
    Standard_True,   // relative
    0.5,             // angular deflection
    Standard_False   // parallel
  );
}

} // namespace io
} // namespace e0

#endif // E0_IO_INTERROGATE_H
