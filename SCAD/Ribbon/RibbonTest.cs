using System;
using System.Collections.Generic;

using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SCAD
{
    public class RibbonTest : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the RibbonTest class.
        /// </summary>
        public RibbonTest()
          : base("RibbonTest", "Nickname",
              "Description",
              "Category", "Subcategory")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Curves", "C", "The bounding curves", GH_ParamAccess.list);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Model", "M1", "Mesh Model", GH_ParamAccess.item);
            pManager.AddMeshParameter("Model", "M2", "Mesh Model", GH_ParamAccess.item);
            pManager.AddPointParameter("Ribbon center", "Rc", "Ribbon center", GH_ParamAccess.list);
            pManager.AddVectorParameter("Ribbon Vector", "Rv", "Ribbon vectors", GH_ParamAccess.list);
        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var surf = new SCAD.Surface<RibbonCompatible>(Domain_Method.Domain_Regular,Parametrization_Method.RadialDistanceFunction, Blending_Method.Special_Side_Blending);
            int resolution = 10;
            double scaling = 20.0;
            double ribbon_length = 0.25;

            List<Curve> curves_ = new List<Curve>();
            DA.GetDataList(0, curves_);

            surf.SetCurves(curves_);
            surf.SetupLoop();
            surf.Update();

            int n = surf.Domain().Size;
            //n = 1;
            int size = n * resolution * 2;
            List<Point3d> pv = new List<Point3d>(); pv.Capacity = size;
            #region Fence
            TriMesh fence_mesh = new TriMesh();
            for (int i = 0; i < n; ++i) {
                for (int j = 0; j < resolution; ++j)
                {
                    double u = (double)j / resolution;
                    Point3d p = surf.Ribbon(i).Curve().PointAt(u);
                    pv.Add(p);
                    p += surf.Ribbon(i).Normal(u) * scaling;
                    pv.Add(p);
                }
            }

            fence_mesh.setPoints(pv);
            int index = 0;
            while (index < size - 2)
            {
                fence_mesh.addTriangle(index, index + 1, index + 2);
                ++index;
                fence_mesh.addTriangle(index, index + 2, index + 1);
                ++index;
            }
            fence_mesh.addTriangle(index, index + 1, 0);
            fence_mesh.addTriangle(index + 1, 1, 0);
            #endregion

            #region Ribbons
            // Ribbons
            TriMesh ribbon_mesh = new TriMesh();
            pv.Clear(); pv.Capacity = (n * (resolution + 1) * 2);
            List<Vector3d> vecList = new List<Vector3d>();
            List<Point3d> ptList = new List<Point3d>();
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= resolution; ++j)
                {
                    double u = (double)j / resolution;
                    pv.Add(surf.Ribbon(i).Curve().PointAt(u));
                    pv.Add(surf.Ribbon(i).Eval(new Point2d(u, ribbon_length)));
                    var pt = surf.Ribbon(i).Curve().PointAt(u);
                    ptList.Add(pt);

                    var vec = surf.Ribbon(i).Eval(new Point2d(u, ribbon_length)) - pt;
                    //vec.Unitize();
                    vecList.Add(vec);
                }
            }
            ribbon_mesh.setPoints(pv);
            index = 0;
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j < resolution; ++j)
                {
                    ribbon_mesh.addTriangle(index, index + 1, index + 2);
                    ++index;
                    ribbon_mesh.addTriangle(index, index + 2, index + 1);
                    ++index;
                }
                index += 2;
            }
            //ribbon_mesh.writeOBJ("../../models/" + filename + "-ribbons.obj");
            #endregion
            DA.SetData(0, ribbon_mesh.Getmesh);
            DA.SetData(1, fence_mesh.Getmesh);    
            DA.SetDataList(2, ptList);     
            DA.SetDataList(3, vecList);
        }

        /// <summary>
        /// Provides an Icon for the component.
        /// </summary>
        protected override System.Drawing.Bitmap Icon
        {
            get
            {
                //You can add image files to your project resources and access them like this:
                // return Resources.IconForThisComponent;
                return null;
            }
        }

        /// <summary>
        /// Gets the unique ID for this component. Do not change this ID after release.
        /// </summary>
        public override Guid ComponentGuid
        {
            get { return new Guid("C3C1A3D2-3E68-4B8F-943B-B93AA02C5813"); }
        }
    }
}