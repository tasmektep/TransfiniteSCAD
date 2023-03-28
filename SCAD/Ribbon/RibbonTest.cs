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

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            var surf = new SCAD.Surface<DomainRegular, Parametrization, RibbonCompatible>();
            int resolution = 10;
            double scaling = 20.0;
            double ribbon_length = 0.25;
            TriMesh fence_mesh = new TriMesh();

            List<Curve> curves_ = new List<Curve>();
            DA.GetDataList(0, curves_);

            surf.setCurves(curves_);
            surf.setupLoop();
            surf.update();

            int n = surf.domain().Size;
            int size = n * resolution * 2;
            List<Point3d> pv = new List<Point3d>(); pv.Capacity = size;

            for(int i = 0; i < n; ++i) {
                for (int j = 0; j < resolution; ++j)
                {
                    double u = (double)j / resolution;
                    Point3d p = surf.ribbon(i).curve().PointAt(u);
                    pv.Add(p);
                    p += surf.ribbon(i).normal(u) * scaling;
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
            //fence_mesh.writeOBJ("../../models/" + filename + "-fence.obj");


            
            // Ribbons
            TriMesh ribbon_mesh = new TriMesh();
            pv.Clear(); pv.Capacity = (n * (resolution + 1) * 2);
            for (int i = 0; i < n; ++i)
            {
                for (int j = 0; j <= resolution; ++j)
                {
                    double u = (double)j / resolution;
                    pv.Add(surf.ribbon(i).curve().PointAt(u));
                    pv.Add(surf.ribbon(i).eval(new Point2d(u, ribbon_length)));
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
            DA.SetData(0, ribbon_mesh.Getmesh);
            DA.SetData(1, fence_mesh.Getmesh);
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