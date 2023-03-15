﻿using System;
using System.Collections.Generic;
using System.Linq;
using Grasshopper.Kernel;
using Rhino.Geometry;

namespace SCAD
{
    public class Test : GH_Component
    {
        /// <summary>
        /// Initializes a new instance of the Test class.
        /// </summary>
        public Test()
          : base("Transfinite Surface", "Transfinite",
              "Description",
             "SCAD", "Transfinite")
        {
        }

        /// <summary>
        /// Registers all the input parameters for this component.
        /// </summary>
        protected override void RegisterInputParams(GH_Component.GH_InputParamManager pManager)
        {
            pManager.AddCurveParameter("Curves", "C", "The bounding curves", GH_ParamAccess.list);
            pManager.AddIntegerParameter("Parametrization", "P", "Parametrization method", GH_ParamAccess.item);
            pManager.AddIntegerParameter("Blending", "B", "Blending Method", GH_ParamAccess.item);
        }

        /// <summary>
        /// Registers all the output parameters for this component.
        /// </summary>
        protected override void RegisterOutputParams(GH_Component.GH_OutputParamManager pManager)
        {
            pManager.AddMeshParameter("Model", "M", "Mesh Model", GH_ParamAccess.item);
            pManager.AddCurveParameter("Domain C", "D", "Domain curves", GH_ParamAccess.list);
            pManager.AddPointParameter("Domain V", "Dv", "Domain vertices", GH_ParamAccess.list);

        }

        /// <summary>
        /// This is the method that actually does the work.
        /// </summary>
        /// <param name="DA">The DA object is used to retrieve from inputs and store in outputs.</param>
        protected override void SolveInstance(IGH_DataAccess DA)
        {
            List<Curve> curves_ = new List<Curve>();
            int p =0, b = 0;
            DA.GetDataList(0, curves_);
            DA.GetData(1, ref p);
            DA.GetData(2, ref b);

            var pm = (Parametrization_Method)p;
            var bm = (Blending_Method)b;
            int resolution = 30;
            Domain dm = new Domain(resolution);
            dm.SetSides(curves_);
            dm.Update();

            var mesh = dm.MeshTopology();
            var uvs = dm.Parameters();
            var msh = mesh.Getmesh;
            var pts = new List<Point3d>();
            var sP = new SurfacePatch(dm, pm, bm);

            for (int i = 0; i < uvs.Count; i++)
            {
                //pts.Add(new Point3d(uvs[i].X, uvs[i].Y, 0));
                Point3d katoout = sP.Kato_Suv(uvs[i].X, uvs[i].Y);
                pts.Add(katoout);
            }

            //var sd = dm.Bounds;
            msh.Vertices.AddVertices(pts);
            msh.VertexColors.SetColors(Enumerable.Repeat(System.Drawing.Color.Silver, pts.Count).ToArray());
            DA.SetData(0, mesh.Getmesh);
            DA.SetDataList(2, dm.Vertices);

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
            get { return new Guid("5587F090-35B5-45DE-8AE8-64D81549853E"); }
        }
    }
}